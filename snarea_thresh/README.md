
# Information
From Regan et al. (2018) appendix: "The snarea_thresh parameter is the maximum threshold of SWE below which the snow-covered-area curve is applied. Computer code was developed to determine, for each HRU, the median of the yearly maximum SWE based on the SWE from [SNODAS](https://nsidc.org/data/g02158)."

Snow-water equivalent (SWE) from SNODAS is provided in units of meters.
The parameter, `snarea_thresh`, is in units of inches. 

The unmasked version of SNODAS, which provides greater coverage than the masked version, is used to generate the `snarea_thresh` parameter. 

**NOTE**: snarea_thresh=0 can cause a divide-by-zero in PRMS


# Requirements
Climate Data Operators (https://code.mpimet.mpg.de/projects/cdo/)

NetCDF Operators (http://nco.sourceforge.net/)

# Workflow
## Create yearly maximum SWE
For each year computed the yearly maximum. In line below replace YYYY with the year to process.
```bash
cdo yearmax SWE_YYYY.nc SWE_max_YYYY.nc
```

## Merge yearly maximum SWE files into one
```bash
ncrcat SWE_max_*.nc SWE_max_2013-2019.nc
```

## Computed the 50th-percentile (median) of yearly maximum SWE
```bash
cdo timpctl,50 SWE_max_2013-2019.nc -timmin SWE_max_2013-2019.nc -timmax SWE_max_2013-2019.nc SWE_median_of_max_yearly.nc
```

## Convert SWE to inches
```bash
ncap2 -O -s SWE=SWE*39.3701; SWE_median_of_max_yearly.nc -o SWE_median_of_max_yearly_inches.nc
ncatted -a units,SWE,m,c,inch SWE_median_of_max_yearly_inches.nc
```

## Restrict NHM bounding box
```bash
ncks -d lat,24.8193,52.90481 -d lon,-124.89666,-66.1607 SWE_median_of_max_yearly_inches.nc -o SWE_median_of_max_yearly_inches_NHM.nc
```

## Reset any zero-value to one
```bash
ncap2 -v -S set_zero_to_one.nco SWE_median_of_max_yearly_inches_NHM.nc SWE_median_of_max_yearly_inches_NHM_adj.nc
```

### Contents of `set_zero_to_one.nc`
```
where (SWE < 1.0)
{
    SWE = 1.0;
}
time=time;
lat=lat;
lon=lon;
crs=crs;
```

## Generate a mask for SNODAS
Use a single timestep from the unmask SNODAS precip variable to create a mask fo r the entire domain. The precip variable is used because it does not have any _FillValues in the model domain. 

Using the unzipped precip raster for a single timestep, the initial mask was created with:
```bash
gdal_translate -of netCDF -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -a_nodata -9999 -a_ullr -130.516666666661 58.2333333333310 -62.2499999999975 24.0999999999990 zz_ssmv01025SlL00T0024TTNATS2013010105DP001.dat tmp1.nc
```

Rename the default `Band1` variable to `prcp`
```bash
ncrename -v Band1,prcp tmp1.nc
```

Use `make_mask.nco` to convert the prcp data to a 0/1 mask
```bash
ncap2 -v -S ../make_mask.nco tmp1.nc tmp3.nc
```

The `make_mask.nco` file contains the following:
```bash
LANDMASK=prcp;
LANDMASK.set_miss(-32767);
where (LANDMASK == -9999.)
{
    LANDMASK = 0;
} elsewhere {
    LANDMASK = 1;
}
lat=lat;
lon=lon;
crs=crs;
```

Create a version that is restricted roughly to the NHM domain. 
__NOTE__: If this is done then the source netcdf data will also need to be restricted to the same bounding box. It is not required to pre-restrict to the NHM domain and would actually be less complicated if you don't. This shown because is what was done when developing the workflow for SNODAS weight generation.
```bash
ncks -d lat,24.8193,52.90481 -d lon,-124.89666,-66.1607 tmp/tmp3.nc -o snodas_landmask_NHM.nc
```

Create a new version of `snodas_landmask_NHM.nc` that includes monotonic arrays with 1) C-order and 2) Fortran-order. Paths and filename are set in the python file below.
```python create_grid_patterns.py```


## Create geopackage of SNODAS domain masked to the NHM domain
This step is done manually in QGIS. Processing steps are shown below. One geopackage is created for each numbered grid order. For creating the SNODAS weights only `grid_horiz` is needed.

Create a geopackage of the full-extent SNODAS grid (not needed when clipping below):
```
processing.run("native:pixelstopolygons", {'INPUT_RASTER':'NETCDF:\"./datasets/SNODAS/unmasked/snodas_test_patterns_NHM.nc\":grid_horiz','RASTER_BAND':1,'FIELD_NAME':'grid_horiz','OUTPUT':'./datasets/SNODAS/unmasked/tmp/snodas_grid_horiz_v1.gpkg'})
```

Clip SNODAS grid to the NHM domain using SAGA script:
__NOTE__: Original `NHM_Outline_v11_1kmBuffer` shapefile used for clipping came from Mike Wieczorek
```
processing.run("saga:cliprasterwithpolygon", {'INPUT':'NETCDF:\"./datasets/SNODAS/unmasked/snodas_test_patterns_NHM.nc\":grid_horiz','POLYGONS':'./datasets/SNODAS/unmasked/GIS/NHM_Outline_v11_1kmBuffer_WGS84.gpkg','OUTPUT':'./datasets/SNODAS/unmasked/GIS/snodas_grid_poly_clip_v2.sdat'})
```

Create a geopackage from the sdat-format clipped raster:
```
processing.run("native:pixelstopolygons", {'INPUT_RASTER':'./datasets/SNODAS/unmasked/GIS/snodas_grid_poly_clip_v2.sdat','RASTER_BAND':1,'FIELD_NAME':'grid_horiz','OUTPUT':'./datasets/SNODAS/unmasked/GIS/snodas_grid_NHMpoly_v1.gpkg'})
```

## Generate weights for SNODAS grid
Run `CalcWeights_SNODAS.py`
__NOTE__: Edit the top of the main function to set paths/directories
On denali this process took about 20 hours. The output from this process is the file `snodas_weights.csv`. This file is used in the next step.

## Compute area-weighted values of SWE for snarea_thresh parameter
Run `process_SNODAS_snarea_thresh.py`
__NOTE__: Edit the top of the main function to set paths/directories
This creates a paramDb-formatted file named `snarea_thresh.csv`




