# ERA5
This workflow was developed by Parker Norton (pnorton@usgs.gov).
## ERA5 information
The ERA5 product that was used for generating these parameters is `ERA5 hourly data on single levels from 1979 to present`. The webpage for this product is at: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
The data download page is at: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
The ERA5 variables of interest for these parameters are `2m_temperature`, `land_sea_mask`, and `precipitation_type`. The land_sea_mask is part of the workflow repository and does not need to be downloaded. 

### ERA5 precipitation_type (ptype)
Value | Description 
:---: | ----------- 
0 | No precipitation
1 | Rain
3 | Freezing rain
5 | Snow
6 | Wet snow
7 | Rain/Snow mix
8 | Ice pellets

## General workflow
### Requirements<a name='reqs'></a>
To complete the workflow described below you will need to install (or have available) the following: 
Climate Data Operators (https://code.mpimet.mpg.de/projects/cdo/)
NetCDF Operators (http://nco.sourceforge.net/)
**Python Libraries**
cdsapi (https://cds.climate.copernicus.eu/api-how-to)
csv
geopandas
multiprocessing
numpy
pandas
shapely
xarray

The Climate Data Operators and NetCDF Operators have further library requirements that are discussed on their websites. The workflow below used nco_4.9.2, cdo_1.9.8, and Python 3.7.


### Download the ERA5 product
Use `request.py` to download the ERA5 product. Edit this script to select the date ranges, product type, variables, and other information as needed. Each year of data is saved to a netcdf file named `download_YYYY.nc` where `YYYY` is the 4-digit year. Note that the script requires a library named cdsapi (see [Requirements](#reqs) for further information). 


### Unpack the data
Run the unpacking script. The unpack was moved to a separate step instead of letting ncap2 auto-unpack because the ncap2 program seems to have a problem with auto-unpack when creating more than one new variable in the output.
```bash
unpack_era5.sh
```

### Convert ptype to NC_SHORT
Unpacking results in a double ptype variable which is total overkill and wrong. I'm not sure why they decided to pack a variable that only has 7 possible values. Run the script below to fix it.
```bash
fix_ptype_era5.sh
```

### Create the tmax_allsnow and tmax_allrain variables
Run the script to process the data. The final product is a netcdf file that contains mean monthly `tmax_allsnow` and `tmax_allrain` for the period of record. Two netcdf files are created: 1) `tmax_all_monthly_BYYY-EYYY.nc` and 2) `tmax_all_ymon_BYYY-EYYY.nc` where `BYYY` and `EYYY` are the starting and ending years, respectively. Only the second file, which contains the mean monthly values for the period, is needed for the next step. The first file is kept in case it is useful to have mean monthly values by year. 
```bash
process_era5_combined.sh
```

### Generate the weights
To create the area-weighted tmax_allsnow and tmax_allrain_offset a set of weights needs to be created which maps the ERA5 grid cells to the model Hydrologic Response Units (HRUs). The land_sea_mask provides the grid cell information and a geospatial fabric (either a shapefile or geodatabase) provides the HRU information. The following script is used to generate the weights file. Variables that may need to be modified are in the top of the main function of the script. This script makes use of multiple cores in a computer. The number of cores used is hardcoded in the `mp_process_weights()` function; if needed it can be modified. This script is based off of code from the [onhm-fetcher-parser](https://github.com/nhm-usgs/onhm-fetcher-parser).
```bash
CalcERA5Weights_multiproc.py
```

### Generate the parameters
Once the weights file has been created the area-weighted averages by HRU can be computed for `tmax_allsnow` and `tmax_allrain_offset`. Variables that may need to be modified are in the top of the main function of the script. This script produces two output files: 1) `tmax_allsnow.csv` and 2) `tmax_allrain_offset.csv`. These can be dropped into the parameter database directory that was used to drive averaging. This script is based off of code from the [onhm-fetcher-parser](https://github.com/nhm-usgs/onhm-fetcher-parser).
```bash
process_era5_tmax_allstar.py
```

