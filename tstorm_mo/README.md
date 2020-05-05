
# General information

The `tstorm_mo` parameter is a monthly indicator of the predominate storm type where a value of zero indicates frontal storms and a value of one indicates convective storms. For the NHM this parameter is derived from the Climate Forecast System Reanalysis model (CFSR) which provides a global reanalysis dataset with horizontal grid spacings as fine as 0.5° x 0.5° and 64 levels in the vertical at sub-daily and monthly intervals (Saha and others, 2010). The CFSR model output was obtained from the Research Data Archive hosted by the National Center for Atmospheric Research (NCAR; http://rda.ucar.edu). Model output from CFSR includes separate convective and frontal precipitation variables in addition to the total precipitation. Convective and frontal precipitation values from CFSR for 1979–2010 were used to derive the `tstorm_mo` parameter for the NHM domain. First mean monthly values of convective and frontal precipitation were computed for 1979-2010. Next a convective mask by month was created; where the convective precipitation divided by the frontal precipitation was greater than 0.9 the mask was set to one, indicating a predominance of convective precipitation events occurring on average for the month, otherwise it was set to zero, indicating a predominance of frontal precipitation events. The `tstorm_mo` parameter is then computed for each month based on the area unweighted maximum value of the CFSR mask. 


# Workflow
Original work was done on tibicus.

## Create the CFSR convective precipitation mask
```bash
./create_monthly_mean_CONUS.sh
ncap2 -v -S tstorms.nco CFSR_precip_months_1979-2010_LTAVG.nc -o crap1.nc
ncatted -a level,,d,, -a long_name,TSTORM_MO,m,c,"Primarily convective precip"  -a standard_name,TSTORM_MO,d,, -a product_description,TSTORM_MO,m,c,"Precipitation which is at least 90% convective (1=convective; 0=largescale)" -a units,TSTORM_MO,m,c,"none" crap1.nc crap2.nc
mv crap2.nc CFSR_tstorm_mo.nc
```

The CFSR product has longitude 0 to 360; convert it to -180 to 180.
```bash
ncks -O --msa -d lon,181.,360. -d lon,0.,180. CFSR_tstorm_mo.nc out1.nc
ncap2 -O -s 'where(lon > 180) lon=lon-360' out1.nc out1.nc
```

## Generate weights for the CFSR domain
**NOTE**: The actual weights are not used in generating the `tstorm_mo` parameter but the grid ids for each HRU are used.
This script produces a weight file named `cfsr_weights.csv`.
```bash
python CalcWeights_CFSR_multiproc.py
```

## Create the `tstorm_mo` parameter from the CFSR mask
```bash
python process_cfsr_tstorm_mo.py
```

# Scripts

## `create_monthly_mean_CONUS.sh`
```bash
#!/bin/bash -e

# Author: Parker Norton (pnorton@usgs.gov)
# Create date: 2015-02-27
# Description: Creates the monthly mean from monthly source dataset(s)

MODEL=CFSR
VAR=precip      # temperature, precip
timevar=time    # name of the time variable

infile=../CFSR_${VAR}_monthly_1979-01_2010-12.nc

start=1       # starting time index
end=383       # ending time index

# Get spatial subset of dataset
ncks -O -d time,${start},${end} ${infile} tmp.nc

# Create the longterm monthly mean
cdo ymonavg -shifttime,-1month tmp.nc ${MODEL}_${VAR}_months_1979-2010_LTAVG.nc
rm tmp.nc

echo "Done!"
```

## `tstorms.nco`
```bash
TSTORM_MO=RAIN;
TSTORM_MO.set_miss(9.96921e+36);

where((RAINC/RAINNC) > .9)
{
    TSTORM_MO = 1;
} elsewhere {
    TSTORM_MO = 0;
}
lat=lat;
lon=lon;
time=time;
```


