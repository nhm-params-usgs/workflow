#!/bin/bash -e

# Create temperatures for tmax_allsnow
SYR=2007
EYR=2019

#for yy in {${SYR}..${EYR}}; do
for yy in $(eval echo "{$SYR..$EYR}"); do
    echo ${yy}

    # Mask temperature when precipitation type is 5 (snow) or 6 (wet snow)
    ncap2 -O -v -S combined_temp_pt1.nco download_${yy}_fixed.nc tmp_1.nc
    ncap2 -O -v -S combined_temp_pt2.nco tmp_1.nc tmp_2.nc

    # Compute mean maximum monthly temperature
    cdo monmean tmp_2.nc tmp_monthly_${yy}.nc
    #cdo monmax tmp2.nc tmp_monthly_{yy}.nc

    rm tmp_1.nc tmp_2.nc
done

ncrcat -O tmp_monthly_*.nc tmax_all_monthly_${SYR}-${EYR}.nc

# Computer mean monthly temperature for entire period
cdo ymonmean tmax_all_monthly_${SYR}-${EYR}.nc tmax_all_ymon_${SYR}-${EYR}.nc

rm tmp_monthly_*.nc


