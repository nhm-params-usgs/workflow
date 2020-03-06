#!/bin/bash -e

# Create temperatures for tmax_allsnow
SYR=2007
EYR=2019

#for yy in {${SYR}..${EYR}}; do
for yy in $(eval echo "{$SYR..$EYR}"); do
    echo ${yy}
    # Convert ptype to NC_SHORT
    ncap2 -O -s "ptype=ptype.convert(NC_SHORT)" download_${yy}_upk.nc download_${yy}_fixed.nc
done
