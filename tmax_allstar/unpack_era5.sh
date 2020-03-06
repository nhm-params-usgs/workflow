#!/bin/bash -e

# Create temperatures for tmax_allsnow
SYR=2007
EYR=2019

#for yy in {${SYR}..${EYR}}; do
for yy in $(eval echo "{$SYR..$EYR}"); do
    echo ${yy}

    ncpdq -O -U download_${yy}.nc download_${yy}_upk.nc
done

