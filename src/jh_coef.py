# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:38:24 2020

@author: markstro
"""

import numpy as np
import pandas as pd
#import fiona 
#import math
#import xarray as xa
#import matplotlib.pyplot as plt
#%matplotlib inline
#from matplotlib import cm
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from datetime import datetime
#import random
#from scipy.interpolate import interp1d
import geopandas as gpd
import sys

#gdb_path = 'c:/Users/markstro/work1.1/GIS/GFv1.1_v2f.gdb'
#markstro_param_path = "c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/"
#path = 'c:/Users/markstro/work1.1/weather/'
#farns_v11_path = 'c:/Users/markstro/work1.1/GIS/farns_gf_1_1_mm_per_month.shp'
#jh_coef_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/jh_coef_new.csv'
#jh_coef_hru_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/jh_coef_hru.csv'
#nrel_shapefile_fn = "c:/Users/markstro/work1.1/GIS/nhm_shapefiles/nhm_hru_gf1_1_nrel_solrad.shp"
#tmax_path = 'c:/Users/markstro/work1.1/cbh_gridmet/tmax.cbh'
#tmin_path = 'c:/Users/markstro/work1.1/cbh_gridmet/tmin.cbh'


def write_jh_coef(fn, vals):
    f = open(fn, "w")
    f.write("$id,jh_coef\n")
    kk = 1
    for ii in range(vals.shape[1]):
        for jj in range(vals.shape[0]):
            f.write(str(kk) + "," + str(vals[jj,ii]) + '\n')
            kk = kk + 1
            
            
# Read the cbh file
def read_cbh(cbh_fn):
#    foo = np.loadtxt(cbh_fn, delimiter=' ', skiprows=3,)
    foo = np.loadtxt(cbh_fn, skiprows=3,)

    dates = foo[:,0:6].astype('int')
    vals = foo[:,6:foo.shape[1]]
    
    return dates, vals

            
def compute_jh_coef(tavg_f_monthly, tavg_c_monthly, v_id, farns_vals, jh_coef_hru_vals, solrad_targets_vals_langleys):
# DO j = 1, Active_hrus
#          elh = (597.3-(0.5653*Tavgc(i)))*2.54
#          Potet(i) = Jh_coef(i, Nowmonth)*(Tavgf(i)-Jh_coef_hru(i))*Swrad(i)/elh
#          IF ( Potet(i)<0.0) Potet(i) = 0.0
#
# so
#          elh = (597.3-(0.5653*Tavgc(iday,ihru)))*2.54
#          Jh_coef(ihru, imon) = farns(ihru, imon) / ((Tavgf(iday,ihru)-Jh_coef_hru(ihru)) / Swrad(imon,ihru) * elh)

    nmo = tavg_f_monthly.shape[0]
    nhru = tavg_f_monthly.shape[1]
    
    jh_coef = np.zeros(nhru * nmo)
    jh_coef.shape = (nhru, nmo)
    
    for imon in range(nmo):
        for ihru in range(nhru):
            v11_id = v_id[ihru] - 1
            if tavg_c_monthly[imon, v11_id] == 0.0:
                tavg_c_monthly[imon, v11_id] = 1.0
            elh = (597.3 - (0.5653 * tavg_c_monthly[imon, v11_id])) * 2.54
            jh_coef[v11_id, imon] = farns_vals[ihru, imon] / (tavg_f_monthly[imon,v11_id] - jh_coef_hru_vals[v11_id]) / solrad_targets_vals_langleys[ihru,imon] * elh

    return jh_coef


def read_temperature(cbh_path):
    
    dates, vals = read_cbh(cbh_path)
    nhru = vals.shape[1]
    ndays = vals.shape[0]

    mon_cnt = 12

    vals_monthly = np.zeros(mon_cnt * nhru)
    vals_monthly.shape = (mon_cnt, nhru)
    monthly_count = np.zeros(mon_cnt * nhru, dtype=np.int)
    monthly_count.shape = (mon_cnt, nhru)
    
    for iday in range(ndays):
        date = datetime(dates[iday,0],dates[iday,1],dates[iday,2])

        imon = date.month - 1
        for ihru in range(nhru):
            vals_monthly[imon,ihru] += vals[iday,ihru]
            monthly_count[imon, ihru] +=1    
            
    vals_monthly = vals_monthly / monthly_count
    del(vals)
    return vals_monthly


def read_farns(shape_fn):
    # Read the v11 Farnsworth shapefile
    # The PET values (cols jan - dec) are in units of mm / month
    farns_v11 = gpd.read_file(farns_v11_path)
    
    # Extract the PET values and convert from mm/month to inches/day
    farns_vals = farns_v11.iloc[:,4:16].values
    farns_vals[farns_vals < 1.0] = 1.0
    
    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    
    # Extract the PET values and convert from mm/month to inches/day
    nhru = farns_vals.shape[0]
    nmonth = farns_vals.shape[1]
    
    for imon in range(nmonth):
        farns_vals[:,imon] = farns_vals[:,imon] / 25.4 / days_in_month[imon]
    #    farns_vals[:,imon] = farns_vals[:,imon] / 25.4
    
    # Extract the Farnsworth v11 HRU order
    farns_nhru_v11 = farns_v11['nhru_v11'].values

    return farns_vals, farns_nhru_v11


def compute_solrad_targets_vals_langleys(nrel_shapefile_fn):
    # Read in the monthly short wave data from NREL. It's in the shapefile that was created by running zonal means
    # on the monthly DNI geotifs downloaded from NREL.
    nrel_solrad = gpd.read_file(nrel_shapefile_fn)
    nhru_v11_id = nrel_solrad["nhru_v11"].values
    
    # Look for nans in the NREL values and subsitute with interpolated value
    for jj in range(8,20):
        nrel_solrad.iloc[:,jj] = nrel_solrad.iloc[:,jj].interpolate()
    
    solrad_targets_vals = nrel_solrad.iloc[:,8:20].values
    
# units of the NREL data are w/m2
# so 1 Langley/day = 0.484583 Watt/m2
# so 1 watt/m2 = 2.06363 Lang/day
    solrad_targets_vals_langleys = solrad_targets_vals * 2.06363
      
    return solrad_targets_vals_langleys, nhru_v11_id


def main(jh_coef_fn, nrel_shapefile_fn, farns_v11_path, jh_coef_hru_fn, \
         tmax_path, tmin_path):
    # Read the  temperature data from the CBHs.
    tmax_c_monthly = read_temperature(tmax_path)
    tmin_c_monthly = read_temperature(tmin_path)
    tavg_c_monthly = (tmax_c_monthly + tmin_c_monthly) / 2.0
    tavg_f_monthly = (tavg_c_monthly * 1.8) + 32
    
# Read the Farnsworth PET target values from the shapefile
    farns_vals, farns_nhru_v11 = read_farns(farns_v11_path)
    
# read the jh_coef_hru values from the csv file produced by jh_coef_hru.ipynb
    jh_coef_hru = pd.read_csv(jh_coef_hru_fn)
    jh_coef_hru_vals = jh_coef_hru['jh_coef_hru'].values
#    jh_coef_hru_idx = jh_coef_hru['$id'].values

# Read in the monthly short wave data from NREL. It's in the shapefile that was created by running zonal means
# on the monthly DNI geotifs downloaded from NREL.
    solrad_targets_vals_langleys, nhru_v11_id = compute_solrad_targets_vals_langleys(nrel_shapefile_fn)
    
    jh_coef = compute_jh_coef(tavg_f_monthly, tavg_c_monthly, nhru_v11_id, \
                              farns_vals, jh_coef_hru_vals, solrad_targets_vals_langleys)
    write_jh_coef(jh_coef_fn, jh_coef)


if __name__ == "__main__":
#    jh_coef_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/jh_coef_new.csv'
#    nrel_shapefile_fn = "c:/Users/markstro/work1.1/GIS/nhm_shapefiles/nhm_hru_gf1_1_nrel_solrad.shp"
#    farns_v11_path = 'c:/Users/markstro/work1.1/GIS/farns_gf_1_1_mm_per_month.shp'
#    jh_coef_hru_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/jh_coef_hru.csv'
#    tmax_path = 'c:/Users/markstro/work1.1/cbh_gridmet/tmax.cbh'
#    tmin_path = 'c:/Users/markstro/work1.1/cbh_gridmet/tmin.cbh'
    
    if (len(sys.argv) == 7):     #test whether any arguments have been passed in
        jh_coef_fn = sys.argv[1]
        nrel_shapefile_fn = sys.argv[2]
        farns_v11_path = sys.argv[3]
        jh_coef_hru_fn = sys.argv[4]
        tmax_path = sys.argv[5]
        tmin_path = sys.argv[6]
    else:
        print("No name passed in")
        
        
    main(jh_coef_fn, nrel_shapefile_fn, farns_v11_path, jh_coef_hru_fn, \
         tmax_path, tmin_path)
