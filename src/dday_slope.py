# -*- coding: utf-8 -*-
"""
Created on Thu May  7 08:53:23 2020

@author: markstro
"""

import numpy as np
import pandas as pd
import geopandas as gpd
#import math
import sys
from datetime import datetime
from scipy.interpolate import interp1d

# File written by dday_slope-radadj.ipynb
#radadj_fn = 'c:/Users/markstro/work1.1/soltab/dday_slope_radadj_GF_v1.1.csv'
#gdb_path = 'c:/Users/markstro/work1.1/GIS/GFv1.1_v2f.gdb'

# Lots of stuff needs to be passed in.
#tmax_fn = 'c:/Users/markstro/work1.1/cbh_gridmet/tmax.cbh'
#prcp_fn = 'c:/Users/markstro/work1.1/cbh_gridmet/prcp.cbh'
#soltab_solt_fn = 'c:/Users/markstro/work1.1/soltab/soltab_solt_GF_v1.1.csv'
#hru_slope_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/hru_slope.csv'
#dday_slope_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/dday_slope_new2.csv'
#tmax_allrain_offset_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/tmax_allrain_offset.csv'
#tmax_allsnow_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/tmax_allsnow.csv'
#dday_intcp_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/dday_intcp.csv'
#nrel_shapefile_fn = "c:/Users/markstro/work1.1/GIS/nhm_shapefiles/nhm_hru_gf1_1_nrel_solrad.shp"


# Read the cbh file
def read_cbh(cbh_fn):
    foo = np.loadtxt(cbh_fn, skiprows=3,)
    dates = foo[:,0:6].astype('int')
    vals = foo[:,6:foo.shape[1]]
    return dates, vals

def save_values(ds, fn):
    f = open(fn, "w")
    f.write("$id,dday_slope\n")

    ii = 1
    for ih in range(ds.shape[0]):
        for im in range(ds.shape[1]):
            f.write(str(ii) + "," + str(ds[ih,im]) + "\n")
            ii, ds[ih,im]
            ii += 1
    f.close()
    
    
#def compute_monthly_dday_slope(dday_slope_daily, dates):
#    nhru = dday_slope_daily.shape[0]
#    ndays = dday_slope_daily.shape[1]
#    nmonths = 12
#
#    dday_slope = np.zeros(nhru*nmonths, dtype=np.float64)
#    dday_slope.shape = (nhru, nmonths)
#    count = np.zeros(nmonths, dtype=np.int)
#    
#    for iday in range(ndays):
#        yr = dates[iday,0]
#        mo = dates[iday,1] 
#        day = dates[iday,2]
#        date = datetime(yr,mo,day)
#        jday = date.timetuple().tm_yday
#        imon = date.month - 1
#
#        for ihru in range(nhru):
#            dday_slope[ihru,imon] += dday_slope_daily[ihru,iday]
#            count[imon] += 1
#            
#    for ihru in range(nhru):
#        for imon in range(0,nmonths):
#            dday_slope[ihru,imon] = dday_slope[ihru,imon] / count[imon]
#            
#    return dday_slope


# This computes the dday slope 
# dday_slope(ihru,imon) = (dday(iday,ihru) - dday_intcp(ihru,imon) - 1.0) / tmax_hru(iday,ihru)
# dday = Dday_slope(j, Nowmonth)*Tmax_hru(j) + Dday_intcp(j, Nowmonth) + 1.0
def compute_dday_slope(dday, dday_intcp, tmax):
    
#(114958, 12)
#(12, 114958)
#(12, 114958)
    foo = dday.transpose()

#    nmonths = tmax.shape[0]
#    nhru = tmax.shape[1]
    
#    dday_slope_monthly = np.zeros(nhru*nmonths, dtype=np.float64)
#    dday_slope_monthly.shape = (nhru,nmonths)
    
#    for imon in range(nmonths):
#        d = dates[iday]
#        yr = d[0]
#        mo = d[1]
#        day = d[2]
#        date = datetime(yr,mo,day)
#        jday = date.dayofyear
#        imon = date.month - 1
        
#        for ihru in range(nhru):
#            dday_slope_monthly[ihru,imon] = (dday[imon,ihru] - dday_intcp[imon,ihru] - 1.0) \
#                                        / tmax[imon,ihru]       

    dday_slope_monthly = (foo - dday_intcp - 1.0) / tmax
    
    return dday_slope_monthly


def compute_radadj(hrus_dni, hru_slope_df, nhru_v11_vals, soltab_df_vals):
# This code copied from jupyter notebook dday_slope-radadj
# Compute radadj for all days-of-the-year and for all HRUs


# For whatever reason, the zonal mean process from QGIS left some NaNs when filling in the Direct Normal Irradiance (dni)
# monthly values in the shapefile (see dni* columns above). Those need to be filled in with real values.
# Without going back to the GIS (which I already ran with the results that are shown above), I am using np.interpolate across
# each column. This is a hack in the sense that the adjacent row do not necessarily mean that the HRUs are adjacent, and it
# is uncertain exactly what is being interpolated, but it is filling in the nan values with real values and allows
# me to move on.
#    count = 0
#    for ii in range(hrus_dni.shape[0]):
#        for jj in range(8,20):
#            if np.isnan(hrus_dni.iloc[ii,jj]):
#                count += 1
#    print(count)
          
    
#    for ii in range(hrus_dni.shape[0]):
#        for jj in range(8,20):
#            hrus_dni.iloc[:,jj] = hrus_dni.iloc[:,jj].interpolate()
    
#    print("subsequent coordinates with nan value")
#    count = 0
#    for ii in range(hrus_dni.shape[0]):
#        for jj in range(8,20):
#            if np.isnan(hrus_dni.iloc[ii,jj]):
#                print(" (", ii,jj, ")", end = '')
#                count += 1
                
# Join the HRU slopes from the PRMS parmaeter file (paramdb) to the features from the shapefile.
    hrus_1 = hrus_dni.set_index('nhru_v11').join(hru_slope_df.set_index('$id'))

# Take the cosine of the hru slopes. The hru_slope must be converted from rise/run to radians with arctan first.
    hrus_1['hru_cossl'] = np.cos(np.arctan(hrus_1['hru_slope']))
    hru_cossl_vals = hrus_1['hru_cossl'].values
    
# Extract the NREL monthly solrad targets from the dataframe into a numpy 2D array for the radadj calculation below.
    solrad_targets_vals= hrus_1[["dni_jan_me", "dni_feb_me", "dni_mar_me", "dni_apr_me", "dni_may_me", "dni_jun_me",
                    "dni_jul_me", "dni_aug_me", "dni_sep_me", "dni_oct_me", "dni_nov_me", "dni_dec_me"]].values

# The NREL solrad values are in units of kWh/m2/Day. PRMS uses Langleys per day. The conversion factor comes from
# https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/GEM/SolarRadConversion.pdf

#                    W-sec    1 KW     1 hour               KW-hours
# 1 Langley = 41868 -------  ------   -------  =   0.6978  ---------
#                     m2     1000 W   60 sec                  m2

    solrad_targets_vals_langleys = solrad_targets_vals / 0.6978

    # create an array (len = 366, number of days in the year) that for any
    # jday, it gives the month index.
#    jan = [0] * 31
#    feb = [1] * 29 # assume leap year to get full table
#    mar = [2] * 31
#    apr = [3] * 30
#    may = [4] * 31
#    jun = [5] * 30
#    jul = [6] * 31
#    aug = [7] * 31
#    sep = [8] * 30
#    octo = [9] * 31
#    nov = [10] * 30
#    dec = [11] * 31
#    month_of_jday = jan + feb + mar + apr + may + jun + jul + aug + sep + octo + nov + dec

    radadj = np.zeros(soltab_df_vals.shape)
    
    print(radadj.shape)
    print(solrad_targets_vals_langleys.shape)
    print(hru_cossl_vals.shape)
    print(soltab_df_vals.shape)
#    (366, 114958)
#    (114958, 12)
#    (114958,)
#    (114958, 12)
    nhru = soltab_df_vals.shape[0]
    nmonths = soltab_df_vals.shape[1]
    
    min_count = 0
    max_count = 0
    for imon in range(nmonths):
        for ihru in range(nhru):
            kk = nhru_v11_vals[ihru] - 1
            try:
#                radadj[jday,ihru] = solrad_targets_vals_langleys[ihru,imon] * hru_cossl_vals[ihru] / soltab_df_vals[jday,kk]
                radadj[kk,imon] = solrad_targets_vals_langleys[ihru,imon] * hru_cossl_vals[kk] / soltab_df_vals[kk,imon]

            except:
                print(ihru, imon, solrad_targets_vals_langleys[ihru,imon], hru_cossl_vals[ihru], soltab_df_vals[kk,imon])
                
    for imon in range(nmonths):
        for ihru in range(nhru):
            if radadj[ihru,imon] < 0.05:
                radadj[ihru,imon] = 0.05
                min_count += 1
                
            if radadj[ihru,imon] > 0.95:
                radadj[ihru,imon] = 0.05
                max_count += 1
            
#    print ((jday * nhru), min_count, max_count)

    return radadj


# compute the dday values for each time step (day) for each HRU
def compute_dday(radadj_solf):
    # Here's the cubic spline interpolation between solf and dday
    solf_pts = np.array([.20, .35, .45, .51, .56, .59, .62, .64, .655, .67, \
                         .682, .69, .70, .71, .715, .72, .722, .724, .726, \
                             .728, .73, .734, .738, .742, .746, .75])
    dday_pts = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, \
                         11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, \
                             19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0])
    f2 = interp1d(solf_pts, dday_pts, kind='cubic')
    
    nmonths = radadj_solf.shape[0]
    nhru = radadj_solf.shape[1]
    
    dday = np.zeros(nmonths*nhru, dtype=np.float64)
    dday.shape = (nmonths, nhru)
    
    for imon in range(nmonths):
        for ihru in range(nhru):
            if radadj_solf[imon, ihru] < 0.2:
                dday[imon,ihru] = 1.0
            elif radadj_solf[imon, ihru] > 0.75:
                dday[imon,ihru] = 26.0
            else:
                dday[imon,ihru] = f2(radadj_solf[imon,ihru])
    return dday
                



def compute_radadj_solf(radadj, pptadj):
#    nmonths = pptadj.shape[0]
#    nhru = pptadj.shape[1]
    
#    radadj_solf = np.zeros(nmonths * nhru)
#    radadj_solf.shape = (nmonths,nhru)
    
#    for imon in range(nmonths):
#        yr = dates[iday,0]
#        mo = dates[iday,1]
#        day = dates[iday,2]
#        date = datetime(yr,mo,day)
#        jday = date.timetuple().tm_yday
        
#        for ihru in range(nhru):
##            print(iday,ihru,jday-1)
#            try:
#                radadj_solf[imon,ihru] = radadj[imon,ihru] * pptadj[imon,ihru]
#            except Exception as inst:
#                print(type(inst))
#                print(inst)
#                print(iday,ihru,jday-1,radadj.shape, pptadj.shape)

    radadj_solf = radadj * pptadj
    return radadj_solf

def compute_pptadj(prcp_in_monthly, prcp_days_monthly, ppt_rad_adj, tmax_f_monthly, tmax_index, radj_sppt, tmax_allrain, radj_wppt, radadj_intcp, radadj_slope):
    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    
    nmonths = tmax_f_monthly.shape[0]
    nhru = tmax_f_monthly.shape[1]
    
    pptadj = np.zeros(nmonths*nhru, dtype=np.float64)
    pptadj.shape = (nmonths,nhru)
    
    for imon in range(nmonths):
        summer_flag = True
        if imon < 3 or imon > 8:
            summer_flag = False
        
        for ihru in range(nhru):
            if tmax_f_monthly[imon,ihru] < tmax_index:
                pptadj[imon,ihru] = radj_sppt
                if tmax_f_monthly[imon,ihru] >= tmax_allrain[imon,ihru]:
                    if not summer_flag:
                        pptadj[imon,ihru] = radj_wppt
                else:
                    pptadj[imon,ihru] = radj_wppt
            else:
                pptadj[imon,ihru] = radadj_intcp + radadj_slope * (tmax_f_monthly[imon,ihru] - tmax_index)

            # The monthly calculation doesn't consider whether the daily values indicate a rain day or not, so
            # determine the decimal fraction of rain days in the month and divide by it to bring the pptadj closer
            # to 1.0.
            if prcp_days_monthly[imon,ihru] == 0:
                prcp_days_fac = 1.0
            else:
                prcp_days_fac = days_in_month[imon] / prcp_days_monthly[imon,ihru]
                
            pptadj[imon,ihru] = pptadj[imon,ihru] * prcp_days_fac
            
            if pptadj[imon,ihru] > 1.0:
                pptadj[imon,ihru] = 1.0
    return pptadj


def read_temp(cbh_path):
    
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


def read_prcp(cbh_path, ppt_rad_adj):
    
    dates, vals = read_cbh(cbh_path)
    nhru = vals.shape[1]
    ndays = vals.shape[0]

    mon_cnt = 12

    vals_monthly = np.zeros(mon_cnt * nhru)
    vals_monthly.shape = (mon_cnt, nhru)
    monthly_count = np.zeros(mon_cnt * nhru, dtype=np.int)
    monthly_count.shape = (mon_cnt, nhru)
    prcp_days_monthly = np.zeros(mon_cnt * nhru, dtype=np.int)
    prcp_days_monthly.shape = (mon_cnt, nhru)
    
    for iday in range(ndays):
        date = datetime(dates[iday,0],dates[iday,1],dates[iday,2])

        imon = date.month - 1
        for ihru in range(nhru):
            vals_monthly[imon,ihru] += vals[iday,ihru]
            if vals[iday,ihru] >= ppt_rad_adj:
                prcp_days_monthly[imon,ihru] += 1
            monthly_count[imon, ihru] +=1    
            
    vals_monthly = vals_monthly / monthly_count
#    prcp_days_monthly = prcp_days_monthly / monthly_count
    
    del(vals)
    
    return vals_monthly, prcp_days_monthly                  

                
def main(radadj_fn, gdb_path, tmax_fn, prcp_fn, soltab_solt_fn, hru_slope_fn, dday_slope_fn, tmax_allrain_offset_fn, tmax_allsnow_fn, dday_intcp_fn, nrel_shapefile_fn):

# Read the necessary parameters
    hru_slope_df = pd.read_csv(hru_slope_fn)
    nhru = hru_slope_df.shape[0]
    nmonths = 12
    
    tmax_allrain_offset_f = pd.read_csv(tmax_allrain_offset_fn)
    tmax_allrain_offset = tmax_allrain_offset_f["tmax_allrain_offset"].values
    
    tmax_allsnow_f = pd.read_csv(tmax_allsnow_fn)
    tmax_allsnow = tmax_allsnow_f["tmax_allsnow"].values
    tmax_allrain = tmax_allsnow + tmax_allrain_offset
    
    #    dday_intcp[imon,ihru]
    dday_intcp_f = pd.read_csv(dday_intcp_fn)
    dday_intcp = dday_intcp_f["dday_intcp"].values
    dday_intcp.shape = (nmonths,nhru)
    
    print(tmax_allrain)
    print(tmax_allrain.shape)
    
    tmax_allrain.shape = (nmonths, nhru)
    
    print(tmax_allrain.shape)
    
    # These 2D parameters are all set to constants
    ppt_rad_adj = 0.02
    tmax_index = 50.0
    radj_sppt = 0.44
    radj_wppt = 0.55
    radadj_intcp = 1.0
    radadj_slope = 0.02    

# Temperature is in C
#dates, tmax = read_cbh(tmax_fn)
    tmax_c_monthly = read_temp(tmax_fn)
    tmax_f_monthly = tmax_c_monthly * 1.8 +32

# prcp is in mm/day
#dates, prcp = read_cbh(prcp_fn)
    ppt_rad_adj_mm = ppt_rad_adj * 25.4
    prcp_mm_monthly, prcp_days_monthly = read_prcp(prcp_fn, ppt_rad_adj_mm)
    prcp_in_monthly = prcp_mm_monthly / 25.4

# set dimension sizes
    nhru = prcp_in_monthly.shape[1]
    nmonths = prcp_in_monthly.shape[0]

# read the soltab values from the csv file produced by solar_table.ipynb
# It's OK to read these values from the table because they won't change
# regardless of changes to the weather values. The file specified has the
# "clear sky" solar radiation values for each HRU of the geospatial fabric
# version 1.1 for each day-of-the-year.
    soltab_df = pd.read_csv(soltab_solt_fn, header=None)
    soltab_df_vals = soltab_df.values
        
# replace negative and super small soltab values.
    soltab_df_vals[soltab_df_vals < 10.0] = 10.0


    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    
    soltab_monthly = np.zeros(nhru * nmonths)
    soltab_count = np.zeros(nhru * nmonths, dtype=np.int)
    soltab_monthly.shape = (nhru, nmonths)
    soltab_count.shape = (nhru, nmonths)

    ll = 0
    for imon in range(nmonths):
        for jj in range(days_in_month[imon]):
            for ihru in range(nhru):
                soltab_monthly[ihru,imon] += soltab_df.iloc[ll,ihru]
                soltab_count[ihru,imon] +=1
        ll += 1
        
    soltab_monthly = soltab_monthly / soltab_count

# Read in the monthly short wave data from NREL. It's in the shapefile that
# was created by running zonal means (to map to GF v1.1 HRUs) on the monthly
# DNI geotifs downloaded from NREL.
    hrus_dni = gpd.read_file(nrel_shapefile_fn)
    
    for jj in range(8,20):
        hrus_dni.iloc[:,jj] = hrus_dni.iloc[:,jj].interpolate()
    
    count = 0
    for ii in range(hrus_dni.shape[0]):
        for jj in range(8,20):
            if np.isnan(hrus_dni.iloc[ii,jj]):
    #            print(" (", ii,jj, ")", end = '')
                count += 1

    hrus_dni_vals = hrus_dni.iloc[:,8:20].values
    hru_dni_nhru_v11_vals = hrus_dni['nhru_v11'].values
    
# Get the order of the nhru_v11 IDs for mapping the soltab values.
    nhru_v11_vals = hrus_dni["nhru_v11"]
    
# Start the code to compute pptadj
# The code from dday_slope-pptadj.ipynb is copied in here
# I copied in the code, as opposed to reading these already computed values
# from the csv file because they will change depending on the source, length
# of record, etc. of the values in the CBH files. I.e., these need to be
# recomputed if the weather data changes!
    pptadj = compute_pptadj(prcp_in_monthly, prcp_days_monthly, ppt_rad_adj, tmax_f_monthly, \
                        tmax_index, radj_sppt, tmax_allrain, \
                        radj_wppt, radadj_intcp, radadj_slope)

   
# Compute the values for radadj
    radadj = compute_radadj(hrus_dni, hru_slope_df, nhru_v11_vals, soltab_monthly)    
    
#radadj_solf
    radadj_solf = compute_radadj_solf(radadj, pptadj.transpose())
    
# compute the dday values for each time step (month) for each HRU
    dday = compute_dday(radadj_solf)
    
# This computes the dday slope 
# dday_slope(ihru,imon) = (dday(iday,ihru) - dday_intcp(ihru,imon) - 1.0) / tmax_hru(iday,ihru)
# dday = Dday_slope(j, Nowmonth)*Tmax_hru(j) + Dday_intcp(j, Nowmonth) + 1.0
    dday_slope_monthly = compute_dday_slope(dday, dday_intcp, tmax_f_monthly)
    
# save values to file
    save_values(dday_slope_monthly, dday_slope_fn)


if __name__ == "__main__":
    # File written by dday_slope-radadj.ipynb
#radadj_fn = 'c:/Users/markstro/work1.1/soltab/dday_slope_radadj_GF_v1.1.csv'
#gdb_path = 'c:/Users/markstro/work1.1/GIS/GFv1.1_v2f.gdb'

# Lots of stuff needs to be passed in.
#tmax_fn = 'c:/Users/markstro/work1.1/cbh_gridmet/tmax.cbh'
#prcp_fn = 'c:/Users/markstro/work1.1/cbh_gridmet/prcp.cbh'
#soltab_solt_fn = 'c:/Users/markstro/work1.1/soltab/soltab_solt_GF_v1.1.csv'
#hru_slope_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/hru_slope.csv'
#dday_slope_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/dday_slope_new2.csv'
#tmax_allrain_offset_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/tmax_allrain_offset.csv'
#tmax_allsnow_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/tmax_allsnow.csv'
#dday_intcp_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/dday_intcp.csv'
#nrel_shapefile_fn = "c:/Users/markstro/work1.1/GIS/nhm_shapefiles/nhm_hru_gf1_1_nrel_solrad.shp"

    
    if (len(sys.argv) == 12):     #test whether any arguments have been passed in
        radadj_fn = sys.argv[1]
        gdb_path = sys.argv[2]
        tmax_fn = sys.argv[3]
        prcp_fn = sys.argv[4]
        soltab_solt_fn = sys.argv[5]
        hru_slope_fn = sys.argv[6]
        dday_slope_fn = sys.argv[7]
        tmax_allrain_offset_fn = sys.argv[8]
        tmax_allsnow_fn = sys.argv[9]
        dday_intcp_fn = sys.argv[10]
        nrel_shapefile_fn = sys.argv[11]
    else:
        print("No name passed in")
        
        
    main(radadj_fn, gdb_path, tmax_fn, prcp_fn, soltab_solt_fn, hru_slope_fn, \
         dday_slope_fn, tmax_allrain_offset_fn, tmax_allsnow_fn, \
             dday_intcp_fn, nrel_shapefile_fn)
