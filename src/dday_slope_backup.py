# -*- coding: utf-8 -*-
"""
Created on Thu May  7 08:53:23 2020

@author: markstro
"""

import numpy as np
import pandas as pd
#import sys
from datetime import datetime
from scipy.interpolate import interp1d
import geopandas as gpd

# Set the file names. These need to be passed in.
tmax_fn = 'c:/Users/markstro/work1.1/cbh_gridmet/tmax.cbh'
prcp_fn = 'c:/Users/markstro/work1.1/cbh_gridmet/prcp.cbh'
param_path = "c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/"

# This file produced by jupyter notebook "solar_table". It produces the
# same values as the soltab.f90 PRMS module.
soltab_solt_fn = 'c:/Users/markstro/work1.1/soltab/soltab_solt_GF_v1.1.csv'

# This is the hru_slope value from the PRMS parameter file
hru_slope_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/hru_slope.csv'
dday_slope_fn = 'c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_master/dday_slope_new.csv'
# This is a shapefile that was created by running zonal means (to map to
# GF v1.1 HRUs) on the monthly DNI geotifs downloaded from NREL.

nrel_shapefile_fn = "c:/Users/markstro/work1.1/GIS/nhm_shapefiles/nhm_hru_gf1_1_nrel_solrad.shp"

# var_name = "tmin"
# frmt = "%.1f"


# Read the cbh file
def read_cbh(cbh_fn):
    foo = np.loadtxt(cbh_fn, skiprows=3,)
    dates = foo[:,0:6].astype('int')
    vals = foo[:,6:foo.shape[1]]
    return dates, vals


def main():

# This section computes the pptadj values.
    # Read the weather data to compute pptadj
    # Temperature is in C
    dates, tmax = read_cbh(tmax_fn)
    tmax_f = tmax * 1.8 +32
    
    # Destroy reference to free memory
    del tmax
    
    # prcp is in mm/day
    dates, prcp = read_cbh(prcp_fn)
    prcp_in = prcp / 25.4
    
    # Destroy reference to free memory
    del prcp
    
    # set dimension sizes
    ndays = prcp_in.shape[0]
    nhru = prcp_in.shape[1]
    nmonths = 12

# Read the necessary parameters to compute pptadj
    tmax_allrain_offset_f = pd.read_csv(param_path + "tmax_allrain_offset.csv")
    tmax_allrain_offset = tmax_allrain_offset_f["tmax_allrain_offset"].values
    
    tmax_allsnow_f = pd.read_csv(param_path + "tmax_allsnow.csv")
    tmax_allsnow = tmax_allsnow_f["tmax_allsnow"].values
    tmax_allrain = tmax_allsnow + tmax_allrain_offset
    
#    dday_intcp[imon,ihru]
    dday_intcp_f = pd.read_csv(param_path + "dday_intcp.csv")
    dday_intcp = dday_intcp_f["dday_intcp"].values
    
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

# Start the code to compute pptadj
# The code from dday_slope-pptadj.ipynb is copied in here
# I copied in the code, as opposed to reading these already computed values
# from the csv file because they will change depending on the source, length
# of record, etc. of the values in the CBH files. I.e., these need to be
# recomputed if the weather data changes!
    pptadj = np.zeros(ndays*nhru, dtype=np.float64)
    pptadj.shape = (ndays,nhru)
    
    for iday in range(ndays):
        yr = dates[iday,0]
        mo = dates[iday,1]
        day = dates[iday,2]
        date = datetime(yr,mo,day)
        jday = date.timetuple().tm_yday
        imon = date.month - 1
        summer_flag = True
        if jday < 79 or jday > 265:
            summer_flag = False
    #    print(iday, yr, mo, day, date, jday, imon, summer_flag)
        
        for ihru in range(nhru):
            pptadj[iday,ihru] = 1.0
    #        print(prcp_in[iday,ihru], ppt_rad_adj)
            if prcp_in[iday,ihru] > ppt_rad_adj:
    #            print(tmax[iday,ihru], tmax_index)
                if tmax_f[iday,ihru] < tmax_index:
                    pptadj[iday,ihru] = radj_sppt
                    if tmax_f[iday,ihru] >= tmax_allrain[imon,jday]:
                        if not summer_flag:
                            pptadj[iday,ihru] = radj_wppt
                    else:
                        pptadj[iday,ihru] = radj_wppt
                else:
                    pptadj[iday,ihru] = radadj_intcp + radadj_slope * (tmax_f[iday,ihru] - tmax_index)
    #                print(pptadj[iday,ihru])
                    if pptadj[iday,ihru] > 1.0:
                        pptadj[iday,ihru] = 1.0

# End of the calculation of pptadj

# Compute the values for radadj
# This code copied from jupyter notebook dday_slope-radadj
# Compute radadj for all days-of-the-year and for all HRUs

# read the soltab values from the csv file produced by solar_table.ipynb
# It's OK to read these values from the table because they won't change
# regardless of changes to the weather values. The file specified has the
# "clear sky" solar radiation values for each HRU of the geospatial fabric
# version 1.1 for each day-of-the-year.

    soltab_df = pd.read_csv(soltab_solt_fn, header=None)
    soltab_df_vals = soltab_df.values
    
# Read in the monthly short wave data from NREL. It's in the shapefile that
# was created by running zonal means (to map to GF v1.1 HRUs) on the monthly
# DNI geotifs downloaded from NREL.
    hrus_dni = gpd.read_file(nrel_shapefile_fn)

# Get the order of the nhru_v11 IDs for mapping the soltab values.
    nhru_v11_vals = hrus_dni["nhru_v11"]
# For whatever reason, the zonal mean process from QGIS left some NaNs when filling in the Direct Normal Irradiance (dni)
# monthly values in the shapefile (see dni* columns above). Those need to be filled in with real values.
# Without going back to the GIS (which I already ran with the results that are shown above), I am using np.interpolate across
# each column. This is a hack in the sense that the adjacent row do not necessarily mean that the HRUs are adjacent, and it
# is uncertain exactly what is being interpolated, but it is filling in the nan values with real values and allows
# me to move on.
    count = 0
    for ii in range(hrus_dni.shape[0]):
        for jj in range(8,20):
            if np.isnan(hrus_dni.iloc[ii,jj]):
                count += 1
    print(count)
          
    
    for ii in range(hrus_dni.shape[0]):
        for jj in range(8,20):
            hrus_dni.iloc[:,jj] = hrus_dni.iloc[:,jj].interpolate()
    
    print("subsequent coordinates with nan value")
    count = 0
    for ii in range(hrus_dni.shape[0]):
        for jj in range(8,20):
            if np.isnan(hrus_dni.iloc[ii,jj]):
                print(" (", ii,jj, ")", end = '')
                count += 1
                
# Join the HRU slopes from the PRMS parmaeter file (paramdb) to the features from the shapefile.
    hru_slope_df = pd.read_csv(hru_slope_fn)
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
    jan = [0] * 31
    feb = [1] * 29 # assume leap year to get full table
    mar = [2] * 31
    apr = [3] * 30
    may = [4] * 31
    jun = [5] * 30
    jul = [6] * 31
    aug = [7] * 31
    sep = [8] * 30
    octo = [9] * 31
    nov = [10] * 30
    dec = [11] * 31
    month_of_jday = jan + feb + mar + apr + may + jun + jul + aug + sep + octo + nov + dec

# Finally, start the calculation of radadj
    radadj = np.zeros(soltab_df.shape)
    
    min_count = 0
    max_count = 0
    for iday in range(soltab_df.shape[0]):
        imon = month_of_jday[jday]
        for ihru in range(nhru):
            kk = nhru_v11_vals[ihru] - 1
            try:
                radadj[iday,ihru] = solrad_targets_vals_langleys[ihru,imon] * hru_cossl_vals[ihru] / soltab_df_vals[iday,kk]
            except:
                print(iday, ihru, imon, solrad_targets_vals_langleys[ihru,imon], hru_cossl_vals[ihru], soltab_df_vals[iday,kk])
                
            if radadj[iday,ihru] < 0.05:
                radadj[iday,ihru] = 0.05
                min_count += 1
                
            if radadj[iday,ihru] > 0.95:
                radadj[iday,ihru] = 0.05
                max_count += 1
            
#    print ((jday * nhru), min_count, max_count)
    
# End of the calculation of radadj

# hru_id_str = np.empty(nhru, dtype='int')
# for jj in range(nhru):
#     hru_id_str[jj] = jj+1

# # Create the pandas DataFrame. This should be comparable to what has been
# # produced and read in from pptadj_fn
# pptadj_df = pd.DataFrame(foo, columns = hru_id_str)
# pptadj_df["date"] = dates

    
# Here's the cubic spline interpolation between solf and dday
    solf_pts = np.array([.20, .35, .45, .51, .56, .59, .62, .64, .655, .67, \
                         .682, .69, .70, .71, .715, .72, .722, .724, .726, \
                             .728, .73, .734, .738, .742, .746, .75])
    dday_pts = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, \
                         11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, \
                             19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0])
    f2 = interp1d(solf_pts, dday_pts, kind='cubic')
    
#radadj_solf
    radadj_solf = np.zeros(nday * nhru)
    radadj_solf.shape = (nday,nhru)
    
    for iday in range(ndays):
        yr = dates[iday,0]
        mo = dates[iday,1]
        day = dates[iday,2]
        date = datetime(yr,mo,day)
        jday = date.timetuple().tm_yday
        
        for ihru in range(nhru):
#            print(iday,ihru,jday-1)
            try:
                radadj_solf[iday,ihru] = radadj[jday-1,ihru] / pptadj[iday,ihru]
            except Exception as inst:
                print(type(inst))
                print(inst)
                print(iday,ihru,jday-1,radadj.shape, pptadj.shape)
    
# compute the dday values for each time step (day) for each HRU
    dday = np.zeros(ndays*nhru, dtype=np.float64)
    dday.shape = (ndays, nhru)
    
    for iday in range(ndays):
        # yr = dates[iday,0]
        # mo = dates[iday,1]
        # day = dates[iday,2]
        # date = datetime(yr,mo,day)
        # jday = date.dayofyear
        # print(date, jday)
        
        for ihru in range(nhru):
            if radadj_solf[iday, ihru] < 0.2:
                dday[iday,ihru] = 1.0
            elif radadj_solf[iday, ihru] > 0.75:
                dday[iday,ihru] = 26.0
            else:
                dday[iday,ihru] = f2(radadj_solf[iday,ihru])
        
# This computes the dday slope for each day of the record (ie daily values)
# dday_slope(ihru,imon) = (dday(iday,ihru) - dday_intcp(ihru,imon) - 1.0) / tmax_hru(iday,ihru)
# dday = Dday_slope(j, Nowmonth)*Tmax_hru(j) + Dday_intcp(j, Nowmonth) + 1.0
    for iday in range(ndays):
        yr = dates[iday,0]
        mo = dates[iday,1]
        day = dates[iday,2]
        date = datetime(yr,mo,day)
        jday = date.dayofyear
        imon = date.month - 1
        
        for ihru in range(nhru):
            dday_slope_daily[ihru,iday] = (dday[iday,ihru] - dday_intcp[imon,ihru] - 1.0) \
                                        / tmax_hru[iday,ihru]
        
# This computes the monthly values. The parameters are monthly.

    dday_slope = np.zeros(nhru*nmonths, dtype=np.float64)
    dday.shape = (nhru, nmonths)
    count = np.zeros(nmonths, dtype=np.int)
    
    for iday in range(ndays):
        yr = dates[iday,0]
        mo = dates[iday,1]
        day = dates[iday,2]
        date = datetime(yr,mo,day)
        jday = date.dayofyear
        imon = date.month - 1

        for ihru in range(nhru):
            dday_slope[ihru,imon] += dday_slope_daily[ihru,iday]
            count[imon] += 1
            
    for ihru in range(nhru):
        for imon in range(0,nmonths):
            dday_slope[ihru,imon] = dday_slope[ihru,imon] / count[imon]        
    
# save values to file
    f = open(dday_slope_fn, "w")
    f.write("$id,dday_slope\n")

    ii = 1
    for ihru in range(nhru):
        for imon in range(nmonths):
            f.write(str(ii) + "," + str(dday_slope[ihru,imon]) + "\n")
            ii, dday_slope[ihru,imon]
            ii += 1

    f.close()    
    
if __name__ == "__main__":
    
# Lots of stuff needs to be passed in.
# cbh_in_fn='c:/Users/markstro/work1.1/cbh_gridmet/tmin.cbh'
# cbh_out_fn='c:/Users/markstro/work1.1/cbh_gridmet/tmin_out_from_script.cbh'
# nhm_id_fn='c:/Users/markstro/work1.1/cbh_gridmet/nhm_id.txt'
# miss_to_pres_mapping ='c:/Users/markstro/work1.1/cbh_gridmet/miss_to_pres_mapping.csv'
# var_name = "tmin"
# frmt = "%.1f"
    
    # if (len(sys.argv) == 7):     #test whether any arguments have been passed in
    #     cbh_in_fn = sys.argv[1]
    #     cbh_out_fn = sys.argv[2]
    #     nhm_id_fn = sys.argv[3]
    #     miss_to_pres_mapping = sys.argv[4]
    #     var_name = sys.argv[5]
    #     frmt = sys.argv[6]
    # else:
    #     print("No name passed in")
        
        
    main()