# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:38:24 2020

@author: markstro
"""

import numpy as np
import sys

#markstro_param_path = "c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/"
#tmin_path = 'c:/Users/markstro/work1.1/cbh_gridmet/tmin.cbh'

# Read the cbh file
def read_cbh(cbh_fn):
    foo = np.loadtxt(cbh_fn, skiprows=3,)

    dates = foo[:,0:6].astype('int')
    vals = foo[:,6:foo.shape[1]]
    
    return dates, vals


def main(markstro_param_path, tmin_path):
    dates, tmin = read_cbh(tmin_path)
    tmin_f = tmin * 1.8 +32
    
    nhru = tmin_f.shape[1]
 
# DO j = 1, Active_hrus
#          elh = (597.3-(0.5653*Tavgc(i)))*2.54
#          Potet(i) = Jh_coef(i, Nowmonth)*(Tavgf(i)-Jh_coef_hru(i))*Swrad(i)/elh
#          IF ( Potet(i)<0.0) Potet(i) = 0.0

# Based on what is above, the parameter jh_coef_hru will be set the the lowest tmin value for each HRU for the
# POR of the weather data.

    jh_coef_hru = np.zeros(nhru)
    
    for ihru in range(nhru):
        jh_coef_hru[ihru] = min(tmin_f[:,ihru])
        jh_coef_hru[ihru] = jh_coef_hru[ihru] - 1.0
    
    f = open(markstro_param_path + 'jh_coef_hru.csv', "w")
    f.write("$id,jh_coef_hru\n")
    
    for jj in range(nhru):
        f.write(str(jj+1) + "," + str(jh_coef_hru[jj]) + '\n')
        
    f.close()


if __name__ == "__main__":
    
# markstro_param_path = "c:/Users/markstro/work1.1/paramdb_v1.1/paramdb_markstro/"
# tmin_path = 'c:/Users/markstro/work1.1/cbh_gridmet/tmin.cbh'
    
    if (len(sys.argv) == 3):     #test whether any arguments have been passed in
        markstro_param_path = sys.argv[1]
        tmin_path = sys.argv[2]
    else:
        print("No name passed in")
        
        
    main(markstro_param_path, tmin_path)
