#!/usr/bin/env python3

# Create a new netcdf file, adding two grid arrays numbered 1..n in 1) C-order,
# 2) Fortran-order.

import numpy as np
import xarray as xr

workdir = './datasets/SNODAS/unmasked'
src_netcdf = f'{workdir}/snodas_landmask_NHM.nc'
src_mask_var = 'LANDMASK'

dst_netcdf = f'{workdir}/snodas_masks_NHM.nc'

ds2 = xr.open_dataset(src_netcdf, mask_and_scale=False)

# Create a grid array the same size and shape as the source netcdf variable.
# grid array is just a monotonic series
grid_v = np.arange(1, ds2[src_mask_var].values.size+1).reshape((ds2[src_mask_var].values.shape), order='F')
grid_h = np.arange(1, ds2[src_mask_var].values.size+1).reshape((ds2[src_mask_var].values.shape), order='C')

gr_vert = xr.DataArray(grid_v, coords={'lat': ds2['lat'].values, 'lon': ds2['lon'].values},
                       dims=['lat', 'lon'])

gr_horiz = xr.DataArray(grid_h, coords={'lat': ds2['lat'].values, 'lon': ds2['lon'].values},
                        dims=['lat', 'lon'])

ds3 = xr.merge([ds2, xr.Dataset({'grid_vert': gr_vert})])
ds3 = xr.merge([ds3, xr.Dataset({'grid_horiz': gr_horiz})])

ds3.to_netcdf(dst_netcdf, encoding={'lat': {'_FillValue': None}, 'lon': {'_FillValue': None}, })
