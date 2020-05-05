#!/usr/bin/env python3

# This requires geopandas >= 0.70.0
import geopandas as gpd
import pandas as pd
import numpy as np

import xarray as xr
import csv
from shapely.geometry import Polygon
import multiprocessing as mp


def process_weights(shp_chunk, shape_key, src_netcdf):
    results = []

    # Have to set spatial_index here instead of passing it in as an argument,
    # otherwise the intersections don't work (all return empty sets)
    # see: https://github.com/geopandas/geopandas/issues/1259
    spatial_index = src_netcdf.sindex

    for index, row in shp_chunk.iterrows():
        count = 0

        possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))

        if len(possible_matches_index) != 0:
            possible_matches = src_netcdf.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

            if len(precise_matches) != 0:
                res_intersection = gpd.overlay(shp_chunk.loc[[index]], precise_matches, how='intersection')

                for nindex, irow in res_intersection.iterrows():
                    tmpfloat = np.float(res_intersection.area.iloc[nindex] / shp_chunk.loc[[index], 'geometry'].area)
                    results.append([np.int(precise_matches.index[count]), np.int(irow[shape_key]), tmpfloat])
                    count += 1
    return results


def mp_process_weights(cell_geom, src_polys, shape_key):
    cpus = 6

    isect_dist = []

    shp_chunk = np.array_split(src_polys, cpus)

    pool = mp.Pool(processes=cpus)

    chunk_processes = [pool.apply_async(process_weights, args=(chunk, shape_key, cell_geom)) for chunk in shp_chunk]

    intersect_results = [chunk.get() for chunk in chunk_processes]
    isect_dist.append(intersect_results)

    return isect_dist


def main():
    # ************ Adjust these settings ****************
    workdir = './datasets/tstorm_mo'
    src_netcdf = f'{workdir}/out1.nc'
    src_gridspacing = 0.5  # Gridspacing of src_netcdf in degrees
    geodatabase = './Trans-boundary_HRUs/GIS/GFv1.1_v2e.gdb'
    layer_name = 'nhru_v11'  # Layer name if geodatabase otherwise None
    shape_key = 'nhru_v11'  # Name of HRU id attribute
    weight_file = f'{workdir}/cfsr_weights.csv'
    # ************ END of adjustable settings ****************

    # Open the source netcdf file
    ds = xr.open_dataset(src_netcdf, mask_and_scale=False)

    # Change the variable names to match the src_netcdf file
    lathandle = ds['lat']
    lonhandle = ds['lon']
    datahandle = ds['TSTORM_MO']

    print('Reading HRUs')
    gdf = gpd.read_file(geodatabase, layer=layer_name)
    print('  ... done.')

    # The geospatial fabric v1.1 is in Albers, we need to re-project it to WGS84.
    print(f'Re-projecting {layer_name} to WGS84 (lat/lon)')
    gdf = gdf.to_crs(epsg=4326)
    print('  ...done.')

    # Print some information on the data
    print('\n Data attributes, sizes, and coords \n')
    print('\n Data sizes are: \n', datahandle.sizes)
    print('\n Data coords are: \n', datahandle.coords)

    # Plot a single timestep of the data
    datahandle.isel(time=0).plot()

    lon, lat = np.meshgrid(lonhandle, lathandle)

    # For weight generation the datahandle variable is not used for anything
    # but creating the geodataframe of the source netcdf file
    df = pd.DataFrame({'grid': datahandle.isel(time=0).values.flatten()})

    # res is half of the 'resolution' (e.g. gridspacing in degrees)
    res = src_gridspacing / 2.0
    poly = []
    index = []
    count = 0

    # Create polygon features of the grid-cell bounding boxes
    for i in range(np.shape(lat)[0]):
        for j in range(np.shape(lon)[1]):
            lat_point_list = [lat[i, j] - res, lat[i, j] + res, lat[i, j] + res, lat[i, j] - res]
            lon_point_list = [lon[i, j] + res, lon[i, j] + res, lon[i, j] - res, lon[i, j] - res]
            poly.append(Polygon(zip(lon_point_list, lat_point_list)))
            index.append(count)
            count += 1

    ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly, crs='epsg:4326')

    # Use multiple cores to process the weights
    final_res = mp_process_weights(ncfcells, gdf, shape_key)

    # Write the results to a weights file
    with open(weight_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['grid_ids', shape_key, 'w'])

        for pp in final_res:
            for xx in pp:
                for dd in xx:
                    writer.writerow([dd[0], dd[1], dd[2]])


if __name__ == '__main__':
    main()
