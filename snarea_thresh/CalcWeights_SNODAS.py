#!/usr/bin/env python3

# This requires geopandas >= 0.70.0
import geopandas as gpd
import pandas as pd
import numpy as np
import math
import time

from shapely.geometry import Polygon
import resource
import sys


def get_mem_usage():
    mem_size = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    if sys.platform == 'darwin':
        # mac OS returns ru_maxrss in bytes (Linux returns it in kilobytes)
        mem_size *= 1024
    return mem_size


def convert_size(size_kbytes):
    # from https://python-forum.io/Thread-Convert-file-sizes-will-this-produce-accurate-results
    if size_kbytes == 0:
        return "0 KB"

    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_kbytes, 1024)))
    power = math.pow(1024, i)
    size = round(size_kbytes / power, 2)
    return f'{size} {size_name[i-1]}'


def create_polygons(lats, lons, grid_spacing, data):
    t1 = time.time()
    lon, lat = np.meshgrid(lons, lats)
    print(f'meshgrid: {time.time() - t1}')

    # For weight generation the datahandle variable is not used for anything
    # but creating the geodataframe of the source netcdf file
    t1 = time.time()
    df = pd.DataFrame({'grid': data})
    print(f'DataFrame: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    # res is half of the 'resolution' (e.g. gridspacing in degrees)
    res = grid_spacing / 2.0
    poly = []
    index = []
    count = 0

    # Create polygon features of the grid-cell bounding boxes
    t1 = time.time()
    for i in range(np.shape(lat)[0]):
        for j in range(np.shape(lon)[1]):
            lat_point_list = [lat[i, j] - res, lat[i, j] + res, lat[i, j] + res, lat[i, j] - res]
            lon_point_list = [lon[i, j] + res, lon[i, j] + res, lon[i, j] - res, lon[i, j] - res]
            poly.append(Polygon(zip(lon_point_list, lat_point_list)))
            index.append(count)
            count += 1
    print(f'poly: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    t1 = time.time()
    ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly, crs='epsg:4326')
    print(f'GeoDataFrame: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))
    # print(ncfcells.memory_usage(deep=True))
    print('..done.')

    return ncfcells


def create_polygons_v2(lats, lons, grid_spacing, data):

    def get_poly(row, lats, lons):
        ii = row.name * 4
        return Polygon(zip(lons[ii:ii+4], lats[ii:ii+4]))

    t1 = time.time()
    lon, lat = np.meshgrid(lons, lats)
    print(f'meshgrid: {time.time() - t1}')

    # res is half of the 'resolution' (e.g. gridspacing in degrees)
    res = grid_spacing / 2.0

    # For weight generation the datahandle variable is only used for creating
    # the geodataframe from the source netcdf file.
    t1 = time.time()
    df = pd.DataFrame({'grid': data})
    print(f'DataFrame: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    # Create polygon features of the grid-cell bounding boxes
    t1 = time.time()
    lats1 = lat[:, :] - res
    lats2 = lat[:, :] + res

    yy = np.dstack((lats1, lats2, lats2, lats1)).flatten()

    lons1 = lon[:, :] + res
    lons2 = lon[:, :] - res

    xx = np.dstack((lons1, lons1, lons2, lons2)).flatten()
    print(f'numpy: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    t1 = time.time()
    df['geometry'] = df.apply(get_poly, lats=yy, lons=xx, axis=1)
    print(f'df_geometry: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    t1 = time.time()
    ncfcells = gpd.GeoDataFrame(df, crs='epsg:4326')
    print(f'GeoDataFrame: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))
    # print(ncfcells.memory_usage(deep=True))
    print('..done.')

    return ncfcells


def main():
    # ************ Adjust these settings ****************
    workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked/test_on_subset'

    geodatabase = f'{workdir}/nhruv11_sim30_WGS84_subset1.gpkg'
    layer_name = 'nhruv11_sim30_WGS84_subset1'
    shape_key = 'nhru_v11'  # Name of HRU id attribute

    gpkg_filename = f'{workdir}/snodas_grid_subset1.gpkg'
    gpkg_layer_name = 'snodas_grid_subset1'
    gpkg_fld_name = 'grid_horiz'

    weight_file = f'{workdir}/snodas_weights.csv'
    # ************ END of adjustable settings ****************

    print('Starting mem usage:', convert_size(get_mem_usage()))
    print('Reading HRUs')
    gdf = gpd.read_file(geodatabase, layer=layer_name)
    print('  ... done.')
    print('    mem usage:', convert_size(get_mem_usage()))

    # The geospatial fabric v1.1 is in Albers, we need to re-project it to WGS84.
    # print(f'Re-projecting {layer_name} to WGS84 (lat/lon)')
    # gdf = gdf.to_crs(epsg=4326)
    # print('  ...done.')
    # print('    mem usage:', convert_size(get_mem_usage()))

    t1 = time.time()
    print('Reading source grid')
    ncfcells = gpd.read_file(gpkg_filename, layer=gpkg_layer_name)
    print(f'  ncfcells: {time.time() - t1}')
    print('     mem usage:', convert_size(get_mem_usage()))
    print('...done.')

    t1 = time.time()
    spatial_index = ncfcells.sindex
    print(f'ncfcells.sindex: {time.time() - t1}')
    print('    mem usage:', convert_size(get_mem_usage()))

    print('='*40)
    outer_t1 = time.time()
    t1 = time.time()

    sys.stdout.write('\nIndex:       ')
    sys.stdout.flush()

    with open(weight_file, 'w') as f:
        # gdf HRU shapefile
        # ncfcells is source netcdf grid geometry

        f.write(f'grid_ids,{shape_key},w\n')

        for index, row in gdf.iterrows():
            if index % 500 == 0:
                sys.stdout.write(f'\nIndex: {index}     ({time.time() - t1} seconds)')
                sys.stdout.flush()
                t1 = time.time()

            count = 0

            # Restrict matches first to the HRU geometry bounding box
            possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))

            if len(possible_matches_index) != 0:
                possible_matches = ncfcells.iloc[possible_matches_index]
                precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]

                if len(precise_matches) != 0:
                    res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')
                    hru_area = gdf.loc[[index], 'geometry'].area

                    for nindex, irow in res_intersection.iterrows():
                        tmpfloat = np.float(res_intersection.area.iloc[nindex] / hru_area)
                        srcidx = np.int(precise_matches.iloc[count][gpkg_fld_name])

                        f.write(f'{srcidx},{np.int(irow[shape_key])},{tmpfloat}\n')
                        count += 1

    print(f'outer_time: {time.time() - outer_t1}')

    print('Ending mem usage:', convert_size(get_mem_usage()))
    

if __name__ == '__main__':
    main()

