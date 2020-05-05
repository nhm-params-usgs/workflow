#!/usr/bin/env python3

# This requires geopandas >= 0.70.0
import geopandas as gpd
import numpy as np
import math
import time

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


def main():
    # ************ Adjust these settings ****************
    workdir = '/home/pnorton/datasets/snodas/unmasked'

    geodatabase = '/home/pnorton/datasets/snodas/unmasked/GIS/GFv1.1_nhrusim_Topology_LATLON.gpkg'
    layer_name = 'GFv1'  # Layer name if geodatabase otherwise None
    shape_key = 'nhru_v11'  # Name of HRU id attribute

    gpkg_filename = '/home/pnorton/datasets/snodas/unmasked/GIS/snodas_grid_NHMpoly_v1.gpkg'
    gpkg_layer_name = 'snodas_grid_NHMpoly_v1'
    gpkg_fld_name = 'grid_horiz'

    weight_file = f'{workdir}/snodas_weights.csv'
    # ************ END of adjustable settings ****************

    print('Starting mem usage:', convert_size(get_mem_usage()))
    print('Reading HRUs')
    gdf = gpd.read_file(geodatabase, layer=layer_name)
    print('  ... done.')
    print('    mem usage:', convert_size(get_mem_usage()))

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

