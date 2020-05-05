#!/usr/bin/env python3

import numpy as np
import pandas as pd
import xarray as xr
import sys
import time

# from helper import np_get_wval
from numpy.ma import masked


def np_get_wval(ndata, wghts, **kwargs):
    nan_default = kwargs.get('nan', np.nan)
    cgrid_ids = wghts['grid_ids'].tolist()
    cwgts = wghts['w'].tolist()
    mdata = np.ma.masked_array(ndata[cgrid_ids], np.isnan(ndata[cgrid_ids]))

    tmp = np.ma.average(mdata, weights=cwgts)
    if tmp is masked:
        return nan_default

    return tmp


def np_unweighted_max(ndata, wghts, **kwargs):
    nan_default = kwargs.get('nan', np.nan)
    cgrid_ids = wghts['grid_ids'].tolist()
    # cwgts = wghts['w'].tolist()
    mdata = np.ma.masked_array(ndata[cgrid_ids], np.isnan(ndata[cgrid_ids]))

    tmp = np.ma.max(mdata)
    if tmp is masked:
        return nan_default

    return tmp


def read_paramdb_file(filename):
    fhdl = open(filename)
    rawdata = fhdl.read().splitlines()
    fhdl.close()
    it = iter(rawdata)

    data = []

    # Skip header
    next(it)

    # Read the parameter values
    for rec in it:
        try:
            data.append(int(rec.split(',')[1]))
        except ValueError:
            data.append(int(float(rec.split(',')[1])))
    return data


def write_parameter(directory, param_name, param_data, nhm_ids):
    with open(f'{directory}/{param_name}.csv', 'w') as ff:
        outstr = f'$id,{param_name}\n'

        for ii, dd in enumerate(param_data.ravel(order='F').tolist()):
            tmp = '{:<20f}'.format(dd).rstrip('0 ')
            if tmp[-1] == '.':
                tmp += '0'
            outstr += '{},{}\n'.format(nhm_ids[ii], tmp)
            # outstr += '{},{}\n'.format(ii+1, tmp)

        ff.write(outstr)


def write_parameter_int(directory, param_name, param_data):
    with open(f'{directory}/{param_name}.csv', 'w') as ff:
        outstr = f'$id,{param_name}\n'

        for ii, dd in enumerate(param_data.ravel(order='F').tolist()):
            outstr += '{},{}\n'.format(ii+1, int(dd))

        ff.write(outstr)


def proc_loop(nhm_ids, hru_ids, src_data, src_var_name, out_size, nan=None):
    dst_data = np.zeros(out_size)

    # 2020-04-13 PAN: order='F' works with vertical test pattern
    #                 order='C' works with horizontal test pattern
    aa = src_data[src_var_name].values.flatten(order="C")

    t1 = time.time()
    for idx, c_nhm_id in enumerate(nhm_ids):
        if c_nhm_id % 1000 == 0:
            sys.stdout.write(f'\rHRU: {c_nhm_id} ({time.time() - t1} seconds)')
            sys.stdout.flush()
            t1 = time.time()

        # print(f'HRU: {c_nhm_id}')
        # get the rows associated with the current nhm_id
        curr_wgts = hru_ids.get_group(c_nhm_id)

        dst_data[idx] = np_get_wval(aa, curr_wgts, nan=nan)
        # dst_data[c_nhm_id-1] = np_get_wval(aa, curr_wgts, nan=nan)

    print()
    return dst_data


def main():
    # ************ Adjust these settings ****************
    # workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/SNODAS/unmasked'
    # paramdb_dir = '/Users/pnorton/tmp/tmp_paramdb'
    #

    # src_file = f'{workdir}/gridtest_v2.nc'
    # src_field = 'LANDMASK'
    # weights_file = f'{workdir}/snodas_weights.csv'

    workdir = './datasets/SNODAS/unmasked'
    paramdb_dir = './tmp/tmp_paramdb'

    # src_file = f'{workdir}/snodas_test_patterns_NHM.nc'
    src_file = f'{workdir}/SWE_median_of_max_yearly_inches_NHM_adj.nc'
    src_field = 'SWE'
    weights_file = f'{workdir}/snodas_weights.csv'

    # nan_value -> value to use when area weighted avg results in a NaN
    nan_value = 1.0
    out_param_name = 'snarea_thresh'
    # ************ END of adjustable settings ****************

    # Read the weights files
    wgt_df = pd.read_csv(weights_file)
    wgt_key = wgt_df.columns[1]
    hru_ids = wgt_df.groupby(wgt_key)
    nhru = len(hru_ids)

    # Get the nhm_id parameter from the parameter database
    # nhm_ids = read_paramdb_file(f'{paramdb_dir}/nhm_id.csv')
    pdb_nhm_ids = read_paramdb_file(f'{paramdb_dir}/nhm_id.csv')
    nhm_ids = sorted(list(set(wgt_df['nhru_v11'].tolist())))

    if len(pdb_nhm_ids) != len(nhm_ids):
        print(f'WARNING: Number of HRUs in {weights_file} does not match number of NHM HRUs')

    try:
        # Read the source data file
        src_data = xr.open_dataset(src_file, mask_and_scale=True)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(f'Processing {out_param_name}')
        # snarea_thresh = proc_loop(nhm_ids, hru_ids, src_data, 'SWE', (nhru), nan=0)
        data = proc_loop(nhm_ids, hru_ids, src_data, src_field, (nhru), nan=1.0)

        # Write the output file
        write_parameter(workdir, out_param_name, data, nhm_ids)
        print('  ...done.')
    except FileNotFoundError:
        print(f'ERROR: File {src_file} does not exist; cannot continue.')
        exit(1)




if __name__ == '__main__':
    main()
