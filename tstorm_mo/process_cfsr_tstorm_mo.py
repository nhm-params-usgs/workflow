#!/usr/bin/env python3

import numpy as np
import pandas as pd
import xarray as xr
import sys

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


def write_parameter(directory, param_name, param_data):
    with open(f'{directory}/{param_name}.csv', 'w') as ff:
        outstr = f'$id,{param_name}\n'

        for ii, dd in enumerate(param_data.ravel(order='F').tolist()):
            tmp = '{:<20f}'.format(dd).rstrip('0 ')
            if tmp[-1] == '.':
                tmp += '0'
            outstr += '{},{}\n'.format(ii+1, tmp)

        ff.write(outstr)


def write_parameter_int(directory, param_name, param_data):
    with open(f'{directory}/{param_name}.csv', 'w') as ff:
        outstr = f'$id,{param_name}\n'

        for ii, dd in enumerate(param_data.ravel(order='F').tolist()):
            outstr += '{},{}\n'.format(ii+1, int(dd))

        ff.write(outstr)


def proc_loop_bymonth(nhm_ids, hru_ids, src_data, src_var_name, out_size, nan=None):
    dst_data = np.zeros(out_size)

    for month in range(12):
        sys.stdout.write(f'\rMonth {month+1}                 ')
        sys.stdout.flush()
        # print(f'Month: {month+1}')
        aa = src_data[src_var_name][month].values.flatten(order="K")

        for c_nhm_id in nhm_ids:
            # get the rows associated with the current nhm_id
            curr_wgts = hru_ids.get_group(c_nhm_id)

            dst_data[c_nhm_id-1, month] = np_unweighted_max(aa, curr_wgts, nan=nan)

            if c_nhm_id % 10000 == 0:
                sys.stdout.write(f'\rMonth {month+1}; HRU: {c_nhm_id}')
                sys.stdout.flush()
                # print(c_nhm_id)
    print()
    return dst_data


def main():
    # ************ Adjust these settings ****************
    workdir = '/Users/pnorton/Projects/National_Hydrology_Model/datasets/tstorm_mo'
    paramdb_dir = '/Users/pnorton/tmp/tmp_paramdb'

    src_file = f'{workdir}/out1.nc'
    weights_file = f'{workdir}/cfsr_weights.csv'
    # ************ END of adjustable settings ****************

    # Get the nhm_id parameter from the parameter database
    nhm_ids = read_paramdb_file(f'{paramdb_dir}/nhm_id.csv')

    # Read the weights files
    wgt_df = pd.read_csv(weights_file)
    wgt_key = wgt_df.columns[1]
    hru_ids = wgt_df.groupby(wgt_key)
    nhru = len(hru_ids)

    # Read the source data file
    src_data = xr.open_dataset(src_file, mask_and_scale=True)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('Processing tstorm_mo')
    tstorm_mo = proc_loop_bymonth(nhm_ids, hru_ids, src_data, 'TSTORM_MO', (nhru, 12), nan=0)

    # Write the output file
    write_parameter_int(workdir, 'tstorm_mo', tstorm_mo)
    print('  ...done.')


if __name__ == '__main__':
    main()
