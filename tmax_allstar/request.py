#!/usr/bin/env python3

import cdsapi

# To suppress repeated InsecureRequestWarning, export the following before running this script.
# export PYTHONWARNINGS=once

# req_year = 2012
req_years = [2013, 2014, 2015, 2016, 2017, 2018, 2019]

c = cdsapi.Client()

for req_year in req_years:
    print(f'Working on year: {req_year}')
    c.retrieve('reanalysis-era5-single-levels',
               {'product_type': 'reanalysis',
                'variable': ['2m_temperature', 'precipitation_type'],
                'year': req_year,
                'month': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'],
                'day': ['01', '02', '03', '04', '05', '06', '07', '08',
                        '09', '10', '11', '12', '13', '14', '15', '16',
                        '17', '18', '19', '20', '21', '22', '23', '24',
                        '25', '26', '27', '28', '29', '30', '31'],
                'area': [54, -135, 22, -64],  # N, W, S, E. Default: global
                'time': ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                         '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                         '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                         '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
                'format': 'netcdf'},
               f'download_{req_year}.nc')
