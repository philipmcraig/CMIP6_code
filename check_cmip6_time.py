# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:08:22 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import glob
#import cartopy.crs as ccrs
#import cartopy.util as util
#import cartopy.io.shapereader as shpreader
#import cartopy.feature as cfeature
#import cartopy.util as util
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pcraig_funcs as pc

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
cmipdir = ncasdir + 'CMIP6/'

nc1 = xr.open_dataset(ncasdir+'CMIP6/standard_grid/tas/' + \
                        'tas_Amon_CESM2_historical_r1i1p1f1_s1.5_seasmean.nc')
time1 = xr.DataArray(nc1.time)
nc1.close()

nc2 = xr.open_dataset(ncasdir+'CMIP6/standard_grid/tas_1950-2014/' + \
                        'tas_Amon_ACCESS-ESM1-5_historical_r1i1p1f1_s1.5_seasmean.nc')
time2 = xr.DataArray(nc2.time)
nc2.close()

print time1.data[0]
print time2.data[0]