# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 15:10:15 2021

@author: pmcraig
"""

import pylab as pl
import glob
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
from scipy import stats
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pcraig_funcs as pc

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'

VAR = 't2m'

nc_reg = xr.open_dataset(era5dir+'era5_surf_seasmean_1950-2017.nc')
era5lat = xr.DataArray(nc_reg.latitude)
era5lon = xr.DataArray(nc_reg.longitude)
era5time_reg = xr.DataArray(nc_reg.time)
era5data_reg = xr.DataArray(getattr(nc_reg,VAR))
nc_reg.close()

nc_lnd = xr.open_dataset(era5dir+'era5_land_seasmean_1981-2021.nc')
era5time_lnd = xr.DataArray(nc_lnd.time)
era5data_lnd = xr.DataArray(getattr(nc_lnd,VAR))
nc_lnd.close()