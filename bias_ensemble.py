# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 15:36:18 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
#import pandas as pd
import glob
from scipy import stats
import cartopy.crs as ccrs
import cartopy.util as util
import cartopy.feature as cfeature
import pcraig_funcs as pc

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
eobsdir = ncasdir + 'EOBS/'

var = 'tas'
var2 = 'tg'
var3 = 't2m'
season = 'JJA'

#inst = ['AS-RCEC','AWI','BCC','CCCma','CMCC','CNRM-CERFACS','CSIRO',
#        'CSIRO-ARCCSS','EC-Earth-Consortium','FIO-QLNM','MIROC','MOHC','MPI-M',
#        'MRI','NCAR','NCC','NOAA-GFDL']
#
#model = ['TaiESM1','AWI-ESM-1-1-LR','BCC-ESM1','CanESM5','CMCC-ESM2',
#         'CNRM-ESM2-1','ACCESS-ESM1-5','ACCESS-CM2','EC-Earth3','FIO-ESM-2-0',
#         'MIROC6','UKESM1-0-LL','MPI-ESM1-2-HR','MRI-ESM2-0','CESM2',
#         'NorESM2-MM','GFDL-ESM4']

allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014/*')

meanfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_ensmean.nc')
lat = xr.DataArray(meanfile.lat)
lon = xr.DataArray(meanfile.lon)
#time = xr.DataArray(meanfile.time)
cmipmean = xr.DataArray(getattr(meanfile,var))
meanfile.close()

eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
eobsdata = xr.DataArray(getattr(eobsfile,var2))
eobsfile.close()

era5file = xr.open_dataset(era5dir+'era5_surf_seasmean_s1.5.nc')
era5data = xr.DataArray(getattr(era5file,var3))
#era5time = xr.DataArray(era5file.time)
era5file.close()

if season == 'JJA':
    cmip_tm = pl.nanmean(cmipmean[1::4,:,:],axis=0)
    eobs_mn = pl.mean(eobsdata[1::4],axis=0) + 273.15
    era5_mn = pl.mean(era5data[2:-14:4],axis=0)
elif season == 'DJF':
    cmip_tm = pl.nanmean(cmipmean[3::4,:,:],axis=0)
    eobs_mn = pl.mean(eobsdata[3::4],axis=0) + 273.15
    era5_mn = pl.mean(era5data[4:-14:4],axis=0)

BIAS_mean = cmip_tm - era5_mn
BIAS_mods = pl.zeros([len(allfiles),cmip_tm.shape[0],cmip_tm.shape[1]])
runtot = pl.zeros_like(BIAS_mean) # running total
CHECK = pl.zeros_like(BIAS_mean); CHECK[:,:] = pl.float32('nan')

for m in range(len(allfiles)):
    filepath = allfiles[m]
    ncfile = xr.open_dataset(filepath)
    cmipdata = xr.DataArray(getattr(ncfile,var))
    ncfile.close()
    
    if season == 'JJA':
        model_tm = pl.nanmean(cmipdata[1::4,:,:],axis=0)
    elif season == 'DJF':
        model_tm = pl.nanmean(cmipdata[3::4,:,:],axis=0)
    
    BIAS_mods[m] = model_tm - era5_mn


for i  in range(BIAS_mean.shape[0]):
    for j in range(BIAS_mean.shape[1]):
        point = BIAS_mods[:,i,j]
        if pl.all(pl.isnan(point)) == True:
            continue
        else:
            pos_count = len(list(filter(lambda x: (x > 0), point)))
            neg_count = len(list(filter(lambda x: (x < 0), point)))
            zer_count = len(list(filter(lambda x: (x == 0), point)))
            
            A = pl.array([pos_count,neg_count,zer_count])
            
            if A.max() > 17.*(2./3.):
                CHECK[i,j] = 1.

CHECK_cyc, lon_cyc = util.add_cyclic_point(CHECK, coord=lon)

proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                                           '50m',edgecolor='grey',
                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')

ax = pl.axes(projection=proj,extent=ext)
ax.coastlines(linewidth=0.5,resolution='50m')
ax.add_feature(ocean_50m,alpha=1,zorder=5)

ax.contourf(lon_cyc,lat,CHECK_cyc,transform=ccrs.PlateCarree(),hatches='...',
            colors='none')