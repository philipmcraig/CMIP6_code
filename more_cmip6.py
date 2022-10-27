# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:03:06 2021

@author: pmcraig
"""
import pylab as pl
import glob
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
from scipy import stats
import cartopy.feature as cfeature
import cartopy.util as util
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pcraig_funcs as pc

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
era5dir = ncasdir + 'ERA5/'

#inst = 'CMCC'
#model = 'CMCC-ESM2'
#mod_lower = model.lower()
#mod_lower = mod_lower.replace('-','')
#variant = 'r1i1p1f1'
#table = 'Amon'
VAR = 'psl'
VAR2 = 'msl'

ensfile = xr.open_dataset(cmipdir+'standard_grid/'+VAR+'_cmip6_ensmean.nc')
cmiplat = xr.DataArray(ensfile.lat)
cmiplon = xr.DataArray(ensfile.lon)
cmiptime = xr.DataArray(ensfile.time)
cmipdata = xr.DataArray(getattr(ensfile,VAR))
ensfile.close()

erafile = xr.open_dataset(era5dir+'/era5_msl_seasmean_s1.5_full.nc')
era5data = xr.DataArray(getattr(erafile,VAR2))
era5time = xr.DataArray(erafile.time)
erafile.close()

era5_jja_mn = pl.nanmean(era5data[2::4],axis=0)/100
era5_djf_mn = pl.nanmean(era5data[4::4],axis=0)/100

allfiles = glob.glob(cmipdir+'standard_grid/'+VAR+'_1950-2014_new/*')

modeldata = pl.zeros([len(allfiles),len(cmiptime),cmiplat.size,cmiplon.size])
modnames = []

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,VAR))
    ncfile.close()
    
    split1 = allfiles[nci].split('/')
    split2 = split1[-1].split('_')[2]
    modnames.append(split2)
    
modnames[1] = modnames[1][:-2]
modnames[2] = modnames[2][:7]
modnames[6] = modnames[6][:-2]
modnames[9] = modnames[9][:-4]
modnames[12] = modnames[12][:-5]
modnames[13] = modnames[13][:-2]
modnames[14] = modnames[14][:-3]
modnames[16] = modnames[16][:6]

mods_jja = pl.nanmean(modeldata[:,1::4],axis=0)/100
mods_djf = pl.nanmean(modeldata[:,3::4],axis=0)/100


data_jja_mn = pl.nanmean(cmipdata[1:-1:4,:,:],axis=0)/100
data_djf_mn = pl.nanmean(cmipdata[3:-1:4,:,:],axis=0)/100

data_jja_cyc, lon_cyc = util.add_cyclic_point(data_jja_mn, coord=cmiplon.data)
data_djf_cyc = util.add_cyclic_point(data_djf_mn)

mods_jja_cyc = util.add_cyclic_point(mods_jja,axis=2)
mods_djf_cyc = util.add_cyclic_point(mods_djf,axis=2)

era5_jja_cyc = util.add_cyclic_point(era5_jja_mn)
era5_djf_cyc = util.add_cyclic_point(era5_djf_mn)

#fname1 = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/'+VAR+'_'+table+'*')[0]
#
#splits = fname1.split('/')
#gridtype = splits[-1].split('_')[-2]
#
#nc1 = xr.open_dataset(fname1)
#cmiplat = xr.DataArray(nc1.lat)
#cmiplon = xr.DataArray(nc1.lon)
#cmiptime = xr.DataArray(nc1.time)
#cmippres = xr.DataArray(getattr(nc1,VAR))
#nc1.close()

proj = ccrs.PlateCarree()
ext = [-65,45,30,75]
ext2 = [-15,42,35,70]
borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                                           '50m',edgecolor='grey',
                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')

#ax = pl.axes(projection=proj,extent=ext)
#ax.coastlines(linewidth=0.5,resolution='50m')
#
#cs = ax.contourf(lon_cyc,cmiplat.data,data_djf_cyc)
#pl.colorbar(cs)

fig, ax  = pl.subplots(6,3,figsize=(9,9.5))

ax1 = pl.subplot(6,3,18,projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m')
cs = ax1.contourf(lon_cyc,cmiplat.data,data_djf_cyc,norm=pl.Normalize(988,1028),
            extend='both',levels=pl.linspace(988,1028,11))

ax1.annotate('DJF CMIP6 mean',(-62.0,74.5),
             bbox={'facecolor':'w'},zorder=7,size=8)

f = pl.gcf()
colax = f.add_axes([0.68,0.04,0.30,0.02])
ca = pl.colorbar(cs,orientation='horizontal',cax=colax)
ca.set_label('hPa')
ca.set_ticks(pl.linspace(988,1028,6))
ca.ax.tick_params(direction='in')

for i in range(0,17):
    axx = pl.subplot(6,3,i+1,projection=proj,extent=ext)
    axx.coastlines(linewidth=0.5,resolution='50m')
    cf = axx.contourf(lon_cyc,cmiplat.data,mods_djf_cyc[i]-data_djf_cyc,
                 cmap='PRGn_r',norm=pl.Normalize(-3,3),extend='both',
                    levels=pl.linspace(-3,3,13))
    
    axx.annotate(modnames[i],(-62.0,74.5),
             bbox={'facecolor':'w'},zorder=7,size=8)

f = pl.gcf()
colax = f.add_axes([0.03,0.04,0.6,0.02])
cb = pl.colorbar(cf,orientation='horizontal',cax=colax)
cb.set_label('hPa')
cb.set_ticks(pl.linspace(-3,3,7))
#cb.set_ticklabels(pl.linspace(-3,3,13).astype('int'))
cb.ax.tick_params(direction='in')

pl.tight_layout()
pl.subplots_adjust(bottom=0.07,hspace=0.06)

#pl.savefig(indecis+'figures/cmip6_djf_'+VAR+'_mod_anoms_ensmean.png')