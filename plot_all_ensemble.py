# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 13:29:50 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import pandas as pd
import glob
from scipy import stats
import cartopy.crs as ccrs
import cartopy.util as util
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import pcraig_funcs as pc

def GridLines(ax,top,bottom,left,right):
    """
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
                  zorder=6)
    gl.xlabels_top = top; gl.xlabels_bottom = bottom
    gl.ylabels_left = left; gl.ylabels_right = right
    gl.xlocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20,30,40,50])
    gl.ylocator = mticker.FixedLocator([20,30,40,50,60,70,80])
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'k','size':10}
    gl.ylabel_style = {'color': 'k','size':10}
    
    return None

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
eobsdir = ncasdir + 'EOBS/'

var = 'pr'
var2 = 'rr'
var3 = 'tp'
season = 'JJA'

meanfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_ensmean.nc')
lat = xr.DataArray(meanfile.lat)
lon = xr.DataArray(meanfile.lon)
time = xr.DataArray(meanfile.time)
ensmean = xr.DataArray(getattr(meanfile,var))
meanfile.close()

eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
eobsdata = xr.DataArray(getattr(eobsfile,var2))
eobsfile.close()

era5file = xr.open_dataset(era5dir+'era5_surf_seasmean_s1.5.nc')
era5data = xr.DataArray(getattr(era5file,var3))
era5time = xr.DataArray(era5file.time)
era5file.close()

allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014/*')

modeldata = pl.zeros([len(allfiles),len(time),lat.size,lon.size])
modnames = []

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    cmipdata = xr.DataArray(getattr(ncfile,var))
    modeldata[nci,:,:,:] = cmipdata.data
    ncfile.close
    
    name = allfiles[nci].split('/')[10].split('_')[2]
    modnames.append(name)

if season == 'JJA':
    data_tm = pl.mean(modeldata[:,1::4,:,:],axis=1)*86400
    ensmean = pl.mean(ensmean[1::4,:,:],axis=0)*86400
    eobs_mn = pl.mean(eobsdata[1::4],axis=0)# + 273.15
    era5_mn = pl.mean(era5data[2:-14:4],axis=0)*1000
elif season == 'DJF':
    data_tm = pl.mean(modeldata[:,3::4,:,:],axis=1)*86400
    ensmean = pl.mean(ensmean[3::4,:,:],axis=0)*86400
    eobs_mn = pl.mean(eobsdata[3::4],axis=0)# + 273.15
    era5_mn = pl.mean(era5data[4:-14:4],axis=0)*1000

if var == 'tas':
    meanlevs = pl.linspace(280,300,11); cmap_a = 'RdYlBu_r'
    sprdlevs = pl.linspace(0.5,5,10)
    biaslevs = pl.linspace(-5,5,11); cmap_cd = 'seismic'
    snlevs = pl.linspace(0.6,2,8); cmap_e = 'YlOrRd'
    units = 'K'
elif var == 'pr':
    #data_es_tm = (data_es_tm/data_em_tm)*100 # coefficient of variation
    meanlevs = pl.linspace(0,5,11); cmap_a = 'YlGnBu'
    units = 'mm day$^{-1}$'
    sprdlevs = pl.linspace(0,2,11)#pl.linspace(0.0,2.0,9)
    biaslevs = pl.linspace(-3,3,13); cmap_cd = 'seismic_r'
    snlevs = pl.linspace(-0.5,0.5,11); cmap_e = 'RdBu'  

fig, ax = pl.subplots(4,5,figsize=(17,10))

proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
#borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
#                                           '50m',edgecolor='grey',
#                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')

ensmean_cyc, lon_cyc = util.add_cyclic_point(ensmean.data, coord=lon)
eobs_mn_cyc = util.add_cyclic_point(eobs_mn.data)
era5_mn_cyc = util.add_cyclic_point(era5_mn.data)

for i in range(len(allfiles)):
    
    BIAS = data_tm[i] - era5_mn
    BIAS_cyc = util.add_cyclic_point(BIAS.data)
    #data_tm_cyc = util.add_cyclic_point(data_tm[i])
    
    axx = pl.subplot(4,5,i+1,projection=proj,extent=ext)
    axx.coastlines(linewidth=0.5,resolution='50m')
    axx.add_feature(ocean_50m,alpha=1,zorder=5)
    
    cs = axx.contourf(lon_cyc,lat.data,BIAS_cyc,cmap=cmap_cd,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=biaslevs,
                extend='both')
    
    axx.annotate(modnames[i],(-13.8,68.8),bbox={'facecolor':'w'},zorder=7)
    GridLines(axx,False,False,False,False)

###############################################################################

ax18 = pl.subplot(4,5,18,projection=proj,extent=ext)
ax18.coastlines(linewidth=0.5,resolution='50m')
ax18.add_feature(ocean_50m,alpha=1,zorder=5)

ax18.contourf(lon_cyc,lat.data,ensmean_cyc,cmap=cmap_a,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=meanlevs,
                extend='max')

ax18.annotate('CMIP6 '+season+' '+var+' ensemble mean',(-13.8,68.8),
              bbox={'facecolor':'w'},zorder=7)
GridLines(ax18,False,False,False,False)
###############################################################################

ax19 = pl.subplot(4,5,19,projection=proj,extent=ext)
ax19.coastlines(linewidth=0.5,resolution='50m')
ax19.add_feature(ocean_50m,alpha=1,zorder=5)

ax19.contourf(lon_cyc,lat.data,eobs_mn_cyc,cmap=cmap_a,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=meanlevs,
                extend='max')

ax19.annotate('E-OBS '+season+' '+var2+' mean',(-13.8,68.8),
              bbox={'facecolor':'w'},zorder=7)
GridLines(ax19,False,False,False,False)
###############################################################################

ax20 = pl.subplot(4,5,20,projection=proj,extent=ext)
ax20.coastlines(linewidth=0.5,resolution='50m')
ax20.add_feature(ocean_50m,alpha=1,zorder=5)

ax20.contourf(lon_cyc,lat.data,era5_mn_cyc,cmap=cmap_a,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=meanlevs,
                extend='max')

ax20.annotate('ERA5 '+season+' '+var3+' mean',(-13.8,68.8),
              bbox={'facecolor':'w'},zorder=7)
GridLines(ax20,False,False,False,False)
###############################################################################

f = pl.gcf()
colax = f.add_axes([0.29,0.07,0.42,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=12,labelpad=1)
cb.set_ticks(biaslevs)
cb.set_ticklabels(biaslevs)
cb.ax.tick_params(labelsize=10,pad=2,direction='in',length=0)

pl.subplots_adjust(left=0.01,right=0.99,top=0.99,hspace=0.12,wspace=0.03)

pl.savefig(indecis+'figures/'+var+'_'+season+'_ens_all_panels.png',dpi=375)