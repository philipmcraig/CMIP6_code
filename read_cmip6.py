# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 15:51:42 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import glob
import cartopy.crs as ccrs
import cartopy.util as util
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import cartopy.util as util
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pcraig_funcs as pc

def GridLines(ax,top,bottom,left,right):
    """
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
                  zorder=3)
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

inst = 'EC-Earth-Consortium'
model = 'EC-Earth3'
mod_lower = model.lower()
mod_lower = mod_lower.replace('-','')
variant = 'r1i1p1f1'
season = 'JJA'

cmip6files = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/*')
#splits = cmip6files[-1].split('/')
#model = splits[9]
#var = splits[-1][:2]

VAR1 = 'tas'
fname = [s for s in cmip6files if VAR1 in s][0]

VAR2 = 'tg'

smth = '1.0'

###############################################################################
ncfile1 = xr.open_dataset(fname)
cmip6lat = xr.DataArray(ncfile1.lat)
cmip6lon = xr.DataArray(ncfile1.lon)
cmip6tas = xr.DataArray(getattr(ncfile1,VAR1))
cmip6time = xr.DataArray(ncfile1.time) # starts with January
ncfile1.close()
print cmip6tas.units

year0 = int(str(cmip6time.data[0])[:4])

#lonswap = pl.zeros_like(cmip6lon)
#lonswap[:cmip6lon.size/2] = cmip6lon[cmip6lon.size/2:]-360
#lonswap[cmip6lon.size/2:] = cmip6lon[:cmip6lon.size/2]
#cmip6lon = lonswap
#
#tas_swap = pl.zeros_like(cmip6tas)
#tas_swap[:,:,:cmip6lon.size/2] = cmip6tas.data[:,:,cmip6lon.size/2:]
#tas_swap[:,:,cmip6lon.size/2:] = cmip6tas.data[:,:,:cmip6lon.size/2]
#cmip6tas = tas_swap

if year0 == 1950 and season == 'JJA':
    start = 2 # slice array from JJA 1950 (index 2 of time dimension)
elif year0 < 1950 and season == 'JJA':
    start = (1950-year0)*4 + 2
elif year0 == 1950 and season == 'DJF':
    start = 4
elif year0 < 1950 and season == 'DJF':
    start = (1950-year0)*4 + 4

tas_jja_mean = pl.mean(cmip6tas[start::4],axis=0)
#tas_mean = util.add_cyclic_point(tas_mean,coord=cmip6lon.data)

###############################################################################
#ncfile2 = xr.open_dataset(indecis+'ncfiles/gtg_season.nc')
#indlon = xr.DataArray(ncfile2.longitude)
#indlat = xr.DataArray(ncfile2.latitude)
#indgtg = xr.DataArray(ncfile2.gtg)
#ncfile2.close()
#
#gtg_jja_mean = pl.mean(indgtg[2:65*4:4],axis=0)

###############################################################################
ncfile3 = xr.open_dataset(eobsdir+VAR2+'_seasmean_remapbil'+smth+'_v23.0e_'+mod_lower+'.nc')
eobsrg_lon = xr.DataArray(ncfile3.lon)
eobsrg_lat = xr.DataArray(ncfile3.lat)
eobsrg_tg = xr.DataArray(getattr(ncfile3,VAR2))
eobsrg_time = xr.DataArray(ncfile3.time) # starts with April
ncfile3.close()
print eobsrg_tg.units

#lonswap = pl.zeros_like(indlon_n96)
#lonswap[:indlon_n96.size/2] = indlon_n96[indlon_n96.size/2:]-360
#lonswap[indlon_n96.size/2:] = indlon_n96[:indlon_n96.size/2]
#indlon_n96 = lonswap
#
#indgtg_n96 = pl.transpose(indgtg_n96.data,axes=(0,2,1))
#tas_swap = pl.zeros_like(indgtg_n96)
#tas_swap[:,:,:indlon_n96.size/2] = indgtg_n96[:,:,indlon_n96.size/2:]
#tas_swap[:,:,indlon_n96.size/2:] = indgtg_n96[:,:,:indlon_n96.size/2]
#indgtg_n96 = tas_swap
if season == 'JJA':
    begin = 1
elif season == 'DJF':
    begin = 3
eobsrg_jja_mean = pl.nanmean(eobsrg_tg[begin::4],axis=0)

###############################################################################
ncfile4 = xr.open_dataset(eobsdir+VAR2+'_seasmean_0.25deg_reg_v23.0e.nc')
eobs_lon = xr.DataArray(ncfile4.longitude)
eobs_lat = xr.DataArray(ncfile4.latitude)
eobs_tg = xr.DataArray(getattr(ncfile4,VAR2)) # starts with April
eobs_time = xr.DataArray(ncfile4.time)
ncfile4.close()
print eobs_tg.units

eobs_jja_mean = pl.nanmean(eobs_tg[begin::4],axis=0)
###############################################################################

if VAR1 == 'tas' and season == 'JJA':
    BIAS = tas_jja_mean - 273.15 - eobsrg_jja_mean
    units = '$^\circ$C'
    levels = pl.linspace(4,32,8)
    blevs = pl.linspace(-5,5,11)
    bcmap = 'seismic'
elif VAR1 == 'tas' and season == 'DJF':
    BIAS = tas_jja_mean - 273.15 - eobsrg_jja_mean
    units = '$^\circ$C'
    levels = pl.linspace(-8,20,8)
    blevs = pl.linspace(-5,5,11)
    bcmap = 'seismic'
elif VAR1 == 'pr':# and season == 'JJA':
    BIAS = tas_jja_mean*86400 - eobsrg_jja_mean
    units = 'mm day$^{-1}$'
    levels = pl.linspace(0,5,11)
    blevs = pl.linspace(-3,3,13)
    bcmap = 'seismic_r'
#elif VAR1 == 'pr' and season == 'DJF':
#    BIAS = tas_jja_mean*86400 - eobsrg_jja_mean
#    units = 'mm day$^{-1}$'
#    levels = pl.linspace(0,5,11)
#    blevs = pl.linspace(-3,3,13)
#    bcmap = 'seismic_r'

tas_mn_cyc, c6lon_cyc = util.add_cyclic_point(tas_jja_mean.data, coord=cmip6lon)
eobsrg_mn_cyc, erg_lon_cyc = util.add_cyclic_point(eobsrg_jja_mean, coord=eobsrg_lon)
BIAS_cyc = util.add_cyclic_point(BIAS.data)

###############################################################################
proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                                           '50m',edgecolor='grey',
                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')
#shpfilename = shpreader.natural_earth(resolution='110m',
#                                      category='cultural',
#                                      name='admin_0_countries')
#reader = shpreader.Reader(shpfilename)
#countries = reader.records()
#
#subunits = shpreader.natural_earth(resolution='10m',category='cultural',name='admin_0_map_subunits')
#reader = shpreader.Reader(subunits)
#things = reader.records()
#X = []
#for stuff in things:
#    #X.append(stuff.attributes['ADMIN'])
#    if stuff.attributes['ADMIN'] == 'Russia' and stuff.attributes['NAME'] == 'Russia':
#        print stuff.geometry.bounds#['NAME']

#X = sorted(X)

#levels = pl.linspace(4,32,8)#pl.linspace(0,5,11)#

fig, ax = pl.subplots(2,2,figsize=(10,7.5))

ax1 = pl.subplot(221,projection=proj,extent=ext) # cmip6 model
ax1.coastlines(linewidth=0.5,resolution='50m',zorder=20)
#ax1.add_feature(ocean_50m)

cs = ax1.contourf(c6lon_cyc,cmip6lat.data,tas_mn_cyc-273.15,cmap='viridis',
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=levels,
                extend='both')
#
#blanks = ['TUR','SYR','TUN','GEO','MAR','DZA','IRQ','CYN']
#for country in countries:
#    if country.attributes['ADM0_A3'] in blanks:
#        ax1.add_geometries(country.geometry, ccrs.PlateCarree(),
#                          facecolor='w',edgecolor='w',#lw=0.3,
#                          label=country.attributes['ADM0_A3'])
#
#for stuff in things:
#    if stuff.attributes['NAME'] == 'Russia' and stuff.attributes['NAME'] == 'Russia':
#        ax1.add_geometries(stuff.geometry,ccrs.PlateCarree(),
#                           facecolor='w',edgecolor='k',lw=0.0)

ax1.annotate(model+' 1950-2014 '+season+' '+VAR1,(-13.8,69.3),bbox={'facecolor':'w'})
GridLines(ax1,True,False,True,False)
cb = pl.colorbar(cs,orientation='horizontal',pad=0.01)
cb.set_label(units,fontsize=11,labelpad=-2)

ax2 = pl.subplot(222,projection=proj,extent=ext) # 0.25 eobs
ax2.coastlines(linewidth=0.5,resolution='50m')
cs = ax2.contourf(eobs_lon,eobs_lat,eobs_jja_mean,cmap='viridis',
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=levels,
                extend='both')
ax2.annotate('EOBS '+season+' 1950-2014 '+VAR2,(-13.8,69.3),bbox={'facecolor':'w'})
GridLines(ax2,True,False,False,True)
cb = pl.colorbar(cs,orientation='horizontal',pad=0.01)
cb.set_label(units,fontsize=11,labelpad=-2)

ax3 = pl.subplot(223,projection=proj,extent=ext) # regridded eobs

pl.subplots_adjust(left=0.05,right=0.95,top=0.97,bottom=0.11,hspace=0.1)
ax3.coastlines(linewidth=0.5,resolution='50m')
cs = ax3.contourf(erg_lon_cyc,eobsrg_lat,eobsrg_mn_cyc,cmap='viridis',
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=levels,
                extend='both')
ax3.annotate('regridded EOBS '+season+' '+VAR2+' \n('+smth+'$\!^\circ\!$ smoothing)',
                             (-13.8,66.8),bbox={'facecolor':'w'})
GridLines(ax3,False,True,True,False)
cb = pl.colorbar(cs,orientation='horizontal',pad=0.08)
cb.set_label(units,fontsize=11,labelpad=-2)

ax4 = pl.subplot(224,projection=proj,extent=ext) # bias
ax4.coastlines(linewidth=0.5,resolution='50m')
cs = ax4.contourf(erg_lon_cyc,eobsrg_lat,BIAS_cyc,cmap=bcmap,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=blevs,
                extend='both')
ax4.add_feature(ocean_50m)
ax4.annotate('model minus observations',(-13.8,69.3),bbox={'facecolor':'w'})
GridLines(ax4,False,True,False,True)
cb = pl.colorbar(cs,orientation='horizontal',pad=0.08)
if VAR1 == 'pr':
    cb.set_ticks(blevs[::2])
elif VAR1 == 'tas':
    cb.set_ticks(blevs)
#cb.set_ticklabels([-3,-2,-1,0,1,2,3])
cb.set_label(units,fontsize=11,labelpad=-0.05)

pl.subplots_adjust(bottom=0.02,hspace=0.05,wspace=0.08)

pl.savefig(indecis+'figures/'+model+'_'+VAR1+'_'+'eobs_smth'+smth+'_bias_'+season+'.png',dpi=375)
#f = pl.gcf()
#colax = f.add_axes([0.15,0.05,0.7,0.025])
#cb = pl.colorbar(cs,orientation='horizontal',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_label('$^\circ$C',fontsize=12,labelpad=-0.05)