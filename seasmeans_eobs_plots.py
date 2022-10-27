# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 11:25:23 2021

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

def GridLines(ax,top,bottom,left,right):
    """
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
                  zorder=15)
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
eobsdir = ncasdir + 'EOBS/'
indecis = ncasdir + 'INDECIS/'

VAR = ['tg','rr']

nc1 = xr.open_dataset(eobsdir+VAR[0]+'_seasmean_0.25deg_reg_v23.0e.nc')
eobslon = xr.DataArray(nc1.longitude)
eobslat = xr.DataArray(nc1.latitude)
eobstime = xr.DataArray(nc1.time)
eobs_temp = xr.DataArray(getattr(nc1,VAR[0]))
nc1.close()

nc2 = xr.open_dataset(eobsdir+VAR[1]+'_seasmean_0.25deg_reg_v23.0e.nc')
eobs_prec = xr.DataArray(getattr(nc2,VAR[1]))
nc2.close()

eobstemp_jja = pl.nanmean(eobs_temp[1::4,:,:],axis=0) + 273.15
eobstemp_djf = pl.nanmean(eobs_temp[3::4,:,:],axis=0) + 273.15

eobsprec_jja = pl.nanmean(eobs_prec[1::4,:,:],axis=0)
eobsprec_djf = pl.nanmean(eobs_prec[3::4,:,:],axis=0)


fig, ax = pl.subplots(2,2,figsize=(8,5))

###############################################################################
proj = ccrs.PlateCarree(central_longitude=0)
ext = [-15,42,35,70]
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        #edgecolor='face',
                                        facecolor='w')
###############################################################################

ax1 = pl.subplot(221,projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax1.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax1.contourf(eobslon.data,eobslat.data,eobstemp_jja,transform=proj,
                  cmap='YlOrRd',alpha=0.8,levels=pl.linspace(260,310,11),extend='both')
GridLines(ax1,True,False,True,False)
ax1.annotate('JJA E-OBS '+VAR[0],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)
###############################################################################

ax2 = pl.subplot(222,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax2.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax2.contourf(eobslon.data,eobslat.data,eobstemp_djf,transform=proj,
                  cmap='YlOrRd',alpha=0.8,levels=pl.linspace(260,310,11),extend='both')
GridLines(ax2,True,False,False,False)
ax2.annotate('DJF E-OBS '+VAR[0],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

f = pl.gcf()
colax = f.add_axes([0.91,0.51,0.02,0.42])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('K',fontsize=11)
###############################################################################

ax3 = pl.subplot(223,projection=proj,extent=ext)
ax3.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax3.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax3.contourf(eobslon.data,eobslat.data,eobsprec_jja,transform=proj,
                  cmap='YlGnBu',alpha=0.8,levels=pl.linspace(0,5,11),extend='max')
GridLines(ax3,False,True,True,False)
ax3.annotate('JJA E-OBS '+VAR[1],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

###############################################################################

ax4 = pl.subplot(224,projection=proj,extent=ext)
ax4.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax4.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax4.contourf(eobslon.data,eobslat.data,eobsprec_djf,transform=proj,
                  cmap='YlGnBu',alpha=0.8,levels=pl.linspace(0,5,11),extend='max')
GridLines(ax4,False,True,False,False)
ax4.annotate('DJF E-OBS '+VAR[1],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

f = pl.gcf()
colax = f.add_axes([0.91,0.05,0.02,0.42])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('mm day$^{-1}$',fontsize=11,labelpad=6)
###############################################################################

pl.subplots_adjust(top=0.96,bottom=0.03,left=0.06,right=0.90,wspace=0.05,hspace=0.0)

#pl.savefig(indecis+'figures/eobs_'+VAR[0]+'_'+VAR[1]+'_seasmeans.png',
 #          dpi=400)