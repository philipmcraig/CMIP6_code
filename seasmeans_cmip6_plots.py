# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 13:56:08 2021

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
cmipdir = ncasdir + 'CMIP6/'
indecis = ncasdir + 'INDECIS/'

inst = 'MPI-M'
model = 'MPI-ESM1-2-HR'
mod_lower = model.lower()
mod_lower = mod_lower.replace('-','')
variant = 'r1i1p1f1'
table = ['Amon','Lmon']
VAR = ['tas','pr','mrsos']

fname1 = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/'+VAR[0]+'_'+table[0]+'*')[0]

splits = fname1.split('/')
gridtype = splits[-1].split('_')[-2]

nc1 = xr.open_dataset(fname1)
cmiplat = xr.DataArray(nc1.lat)
cmiplon = xr.DataArray(nc1.lon)
cmiptime = xr.DataArray(nc1.time)
cmiptemp = xr.DataArray(getattr(nc1,VAR[0]))
nc1.close()

fname2 = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/'+VAR[1]+'_'+table[0]+'*')[0]

nc2 = xr.open_dataset(fname2)
cmipprec = xr.DataArray(getattr(nc2,VAR[1]))
nc2.close()

fname3 = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/'+VAR[2]+'_'+table[1]+'*')[0]

nc3 = xr.open_dataset(fname3)
cmipsomo = xr.DataArray(getattr(nc3,VAR[2]))
nc3.close()

year0 = int(str(cmiptime.data[0])[:4])

if year0 == 1950:
    start = 0
elif year0 < 1950:
    start = (1950-year0)*4

cmiptemp_jja = pl.mean(cmiptemp[start+2:-1:4,:,:],axis=0)
cmiptemp_djf = pl.mean(cmiptemp[start:-1:4,:,:],axis=0)

cmipprec_jja = pl.mean(cmipprec[start+2:-1:4,:,:],axis=0)*86400
cmipprec_djf = pl.mean(cmipprec[start:-1:4,:,:],axis=0)*86400

cmipsomo_jja = pl.mean(cmipsomo[start+2:-1:4,:,:],axis=0)/(0.1*1000)
cmipsomo_djf = pl.mean(cmipsomo[start:-1:4,:,:],axis=0)/(0.1*1000)


temp_jja_cyc, c6lon_cyc = util.add_cyclic_point(cmiptemp_jja.data, coord=cmiplon.data)
temp_djf_cyc = util.add_cyclic_point(cmiptemp_djf.data)

prec_jja_cyc = util.add_cyclic_point(cmipprec_jja.data)
prec_djf_cyc = util.add_cyclic_point(cmipprec_djf.data)

somo_jja_cyc = util.add_cyclic_point(cmipsomo_jja.data)
somo_djf_cyc = util.add_cyclic_point(cmipsomo_djf.data)

fig, ax = pl.subplots(3,2,figsize=(8,8))

###############################################################################
proj = ccrs.PlateCarree(central_longitude=0)
ext = [-15,42,35,70]
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        #edgecolor='face',
                                        facecolor='w')
###############################################################################

ax1 = pl.subplot(321,projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax1.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax1.contourf(c6lon_cyc,cmiplat.data,temp_jja_cyc,transform=proj,
                  cmap='YlOrRd',alpha=0.8,levels=pl.linspace(260,310,11),extend='both')
GridLines(ax1,True,False,True,False)
ax1.annotate('JJA '+model+' '+VAR[0],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)
###############################################################################

ax2 = pl.subplot(322,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax2.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax2.contourf(c6lon_cyc,cmiplat.data,temp_djf_cyc,transform=proj,
                  cmap='YlOrRd',alpha=0.8,levels=pl.linspace(260,310,11),extend='both')
GridLines(ax2,True,False,False,False)
ax2.annotate('DJF '+model+' '+VAR[0],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

f = pl.gcf()
colax = f.add_axes([0.91,0.67,0.02,0.27])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('K',fontsize=11)
###############################################################################

ax3 = pl.subplot(323,projection=proj,extent=ext)
ax3.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax3.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax3.contourf(c6lon_cyc,cmiplat.data,prec_jja_cyc,transform=proj,
                  cmap='YlGnBu',alpha=0.8,levels=pl.linspace(0,5,11),extend='max')
GridLines(ax3,False,False,True,False)
ax3.annotate('JJA '+model+' '+VAR[1],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

###############################################################################

ax4 = pl.subplot(324,projection=proj,extent=ext)
ax4.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax4.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax4.contourf(c6lon_cyc,cmiplat.data,prec_djf_cyc,transform=proj,
                  cmap='YlGnBu',alpha=0.8,levels=pl.linspace(0,5,11),extend='max')
GridLines(ax4,False,False,False,False)
ax4.annotate('DJF '+model+' '+VAR[1],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

f = pl.gcf()
colax = f.add_axes([0.91,0.36,0.02,0.27])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('mm day$^{-1}$',fontsize=11,labelpad=-2)
###############################################################################

ax5 = pl.subplot(325,projection=proj,extent=ext)
ax5.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax5.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax5.contourf(c6lon_cyc,cmiplat.data,somo_jja_cyc,transform=proj,
                  cmap='plasma_r',alpha=0.8,levels=pl.linspace(0,0.6,13),extend='max')
GridLines(ax5,False,True,True,False)
ax5.annotate('JJA '+model+' '+VAR[2],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)
###############################################################################

ax6 = pl.subplot(326,projection=proj,extent=ext)
ax6.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ax6.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax6.contourf(c6lon_cyc,cmiplat.data,somo_djf_cyc,transform=proj,
                  cmap='plasma_r',alpha=0.8,levels=pl.linspace(0,0.6,13),extend='max')
GridLines(ax6,False,True,False,False)
ax6.annotate('DJF '+model+' '+VAR[2],(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

f = pl.gcf()
colax = f.add_axes([0.91,0.05,0.02,0.27])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('m$^{3}$ m$^{-3}$',fontsize=11)
###############################################################################

pl.subplots_adjust(top=0.96,bottom=0.03,left=0.06,right=0.90,wspace=0.05,hspace=0.0)

pl.savefig(indecis+'figures/'+model+'_'+VAR[0]+'_'+VAR[1]+'_'+VAR[2]+'_seasmeans.png',
           dpi=400)