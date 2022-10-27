# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 11:32:42 2020

@author: pmcraig
"""

import pylab as pl
import glob
import xarray as xr
import cartopy.crs as ccrs
import pandas as pd
import cartopy.feature as cfeature
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

teleind = ['nao','ea','scand','eawr']

ncfile1 = xr.open_dataset(ncasdir+'ERA5/era5_surf_mm_1950-1978.nc')
era5lat = xr.DataArray(ncfile1.latitude)
era5lon = xr.DataArray(ncfile1.longitude)
mslp_early = xr.DataArray(ncfile1.msl)
ncfile1.close()

ncfile2 = xr.open_dataset(ncasdir+'ERA5/era5_surf_mm_1979_2017.nc')
mslp_late = xr.DataArray(ncfile2.msl)
ncfile2.close()

#ncfile = xr.open_dataset(era5dir+'era5_surf_mm_1950-2017.nc')
##    if i == 0:
#lat = xr.DataArray(ncfile.latitude)
#lon = xr.DataArray(ncfile.longitude)
#time = xr.DataArray(ncfile.time)
##t2m = xr.DataArray(ncfile.t2m)
#mslp = xr.DataArray(ncfile.msl)
#ncfile.close()

mslp = pl.zeros([mslp_early.shape[0]+mslp_late.shape[0],era5lat.size,
                                                        era5lon.size])

mslp[:mslp_early.shape[0],:,:] = mslp_early.values
mslp[mslp_early.shape[0]:,:,:] = mslp_late.values

mslp = pl.reshape(mslp,newshape=(mslp.shape[0]/12,12,mslp.shape[1],mslp.shape[2]))
mslp_jfm = pl.mean(mslp[:,0:3,:,:],axis=1)
mslp_ao = pl.mean(mslp[:,3:10,:,:],axis=1)
  
tele_jfm = pl.zeros([len(teleind),mslp.shape[0]])
tele_ao = tele_jfm.copy()
for i in range(len(teleind)):
    tele = pd.read_csv(ncasdir+teleind[i]+'_index.tim',header=5,delim_whitespace=True)
    tele = pl.asarray(tele)
    tele = pl.reshape(tele[:,-1],newshape=(tele.shape[0]/12,12))
#    tele_sns = pl.zeros([tele.shape[0],4])
#    tele_sns[:,1] = pl.mean(tele[:,2:5],axis=1) # MAM
#    tele_sns[:,2] = pl.mean(tele[:,5:8],axis=1) # JJA
#    tele_sns[:,3] = pl.mean(tele[:,8:11],axis=1) # SON
#    tele_sns[1:,0] = (tele[:-1,-1]+tele[1:,0]+tele[1:,1])/3 # DJF
#    tele_sns[0,0] = pl.float32('nan')
    tele_jfm[i,:] = pl.mean(tele[:,:3],axis=1)
    tele_ao[i,:] = pl.mean(tele[:,3:11],axis=1)

jfm_std = pl.nanstd(tele_jfm,axis=1)
jfm_mn = pl.nanmean(tele_jfm,axis=1)
ao_std = pl.std(tele_ao,axis=1)
ao_mn = pl.mean(tele_ao,axis=1)

gtr = []; lsr = []
above = pl.zeros([len(teleind),era5lat.size,era5lon.size]); below = above.copy()
pos_anom = above.copy(); neg_anom = below.copy()
for i in range(len(teleind)):
    gtr.append(pl.where(tele_ao[i]>ao_mn[i]+ao_std[i])[0])
    lsr.append(pl.where(tele_ao[i]<ao_mn[i]-ao_std[i])[0])

    above[i] = pl.nanmean(mslp_ao[gtr[i],:,:],axis=0)
    below[i] = pl.nanmean(mslp_ao[lsr[i],:,:],axis=0)
    
    pos_anom[i] = pl.nanmean(mslp_ao[gtr[i],:,:],axis=0) - pl.nanmean(mslp_ao,axis=0)
    neg_anom[i] = pl.nanmean(mslp_ao[lsr[i],:,:],axis=0) - pl.nanmean(mslp_ao,axis=0)

proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                                           '50m',edgecolor='grey',
                                        facecolor='none')
#cmap = 'seismic'
#levs = pl.linspace(-1,1,21)

fig, ax = pl.subplots(4,2,figsize=(8,9.5))

for i in range(8):
    axx = pl.subplot(4,2,i+1,projection=proj,extent=ext)
    axx.coastlines(linewidth=0.5,resolution='50m')
    axx.add_feature(borders_50m,linewidth=0.5,zorder=5)
    #levs = pl.linspace(-20,20,41)#[0,31,59,90,120,151]
    if i in (0,2,4,6):
        #cn = axx.contour(era5lon,era5lat,pos_anom[i/2]/100,transform=ccrs.PlateCarree(),#levels=levs,
        #              colors='k',extend='both',norm=pl.Normalize(-8,8))
        cs = axx.contourf(era5lon,era5lat,pos_anom[i/2]/100,transform=ccrs.PlateCarree(),
                          cmap='PuOr_r',norm=pl.Normalize(-4,4),extend='both',
                            levels=pl.linspace(-4,4,9))
    elif i in (1,3,5,7):
        #cn = axx.contour(era5lon,era5lat,neg_anom[(i-1)/2]/100,transform=ccrs.PlateCarree(),#levels=levs,
        #              colors='k',extend='both',norm=pl.Normalize(-8,8))
        cs = axx.contourf(era5lon,era5lat,neg_anom[(i-1)/2]/100,transform=ccrs.PlateCarree(),
                          cmap='PuOr_r',norm=pl.Normalize(-4,4),extend='both',
                            levels=pl.linspace(-4,4,9))
    #pl.clabel(cn,inline=True,fmt="%.1f",zorder=3,inline_spacing=5,manual=False)
    
    if i == 0:
        GridLines(axx,True,False,True,False)
        axx.annotate('a. $>1$ std GS NAO',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 1:
        GridLines(axx,True,False,False,True)
        axx.annotate('b. $<1$ std GS NAO',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 2:
        GridLines(axx,False,False,True,False)
        axx.annotate('c. $>1$ std GS EA',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 3:
        GridLines(axx,False,False,False,True)
        axx.annotate('d. $<1$ std GS EA',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 4:
        GridLines(axx,False,False,True,False)
        axx.annotate('e. $>1$ std GS SCA',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 5:
        GridLines(axx,False,False,False,True)
        axx.annotate('f. $<1$ std GS SCA',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 6:
        GridLines(axx,False,True,True,False)
        axx.annotate('g. $>1$ std GS EAWR',(-13.8,67.3),bbox={'facecolor':'w'})
    elif i == 7:
        GridLines(axx,False,True,False,True)
        axx.annotate('h. $<1$ std GS EAWR',(-13.8,67.3),bbox={'facecolor':'w'})

pl.subplots_adjust(left=0.005,right=0.995,top=0.97,bottom=0.10,wspace=-0.15,hspace=0.10)

f = pl.gcf()
colax = f.add_axes([0.15,0.05,0.7,0.025])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)#pad=0.05,fraction=0.10,
cb.set_label('hPa',fontsize=12)

#pl.savefig(indecis+'figures/era5_mslp_telecon_anoms_gs.png',dpi=400)