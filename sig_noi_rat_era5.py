# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:03:12 2021

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

def CombineEarlyLate(early,late,era5lat,era5lon):
    """
    """
    varfull = pl.zeros([early.shape[0]+late.shape[0],era5lat.size,era5lon.size])
    
    varfull[:early.shape[0],:,:] = early
    varfull[early.shape[0]:,:,:] = late
    
    varfull = pl.reshape(varfull,newshape=(varfull.shape[0]/12,12,
                                           varfull.shape[1],varfull.shape[2]))

    return varfull
    
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

def S_N_ratio(data_jja,data_djf,lat,lon,gmst_smth_all):
    """
    """
    data_all = [data_jja,data_djf[1:]]
    # apply the linear model
    # L(t) = alpha*G(t) + beta
    L = [pl.zeros_like(data_jja[:,:,:]),pl.zeros_like(data_djf[1:,:,:])]
    # residuals array
    R = [pl.zeros_like(data_jja[:,:,:]),pl.zeros_like(data_djf[1:,:,:])]
    
    for season in range(2):
        for i in range(lat.size):# loop over lon:
            for j in range(lon.size):# loop over lat:
                # linear regression at each grid point
                GMST = gmst_smth_all[season]
                data_season = data_all[season]
                out = stats.linregress(GMST,data_season[:,i,j])
                slope = out[0]
                intercept = out[1]
                # GMST regressed onto the variable
                
                # calculate L(t) at each grid point
                # alpha & beta are the regresion coefficients
                # G(t) is GMST
                L[season][:,i,j] = slope*GMST + intercept
                
                # calculate the residuals at each grid point
                # variable minus L(t)
                R[season][:,i,j] = data_season[:,i,j] - L[season][:,i,j]
    
    
    # calculate the signal
    # S = L[timeN] - L[time0]
    signal_jja = L[0][-1] - L[0][0]
    signal_djf = L[1][-1] - L[1][0]
    
    # calculate the noise
    # N = std(residual)
    noise_jja = pl.std(R[0],axis=0)
    noise_djf = pl.std(R[1],axis=0)
    
    # signal-to-noise ratio
    # S/N
    sig_noi_ratio_jja = signal_jja/noise_jja
    sig_noi_ratio_djf = signal_djf/noise_djf
    
    #SN_ratios = pl.array([sig_noi_ratio_jja,sig_noi_ratio_djf])
    
    return sig_noi_ratio_jja, sig_noi_ratio_djf

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'

ctgry = 'land'
VAR = 'swvl1'

ncx = xr.open_dataset(ncasdir+'ERA5/era5_'+ctgry+'_seasmean_1950-2014.nc')
era5lat = xr.DataArray(ncx.latitude)
era5lon = xr.DataArray(ncx.longitude)
era5time = xr.DataArray(ncx.time)
era5data = xr.DataArray(getattr(ncx,VAR))
#t2m = xr.DataArray(ncx.t2m)
#tp = xr.DataArray(ncx.tp)
#swvl1 = xr.DataArray(ncx.swvl1)
ncx.close()

era5data.data[0,:,:] = pl.float64('nan')
era5data.data[-1,:,:] = pl.float64('nan')

era5data_jja = era5data[2:-1:4,:,:]
era5data_djf = era5data[4:-1:4,:,:]

#tp.data[0,:,:] = pl.float64('nan')
#t2m.data[-1,:,:] = pl.float64('nan')
#
#tp_jja = tp[2:-13:4,:,:]
#tp_djf = tp[0:-13:4,:,:]

#swvl1.data[0,:,:] = pl.float64('nan')
#swvl1.data[-1,:,:] = pl.float64('nan')
#
#swvl1_jja = swvl1[2:-27:4,:,:]
#swvl1_djf = swvl1[0:-27:4,:,:]

# calculate seasonal GMST
# use HADCRUT4 monthly time series
# anomalies?
#nc_hadcru = xr.open_dataset(ncasdir+'hadcrut5_seasmean.nc')
#hadcru_tas = xr.DataArray(nc_hadcru.tas_mean)
#hadcru_tm = xr.DataArray(nc_hadcru.time)
#nc_hadcru.close()

nc_gmst = xr.open_dataset(era5dir+'era5_t2m_seasmean_1950-2014_gm.nc')
gmst = xr.DataArray(nc_gmst.t2m)
nc_gmst.close()

gmst.data[0] = pl.float64('nan')
gmst.data[-1] = pl.float64('nan')

gmst_jja = gmst.data[2:-1:4]
gmst_djf = gmst.data[4:-1:4]

# smooth GMST with 15 years low-pass filter
# use Pandas rolling average
df1 = pd.DataFrame(pl.squeeze(gmst_jja))
gmst_jja_smth = df1.rolling(15,min_periods=1,center=True).mean()
gmst_jja_smth = pl.squeeze(pl.asarray(gmst_jja_smth))

df2 = pd.DataFrame(pl.squeeze(gmst_djf))
gmst_djf_smth = df2.rolling(15,min_periods=1,center=True).mean()
gmst_djf_smth = pl.squeeze(pl.asarray(gmst_djf_smth))

gmst_smth_all = [gmst_jja_smth,gmst_djf_smth[1:]]


SN_jja, SN_djf = S_N_ratio(era5data_jja,era5data_djf,era5lat,era5lon,gmst_smth_all)


###############################################################################

fig, ax = pl.subplots(1,2,figsize=(12,4.8))
###############################################################################
proj = ccrs.PlateCarree(central_longitude=0)
ext = [-15,42,35,70]

#pl.figure(2)
ax1 = pl.subplot(121,projection=proj,extent=ext)
#ax = pl.axes(projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m',zorder=10)
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        #edgecolor='face',
                                        facecolor='w')
ax1.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax1.contourf(era5lon.data,era5lat.data,SN_jja,transform=proj,cmap='seismic_r',
                 levels=pl.linspace(-2,2,9),extend='both',alpha=0.8)
cn = ax1.contour(era5lon.data,era5lat.data,SN_jja,transform=proj,
                 levels=[-1,0,1,2,3],colors='k',linestyles='solid')
ax1.clabel(cn,manual=True,fmt="%0.0f")#pl.colorbar(cs,orientation='horizontal')
GridLines(ax1,False,True,True,False)
ax1.annotate('JJA ERA5 '+VAR+' S/N',(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

ax2 = pl.subplot(122,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m',zorder=10)
#ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        #edgecolor='face',
#                                        facecolor='w')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax2.contourf(era5lon.data,era5lat.data,SN_djf,transform=proj,cmap='seismic_r',
                 levels=pl.linspace(-2,2,9),extend='both',alpha=0.8)
cn = ax2.contour(era5lon.data,era5lat.data,SN_djf,transform=proj,
                 levels=[-1,0,1,2,3],colors='k',linestyles='solid')
ax2.clabel(cn,manual=True,fmt="%0.0f")#pl.colorbar(cs,orientation='horizontal')
GridLines(ax2,False,True,False,True)
ax2.annotate('DJF ERA5 '+VAR+' S/N',(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)
           
f = pl.gcf()
colax = f.add_axes([0.15,0.10,0.7,0.028])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)#pad=0.05,fraction=0.10,
cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('S/N',fontsize=11)

pl.subplots_adjust(top=1.05,bottom=0.08,left=0.04,right=0.96,hspace=0.12,wspace=0.06)

#pl.savefig(indecis+'figures/S_N_ratio_ERA5_'+VAR+'.png',dpi=375)