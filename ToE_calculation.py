# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 13:22:53 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import pandas as pd
import glob
from scipy import stats
from scipy import signal
import cartopy.crs as ccrs
import cartopy.util as util
import cartopy.io.shapereader as shpreader
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
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
eobsdir = ncasdir + 'EOBS/'

inst = 'EC-Earth-Consortium'
model = 'EC-Earth3'
variant = 'r1i1p1f1'
table = 'Amon'
#cmipfiles = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/*')

VAR = 'pr'
#fname = [s for s in cmipfiles if VAR in s][0]
fname = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/'+VAR+'_'+table+'*')[0]

splits = fname.split('/')
gridtype = splits[-1].split('_')[-2]

# calculate seasonal GMST
# use HADCRUT4 monthly time series
# anomalies?
#nc_hadcru = xr.open_dataset(ncasdir+'hadcrut5_seasmean.nc')
#hadcru_tas = xr.DataArray(nc_hadcru.tas_mean)
#hadcru_tm = xr.DataArray(nc_hadcru.time)
#nc_hadcru.close()
#
#hadcru_tas.data[0] = pl.float64('nan')
#hadcru_tas.data[-1] = pl.float64('nan')
#
#hadcru_jja = hadcru_tas.data[2:-1:4]
#hadcru_djf = hadcru_tas.data[0:-1:4]

#dates = pl.asarray([str(pd.to_datetime(i))[:7] for i in hadcru_tm.data])
#start = pl.where(dates=='1950-01')[0][0]
#end = pl.where(dates=='2017-12')[0][0]
#dates = dates[start+2:end+2]
#
#hctas_1950_2018 = hadcru_tas.data[start+2:end+2]
#hctas_1950_2018 = pl.reshape(hctas_1950_2018,newshape=(68,12))
#hctas_seasmeans = 

# smooth GMST with 15 years low-pass filter
# use Pandas rolling average
#df1 = pd.DataFrame(hadcru_jja)
#gmst_jja_smth = df1.rolling(15,min_periods=1,center=True).mean()
#gmst_jja_smth = pl.squeeze(pl.asarray(gmst_jja_smth))
#
#df2 = pd.DataFrame(hadcru_djf)
#gmst_djf_smth = df2.rolling(15,min_periods=1,center=True).mean()
#gmst_djf_smth = pl.squeeze(pl.asarray(gmst_djf_smth))
#
#gmst_smth_all = [gmst_jja_smth[100:165],gmst_djf_smth[101:165]]


#fig, ax = pl.subplots(2,2,figsize=(10,7))

#ax1 = pl.subplot(221)
#ax1.plot(pl.linspace(1950,2014,65),hadcru_jja[100:165],label='HADCRUT5 JJA GMST anomalies')
#ax1.plot(pl.linspace(1950,2014,65),gmst_jja_smth[100:165],label='15 year moving average')
#pl.legend(loc=2)
#pl.xlim(1950,2014); pl.ylim(-0.4,0.8)
#ax1.grid(axis='y',ls='--',color='grey')
#pl.ylabel(hadcru_tas.units,labelpad=-5,fontsize=14)
##pl.savefig(indecis+'figures/gmst_smooth.png',dpi=350)
#
#ax2 = pl.subplot(222)
#ax2.plot(pl.linspace(1950,2014,65),hadcru_djf[100:165],label='HADCRUT5 DJF GMST anomalies')
#ax2.plot(pl.linspace(1950,2014,65),gmst_djf_smth[100:165],label='15 year moving average')
#pl.legend(loc=2)
#pl.xlim(1950,2014); pl.ylim(-0.4,0.8)
#ax2.grid(axis='y',ls='--',color='grey')
#ax2.yaxis.set_major_formatter(pl.NullFormatter())

nc_gmst = xr.open_dataset(cmipdir+'/'+inst+'/'+model+'/'+variant+'/tas_'+\
                    model+'_historical_'+variant+'_'+gridtype+'_seasmean_gm.nc')
gmst_ts = xr.DataArray(nc_gmst.tas)
gmst_time = xr.DataArray(nc_gmst.time)
nc_gmst.close()

gmst_ts = pl.squeeze(gmst_ts.data)

year0 = int(str(gmst_time.data[0])[:4])

if year0 == 1950:# and season == 'JJA':
    start_jja = 2 # slice array from JJA 1950 (index 2 of time dimension)
    start_djf = 0
elif year0 < 1950:# and season == 'JJA':
    start_jja = (1950-year0)*4 + 2
    start_djf = (1950-year0)*4
#elif year0 == 1950 and season == 'DJF':
#    start = 4
#elif year0 < 1950 and season == 'DJF':
#    start = (1950-year0)*4 + 4

gmst_djf = gmst_ts[start_jja:-1:4]
gmst_jja = gmst_ts[start_djf:-1:4]

# open up variable/index & calculate seasonal anomalies
nc_cmip = xr.open_dataset(fname)
cmip_lon = xr.DataArray(nc_cmip.lon)
cmip_lat = xr.DataArray(nc_cmip.lat)
#cmip_temp = xr.DataArray(nc_cmip.tas)
cmip_data = xr.DataArray(getattr(nc_cmip,VAR))
cmip_time = xr.DataArray(nc_cmip.time)
nc_cmip.close()

#temp_jja = cmip_temp.data[2:-1:4,:,:]
#temp_djf = cmip_temp.data[0:-1:4,:,:]

########### ------ CALCULATE GLOBALLY AVERAGED TEMPERATURE ----- ##############

#lat_rad = pl.radians(cmip_lat.data[:])
#lon_rad = pl.radians(cmip_lon.data[:])
#
#lat_half = pc.HalfGrid(pl.flipud(lat_rad))
#
#
#nlon = lon_rad.shape[0] # number of longitude points
#delta_lambda = (2*pl.pi)/nlon
#
#areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
#radius = 6.37*(10**6)
## loop over latitude and longitude
#for i in range(lat_half.shape[0]-1): # loops over 256
#    latpair = (lat_half[i],lat_half[i+1])
#    for j in range(cmip_lon.shape[0]): # loops over 512
#        #lonpair = (lon_half[i],lon_half[i+1])
#        areas[i,j] = pc.AreaFullGaussianGrid(radius,delta_lambda,latpair)
#
#gmst = pl.zeros([cmip_temp.shape[0]])
#for i in range(cmip_temp.shape[0]):
#    gmst[i] = pl.average(cmip_temp.data[i],weights=areas)

###############################################################################

#cmip_lon.data[0] = 0.0
#cmip_lon.data[-1] = 360.0
cmip_data.data[0,:,:] = pl.float64('nan')
cmip_data.data[-1,:,:] = pl.float64('nan')

cmip_jja = cmip_data.data[start_jja:-1:4,:,:]
cmip_djf = cmip_data.data[start_djf:-1:4,:,:]
cmip_all = [cmip_jja,cmip_djf[1:]]

# smooth GMST with 15 years low-pass filter
# use Pandas rolling average
df1 = pd.DataFrame(gmst_jja)
gmst_jja_smth = df1.rolling(15,min_periods=1,center=True).mean()
gmst_jja_smth = pl.squeeze(pl.asarray(gmst_jja_smth))

df2 = pd.DataFrame(gmst_djf)
gmst_djf_smth = df2.rolling(15,min_periods=1,center=True).mean()
gmst_djf_smth = pl.squeeze(pl.asarray(gmst_djf_smth))

gmst_smth_all = [gmst_jja_smth,gmst_djf_smth[1:]]

# apply the linear model
# L(t) = alpha*G(t) + beta
L = [pl.zeros_like(cmip_jja[:,:,:]),pl.zeros_like(cmip_djf[1:,:,:])]
# residuals array
R = [pl.zeros_like(cmip_jja[:,:,:]),pl.zeros_like(cmip_djf[1:,:,:])]

for season in range(2):
    for i in range(cmip_lat.size):# loop over lon:
        for j in range(cmip_lon.size):# loop over lat:
            # linear regression at each grid point
            GMST = gmst_smth_all[season]
            cmip_season = cmip_all[season]
            out = stats.linregress(GMST,cmip_season[:,i,j])
            slope = out[0]
            intercept = out[1]
            # GMST regressed onto the variable
            
            # calculate L(t) at each grid point
            # alpha & beta are the regresion coefficients
            # G(t) is GMST
            L[season][:,i,j] = slope*GMST + intercept
            
            # calculate the residuals at each grid point
            # variable minus L(t)
            R[season][:,i,j] = cmip_season[:,i,j] - L[season][:,i,j]


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

#cyclic_snr_jja, cyclic_lons = util.add_cyclic_point(sig_noi_ratio_jja, coord=cmip_lon)
#cyclic_snr_djf = util.add_cyclic_point(sig_noi_ratio_djf)
lonshift = cmip_lon.data - 180

snr_jja_shift = pl.zeros_like(sig_noi_ratio_jja)
snr_jja_shift[:,:cmip_lon.size/2] = sig_noi_ratio_jja[:,cmip_lon.size/2:]
snr_jja_shift[:,cmip_lon.size/2:] = sig_noi_ratio_jja[:,:cmip_lon.size/2]

snr_djf_shift = pl.zeros_like(sig_noi_ratio_djf)
snr_djf_shift[:,:cmip_lon.size/2] = sig_noi_ratio_djf[:,cmip_lon.size/2:]
snr_djf_shift[:,cmip_lon.size/2:] = sig_noi_ratio_djf[:,:cmip_lon.size/2]

#sig_noi_ratio = pl.array([snr_jja_shift,snr_djf_shift])

fig, ax = pl.subplots(1,2,figsize=(10,4))
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

cs = ax1.contourf(lonshift,cmip_lat,snr_jja_shift,transform=proj,cmap='seismic',
                 levels=pl.linspace(-2,2,9),extend='both',alpha=0.8)
cn = ax1.contour(lonshift,cmip_lat,snr_jja_shift,transform=proj,
                 levels=[-1,0,1,2,3],colors='k',linestyles='solid')
ax1.clabel(cn,manual=True,fmt="%0.0f")#pl.colorbar(cs,orientation='horizontal')
GridLines(ax1,False,True,True,False)
ax1.annotate('JJA '+model+' '+VAR+' S/N',(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)

ax2 = pl.subplot(122,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m',zorder=10)
#ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        #edgecolor='face',
#                                        facecolor='w')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax2.contourf(lonshift,cmip_lat,snr_djf_shift,transform=proj,cmap='seismic',
                 levels=pl.linspace(-2,2,9),extend='both',alpha=0.8)
cn = ax2.contour(lonshift,cmip_lat,snr_djf_shift,transform=proj,
                 levels=[-1,0,1,2,3],colors='k',linestyles='solid')
ax2.clabel(cn,manual=True,fmt="%0.0f")#pl.colorbar(cs,orientation='horizontal')
GridLines(ax2,False,True,False,True)
ax2.annotate('DJF '+model+' '+VAR+' S/N',(-13.8,69.3),bbox={'facecolor':'w'},
             xycoords='data',zorder=50)
           
f = pl.gcf()
colax = f.add_axes([0.15,0.12,0.7,0.025])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)#pad=0.05,fraction=0.10,
cb.set_ticks(pl.linspace(-2,2,9))
cb.set_label('S/N',fontsize=11)

pl.subplots_adjust(top=1.00,bottom=0.13,left=0.05,right=0.95,hspace=0.12,wspace=0.06)

pl.savefig(indecis+'figures/S_N_ratio_'+model+'_'+VAR+'.png')