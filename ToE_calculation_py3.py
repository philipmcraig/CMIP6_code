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
import cf
import cfplot as cfp
import pcraig_funcs as pc

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
eobsdir = ncasdir + 'EOBS/'

inst = 'MOHC'
model = 'UKESM1-0-LL'
variant = 'r1i1p1f2'
cmipfiles = glob.glob(cmipdir+'/'+inst+'/'+model+'/'+variant+'/*')

VAR = 'mrso'
fname = [s for s in cmipfiles if 'mrso' in s][0]

# calculate seasonal GMST
# use HADCRUT4 monthly time series
# anomalies?
#nc_hadcru = xr.open_dataset(ncasdir+'hadcrut5_seasmean.nc')
#hadcru_tas = xr.DataArray(nc_hadcru.tas_mean)
#hadcru_tm = xr.DataArray(nc_hadcru.time)
#nc_hadcru.close()

# reading with cf-python
nc_hadcru = cf.read(ncasdir+'hadcrut5_seasmean.nc')
hadcru_tas = nc_hadcru.select('long_name=blended air_temperature_anomaly over land with sea_water_temperature_anomaly')[0]
#hadcru_tm = nc_hadcru.select('time')

#print(hadcru_tas)

#hadcru_tas[0] = pl.float64('nan')
#hadcru_tas[-1] = pl.float64('nan')

hadcru_jja = hadcru_tas[2:-1:4]
hadcru_djf = hadcru_tas[0:-1:4]

#print(hadcru_jja)

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
df1 = pd.DataFrame(hadcru_jja.data)
gmst_jja_smth = df1.rolling(15,min_periods=1).mean()
gmst_jja_smth = pl.squeeze(pl.asarray(gmst_jja_smth))

df2 = pd.DataFrame(hadcru_djf.data)
gmst_djf_smth = df2.rolling(15,min_periods=1).mean()
gmst_djf_smth = pl.squeeze(pl.asarray(gmst_djf_smth))

pl.plot(hadcru_jja[100:165],label='HADCRUT5 JJA GMST anomalies')
pl.plot(gmst_jja_smth[100:165],label='15 year moving average')

pl.ylabel(hadcru_tas.units)
#pl.savefig(indecis+'figures/gmst_smooth.png',dpi=350)

# open up variable/index & calculate seasonal anomalies
nc_cmip = xr.open_dataset(fname)
cmip_lon = xr.DataArray(nc_cmip.lon)
cmip_lat = xr.DataArray(nc_cmip.lat)
#cmip_data = xr.DataArray(nc_cmip.mrso)
#cmip_time = xr.DataArray(nc_cmip.time)
nc_cmip.close()

nc_cmip = cf.read(fname)[0]
cmip_data = nc_cmip.array
#cmip_data = nc_cmip.select('long_name=mass_content_of_water_in_soil')

#cmip_lon.data[0] = 0.0
#cmip_lon.data[-1] = 360.0
#cmip_data.data[0,:,:] = pl.float64('nan')
#cmip_data.data[-1,:,:] = pl.float64('nan')

cmip_jja = cmip_data[2:-1:4]
cmip_djf = cmip_data[0:-1:4]

#print(cmip_djf)

# apply the linear model
# L(t) = alpha*G(t) + beta
L = pl.zeros_like(cmip_djf[1:])
# residuals array
R = pl.zeros_like(cmip_djf[1:])

for i in range(144):# loop over lon:
    for j in range(192):# loop over lat:
        # linear regression at each grid point
        out = stats.linregress(gmst_djf_smth[101:165],cmip_djf[1:,i,j])
        slope = out[0]
        intercept = out[1]
        # GMST regressed onto the variable
        
        # calculate L(t) at each grid point
        # alpha & beta are the regresion coefficients
        # G(t) is GMST
        L[:,i,j] = slope*gmst_djf_smth[101:165] + intercept
        
        # calculate the residuals at each grid point
        # variable minus L(t)
        R[:,i,j] = cmip_djf[1:,i,j] - L[:,i,j]


# calculate the signal
# S = L[timeN] - L[time0]
signal = L[-1] - L[0]

# calculate the noise
# N = std(residual)
noise = pl.std(R,axis=0)

# signal-to-noise ratio
# S/N
sig_noi_ratio = signal/noise


###############################################################################
#proj = ccrs.PlateCarree(central_longitude=0)
#ext = [-15,42,35,70]

#pl.figure(2)
#ax = pl.axes(projection=proj,extent=ext)
#ax.coastlines(linewidth=0.5,resolution='50m',zorder=10)
#ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        #edgecolor='face',
 #                                       facecolor='w')
#ax.add_feature(ocean_50m,alpha=1,zorder=5)

#cs = ax.contourf(cmip_lon,cmip_lat,sig_noi_ratio,transform=proj,cmap='BrBG',
 #                levels=pl.linspace(-1,1,11),extend='both')
#ax.contour(cmip_lon,cmip_lat,sig_noi_ratio,transform=proj,levels=[-1,0,1],
 #          colors='k',linestyles=['dashed','solid','dotted'])
#pl.colorbar(cs)

cfp.mapset(-15,42,35,70)
cfp.cscale('precip4_diff_19lev')
cfp.levs(min=-1,max=1,step=0.2)
cfp.con(f=sig_noi_ratio,x=cmip_lon.data,y=cmip_lat.data,ptype=1)
#pl.show()