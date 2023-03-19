# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 18:36:39 2023

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import pandas as pd
import glob
from scipy import stats
import cartopy.crs as ccrs
import cartopy.util as util
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import matplotlib.path as mplPath
import pcraig_funcs as pc

def S_N_ratio(data,lat,lon,gmst_smth):
    """
    """
    # apply the linear model
    # L(t) = alpha*G(t) + beta
    L = pl.zeros_like(data)
    
    R = pl.zeros_like(data)
    
    pvalue = pl.zeros([lat.size,lon.size])
    
    #for season in range(2):
    for i in range(lat.size):# loop over lon:
        for j in range(lon.size):# loop over lat:
            # linear regression at each grid point
            #GMST = gmst_smth_all[season]
            #data_season = data_all[season]
            out = stats.linregress(gmst_smth,data[:,i,j])
            slope = out[0]
            intercept = out[1]
            pvalue[i,j] = out[3]
            # GMST regressed onto the variable
            
            # calculate L(t) at each grid point
            # alpha & beta are the regresion coefficients
            # G(t) is GMST
            L[:,i,j] = slope*gmst_smth + intercept
            
            # calculate the residuals at each grid point
            # variable minus L(t)
            R[:,i,j] = data[:,i,j] - L[:,i,j]
    
    
    # calculate the signal
    # S = L[timeN] - L[time0]
    signal = L[-1] - L[0]
    #signal_djf = L[1][-1] - L[1][0]
    
    # calculate the noise
    # N = std(residual)
    #noise = pl.std(R,axis=0)
    #noise_djf = pl.std(R[1],axis=0)
    
    # signal-to-noise ratio
    # S/N
    #sig_noi_ratio = signal/noise
    #sig_noi_ratio_djf = signal_djf/noise_djf
    
    #SN_ratios = pl.array([sig_noi_ratio_jja,sig_noi_ratio_djf])
    
    return signal#, pvalue

def AreasCalc2(lon,lat):
    """
    """
    #glat = pl.arange(-89.875,89.876,0.25)
    #glon = pl.arange(-179.875,179.876,0.25)
    
    # Convert lat & lon arrays to radians
    lat_rad = pl.radians(pl.flipud(lat[:]))
    lon_rad = pl.radians(lon[:])
    
    lat_half = pc.HalfGrid(lat_rad)
    nlon = lon_rad.size # number of longitude points
    delta_lambda = (2*pl.pi)/nlon


    #--------------calculate cell areas here, use function from above--------------
    # set up empty array for area, size lat_half X lon_half
    areas = pl.zeros([lon_rad.size,lat_rad.size])
    radius = 6.37*(10**6)
    # loop over latitude and longitude
    for i in range(lon.size): # loops over 256
        for j in range(lat_half.size-1): # loops over 512
            latpair = (lat_half[j+1],lat_half[j])
            areas[i,j] = pc.AreaFullGaussianGrid(radius,delta_lambda,latpair)
    
    #areas_clip = areas[70:326,34:200]
    
    return areas#_clip

def RegionMask(vertices,lon,lat):
    """
    """
    rPath = mplPath.Path(vertices)
    TF = pl.zeros([lon.size,lat.size])
    rmask = pl.zeros([lon.size,lat.size])
    rmask[:] = pl.float32('nan')
    
    for i in range(lon.size):
            for j in range(lat.size):
                X = rPath.contains_point((lon[i],lat[j]))
                TF[i,j] = X
    
    Y = pl.where(TF)
    rmask[Y[0],Y[1]] = 1
    
#    rm0 = pl.zeros_like(rmask)
#    rm0[:120,:] = rmask[120:,:]
#    rm0[120:,:] = rmask[:120,:]
#    rmask = rm0.copy()
#    del rm0
    
    return rmask

def RegionCalc(rmask,lon2,lat2,data,lsmask,areas):
    """
    rmask: region mask (1 inside, 0 outside)
    data: variable to calculate area average of
    lsmask: land-sea mask (1 land, 0 sea)
    areas: each grid cell's area
    """
#    rPath = mplPath.Path(vertices)
#    TF = pl.zeros([lon2.size,lat2.size])
#    rmask = pl.zeros([lon2.size,lat2.size])
#    rmask[:] = pl.float32('nan')
#    
#    for i in range(lon2.size):
#            for j in range(lat2.size):
#                X = rPath.contains_point((lon2[i],lat2[j]))
#                TF[i,j] = X
#    
#    Y = pl.where(TF)
#    rmask[Y[0],Y[1]] = 1
#    
#    rm0 = pl.zeros_like(rmask)
#    rm0[:120,:] = rmask[120:,:]
#    rm0[120:,:] = rmask[:120,:]
#    rmask = rm0.copy()
#    del rm0
    
    #areas = AreasCalc2(lon2,lat2)
    
    rdata = data[:,:]*rmask[:,:] # region data
    rareas = areas*rmask*lsmask # areas of grid cells in region
    
    #Q = pl.ones_like(data) # array of ones, same shape as data
    #f = pl.isnan(data) # where data array has nan values (True/False)
    #d = pl.where(f==True) # indices of nan values in data array
    #Q[d[0],d[1]] = pl.float32('nan') # set grid points in Q to nan
    
    #P = pl.average(rdata[0],weights=pl.nan_to_num(rareas))
    #W = pl.zeros([data.shape[0]])
    #W[0] = pl.float32('nan')
    # sum of product of region data & region areas
    # sum of region areas ()
    W = pl.nansum(rdata*rareas)/(pl.nansum(rareas))
     
    #for i in range(data.shape[0]): # loop over years
        #W[i] = pl.nansum(rdata[i]*rareas)/(pl.nansum(rareas*Q[i]))
    
    return W

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
cmipdir = ncasdir + 'CMIP6/'
maskdir = cmipdir + '/masks/'

nrth_eur = [(-5,48),(-11,51.8),(-11,58),(16.7,71),(28.9,71),(28.9,48)]
sth_eur = [(-11,36),(-11,43.8),(-5,48),(28.9,48),(28.6,41),(25.1,41),(25.1,36),
           (15,36),(10.3,38.1),(2.7,38.1),(-5.6,36)]

var = 'tas'
season = 'JJA'

meanfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_ensmean.nc')
lat = xr.DataArray(meanfile.lat)
lon = xr.DataArray(meanfile.lon)
time = xr.DataArray(meanfile.time)
cmipmean = xr.DataArray(getattr(meanfile,var))
meanfile.close()

areas = AreasCalc2(lon,lat)
nrth_msk = RegionMask(nrth_eur,lon,lat)
sth_msk = RegionMask(sth_eur,lon,lat)
region = 'south'

maskfile = xr.open_dataset(maskdir+'lsmask_cmip6_s1.5.nc')
cmipmask = xr.DataArray(maskfile.topo).data
maskfile.close()

#cm0 = pl.zeros_like(cmipmask) # temporary array
#cm0[:,:120] = cmipmask[:,120:]
#cm0[:,120:] = cmipmask[:,:120]
#cmipmask = cm0.copy()
#del cm0

newlon = lon - 180
if region == 'north':
    msk_in = RegionMask(nrth_eur,newlon,lat)
elif region == 'south':
    msk_in = RegionMask(sth_eur,newlon,lat)

mi0 = pl.zeros_like(msk_in)
mi0[:120,:] = msk_in[120:,:]
mi0[120:,:] = msk_in[:120,:]
msk_in = mi0.copy()
del mi0

allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014_new/*')
tasfiles = glob.glob(cmipdir+'global_means/*')

models = []

modeldata = pl.zeros([len(allfiles),len(time),lat.size,lon.size])
SIGNALS = pl.zeros([len(allfiles),lat.size,lon.size])
SIG_AA = pl.zeros([len(allfiles)])
gmst = pl.zeros([len(tasfiles),len(time)])

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var))
    ncfile.close
    
    tasfile = xr.open_dataset(tasfiles[nci])
    gmst[nci,:] = pl.squeeze(xr.DataArray(tasfile.tas))
    tasfile.close()
    
    #split1 = allfiles[nci].split('/')
    #split2 = split1[-1].split('_')[2]
    models.append(allfiles[nci].split('/')[-1].split('_')[2])

if season == 'JJA':
    data_em_tm = pl.mean(cmipmean[1::4,:,:],axis=0)
    #data_es_tm = pl.mean(cmipsprd[1::4,:,:],axis=0)
    gmst_ts = gmst[:,1::4]
    data_ssn = modeldata[:,1::4]
elif season == 'DJF':
    data_em_tm = pl.mean(cmipmean[3::4,:,:],axis=0)
    #data_es_tm = pl.mean(cmipsprd[3::4,:,:],axis=0)
    gmst_ts = gmst[:,3::4]
    data_ssn = modeldata[:,3::4]

for i in range(len(allfiles)):
    # smooth GMST with 15 years low-pass filter
    # use Pandas rolling average
    df = pd.DataFrame(gmst_ts[i])
    gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
    gmst_smth = pl.squeeze(pl.asarray(gmst_smth))
    
    SIGNALS[i] = S_N_ratio(data_ssn[i],lat,lon,gmst_smth)
    #cmipmask[:,:] = 1
    SIG_AA[i] = RegionCalc(msk_in,lon,lat,SIGNALS[i].T,cmipmask.T,areas)

out = pd.DataFrame(data=SIG_AA,index=models,columns=['S('+season+','+region+')'])

#with pd.ExcelWriter('signals_revisions.xlsx',mode='a') as writer:
#    out.to_excel(writer, sheet_name=season+'_'+region)
#out.to_excel(jashome+'signals_revisions.xlsx',sheet_name=season+'_'+region)