# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 16:16:54 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import pandas as pd
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
import glob
from adjustText import adjust_text
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mplPath
import pcraig_funcs as pc

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
    TF = pl.zeros([lat.size,lon.size])
    rmask = pl.zeros([lat.size,lon.size])
    rmask[:] = pl.float32('nan')
    
    for i in range(lat.size):
            for j in range(lon.size):
                X = rPath.contains_point((lat[i],lon[j]))
                TF[i,j] = X
    
    Y = pl.where(TF)
    rmask[Y[0],Y[1]] = 1
    
#    rm0 = pl.zeros_like(rmask)
#    rm0[:,:lon.size/2] = rmask[:,lon.size/2:]
#    rm0[:,lon.size/2:] = rmask[:,:lon.size/2]
#    rmask = rm0.copy()
#    del rm0
    
    return rmask

def RegionCalc(rmask,lon2,lat2,data,lsmask,areas):
    """
    """
    
    rdata = data[:,:]*rmask[:,:]#None,
    rareas = areas*rmask*lsmask
    
    Q = pl.ones_like(data)
    f = pl.isnan(data)
    d = pl.where(f==True)
    Q[d[0],d[1]] = pl.float32('nan') #,d[2]
    
    #P = pl.average(rdata[0],weights=pl.nan_to_num(rareas))
    #W = pl.zeros([data.shape[0]])
    #W[0] = pl.float32('nan')
    W = pl.nansum(rdata*rareas)/(pl.nansum(rareas*Q))
     
    #for i in range(data.shape[0]): # loop over years
        #W[i] = pl.nansum(rdata[i]*rareas)/(pl.nansum(rareas*Q[i]))
    
    return W

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'

inst = 'EC-Earth-Consortium'
model = 'EC-Earth3'
variant = 'r1i1p1f1'
table = 'Amon'
region = 'north'

VAR = 'tas'
fname = glob.glob(cmipdir+inst+'/'+model+'/'+variant+'/'+VAR+'_'+table+'*')[0]

splits = fname.split('/')
gridtype = splits[-1].split('_')[-2]

nrth_eur = [(48,-5),(51.8,-11),(58,-11,),(71,16.7),(71,28.9),(48,28.9)]
sth_eur = [(-11,36),(-11,43.8),(-5,48),(28.9,48),(28.6,41),(25.1,41),(25.1,36),
           (15,36),(10.3,38.1),(2.7,38.1),(-5.6,36)]



#nc_gmst = xr.open_dataset(cmipdir+'/'+inst+'/'+model+'/'+variant+'/tas_'+\
#                    model+'_historical_'+variant+'_'+gridtype+'_seasmean_gm.nc')
#gmst_ts = xr.DataArray(nc_gmst.tas)
#gmst_time = xr.DataArray(nc_gmst.time)
#nc_gmst.close()
#
#gmst_ts = pl.squeeze(gmst_ts.data)
#
#ensmean_nc = xr.open_dataset(cmipdir+'standard_grid/tas_cmip6_ensmean_gm.nc')
#gmst_em = xr.DataArray(ensmean_nc.tas)
#ensmean_nc.close()
#
#gmst_em = pl.squeeze(gmst_em.data)
#
#era5_nc = xr.open_dataset(era5dir+'era5_t2m_seasmean_1950-2017_gm.nc')
#gmst_e5 = xr.DataArray(era5_nc.t2m)
#e5_time = xr.DataArray(era5_nc.time)
#era5_nc.close()
#
#gmst_e5 = pl.squeeze(gmst_e5.data)
#
#pl.plot(gmst_ts[4:-1:4],label=model)
#pl.plot(gmst_em[4:-1:4],label='CMIP6 ensemble mean')
#pl.plot(gmst_e5[4:-1:4],label='ERA5')
#pl.legend()
#
#pl.tight_layout()

ncfile = xr.open_dataset(fname)
model_tas = xr.DataArray(ncfile.tas)
lon = xr.DataArray(ncfile.lon)
lat = xr.DataArray(ncfile.lat)
ncfile.close()

#newlon = pl.zeros_like(lon.data)
#newlon[:lon.size/2] = lon.data[lon.size/2:] - 360.
#newlon[lon.size/2:] = lon.data[:lon.size/2]
newlon = lon.data - 180.

mt0 = pl.zeros_like(model_tas.data)
mt0[:,:,:lon.size/2] = model_tas.data[:,:,lon.size/2:]
mt0[:,:,lon.size/2:] = model_tas.data[:,:,:lon.size/2]
model_tas = mt0.copy()
del mt0

if region == 'north':
    msk_in = RegionMask(nrth_eur,newlon,lat)
    title_lab = 'North'
    filelab = 'nrth'
elif region == 'south':
    msk_in = RegionMask(sth_eur,lon,lat)
    title_lab = 'South'
    filelab = 'sth'

