# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 14:34:36 2021

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
    
    return signal

def AreasCalc(lons,lats,north,south,west,east):
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
    for i in range(lons.size): # loops over 256
        for j in range(lat_half.size-1): # loops over 512
            latpair = (lat_half[j+1],lat_half[j])
            areas[i,j] = pc.AreaFullGaussianGrid(radius,delta_lambda,latpair)
    
    areas_clip = areas[west:east+1,south:north+1].T
    
    return areas_clip

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

def season_outs(season,gmst,era5_gmst,eobs_gmst,modeldata,era5data,eobsdata):
    """
    """
    if season == 'JJA':
    #    data_em_tm = pl.mean(cmipmean[1:-1:4,:,:],axis=0)#/(0.1*1000)#*86400#
    #    data_es_tm = pl.mean(cmipsprd[1::4,:,:],axis=0)#/(0.1*1000)#*86400
    #    eobs_mn = pl.mean(eobsdata[1:-26:4],axis=0)# + 273.15
    #    era5_mn = pl.mean(era5data[2:-1:4],axis=0)*1000
        gmst_ts = gmst[:,1:-1:4]
        e5_gmst = era5_gmst[2:-1:4]
        eb_gmst = eobs_gmst[402:-26:4] + 273.15
        data_ssn = modeldata[:,1:-1:4]*86400#/(0.1*1000)##
        #c6tas_ssn = modeltemp[:,1:-1:4]
    #    c6pr_ssn = modelprec[:,1:-1:4]*86400
    #    c6sm_ssn = modelsomo[:,1:-1:4]/(0.1*1000)
        era5_ssn = era5data[2:-14:4]*1000
        eobs_ssn = eobsdata[1:-1:4,:,:]# + 273.15
        #e5tmp_ssn = era5temp[2:-14:4]
        #ebtg_ssn = eobstemp[1:-1:4,:,:] + 273.15
    #    e5tp_ssn = era5prec[2:-1:4]*1000
    #    ebrr_ssn = eobsprec[1:-26:4]
    #    e5sm_ssn = era5somo[2:-1:4]
    elif season == 'DJF':
    #    data_em_tm = pl.mean(cmipmean[3:-1:4,:,:],axis=0)/(0.1*1000)#*86400#
    #    data_es_tm = pl.mean(cmipsprd[3::4,:,:],axis=0)#/(0.1*1000)#*86400
    #    eobs_mn = pl.mean(eobsdata[3:-26:4],axis=0)# + 273.15
    #    era5_mn = pl.mean(era5data[4:-14:4,:,:],axis=0)*1000
        gmst_ts = gmst[:,3:-1:4]
        e5_gmst = era5_gmst[4:-1:4]
        eb_gmst = eobs_gmst[404:-26:4] + 273.15
        data_ssn = modeldata[:,3:-1:4]*86400#/(0.1*1000)##
        #c6tas_ssn = modeltemp[:,3:-1:4]
    #    c6pr_ssn = modelprec[:,3:-1:4]*86400
    #    c6sm_ssn = modelsomo[:,3:-1:4]/(0.1*1000)
        era5_ssn = era5data[4:-14:4]*1000
        eobs_ssn = eobsdata[3:-1:4,:,:]# + 273.15
        #e5tmp_ssn = era5temp[4:-14:4]
        #ebtg_ssn = eobstemp[3:-1:4,:,:] + 273.15
    #    e5tp_ssn = era5prec[4:-1:4]*1000
    #    ebrr_ssn = eobsprec[3:-26:4]
     #   e5sm_ssn = era5somo[4:-1:4]

    return gmst_ts, e5_gmst, eb_gmst, data_ssn, era5_ssn, eobs_ssn

def season_outs2(season,gmst,era5_gmst,modeldata,era5data):
    """
    """
    if season == 'JJA':
        gmst_ts = gmst[:,1:-1:4]
        e5_gmst = era5_gmst[2:-1:4]
        data_ssn = modeldata[:,1:-1:4]/(0.1*1000)
        era5_ssn = era5data[2:-1:4]
    elif season == 'DJF':
        gmst_ts = gmst[:,3:-1:4]
        e5_gmst = era5_gmst[4:-1:4]
        data_ssn = modeldata[:,3:-1:4]/(0.1*1000)
        era5_ssn = era5data[4:-1:4]
    
    return gmst_ts, e5_gmst, data_ssn, era5_ssn

def region_outs(region,nrth_eur,sth_eur,newlon,lat):
    """
    """
    if region == 'north':
        msk_in = RegionMask(nrth_eur,newlon,lat)
        title_lab = 'North'
        #filelab = 'nrth'
    elif region == 'south':
        msk_in = RegionMask(sth_eur,newlon,lat)
        title_lab = 'South'
        #filelab = 'sth'
    
    return msk_in, title_lab

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
cmipdir = ncasdir + 'CMIP6/'
era5dir = ncasdir + 'ERA5/'
eobsdir = ncasdir + 'EOBS/'
indecis = ncasdir + 'INDECIS/'
maskdir = cmipdir + '/masks/'

var = 'mrsos'
#var2 = 'rr'
var3 = 'swvl1'
season = ['DJF','DJF']#,'JJA','DJF']#,'JJA','JJA']
region = ['north','south']#,'north','north']#,'south','south']

meanfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_ensmean.nc')
lat = xr.DataArray(meanfile.lat[60:]).data
lon = xr.DataArray(meanfile.lon).data
time = xr.DataArray(meanfile.time).data
cmipmean = xr.DataArray(getattr(meanfile,var)[:,60:,:]).data
meanfile.close()

maskfile = xr.open_dataset(maskdir+'lsmask_cmip6_s1.5.nc')
cmipmask = xr.DataArray(maskfile.topo[60:,:]).data
maskfile.close()

cm0 = pl.zeros_like(cmipmask) # temporary array
cm0[:,:120] = cmipmask[:,120:]
cm0[:,120:] = cmipmask[:,:120]
cmipmask = cm0.copy()
del cm0

allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014_new/*')
tempfiles = glob.glob(cmipdir+'standard_grid/tas_1950-2014_new/*')
tasfiles = glob.glob(cmipdir+'global_means/*')

precfiles = glob.glob(cmipdir+'standard_grid/pr_1950-2014/*')
somofiles = glob.glob(cmipdir+'standard_grid/mrsos_1950-2014/*')

modeldata = pl.zeros([len(allfiles),len(time),lat.size,lon.size])
#modeltemp = modeldata.copy()
#modelprec = modeldata.copy()
#modelsomo = modeldata.copy()
modnames = []
nums = pl.linspace(1,17,17).astype(int).astype(str)
gmst = pl.zeros([len(tasfiles),len(time)])

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var)[:,60:,:]).data
    ncfile.close
    
#    precfile = xr.open_dataset(precfiles[nci])
#    modelprec[nci,:,:,:] = xr.DataArray(precfile.pr)
#    precfile.close()
    
#    tempfile = xr.open_dataset(tempfiles[nci])
#    modeltemp[nci,:,:,:] = xr.DataArray(tempfile.tas[:,60:,:]).data
#    tempfile.close()
#    
#    somofile = xr.open_dataset(somofiles[nci])
#    modelsomo[nci,:,:,:] = xr.DataArray(somofile.mrsos[:,60:,:]).data
#    somofile.close()
    
    split1 = allfiles[nci].split('/')
    split2 = split1[-1].split('_')[2]
    modnames.append(split2)
    
    tasfile = xr.open_dataset(tasfiles[nci])
    gmst[nci,:] = pl.squeeze(xr.DataArray(tasfile.tas)).data
    tasfile.close()

modnames[1] = modnames[1][:-2]
modnames[2] = modnames[2][:7]
modnames[6] = modnames[6][:-2]
modnames[9] = modnames[9][:-4]
modnames[12] = modnames[12][:-5]
modnames[13] = modnames[13][:-2]
modnames[14] = modnames[14][:-3]
modnames[16] = modnames[16][:6]

era5file = xr.open_dataset(era5dir+'era5_land_seasmean_s1.5.nc')
era5data = xr.DataArray(getattr(era5file,var3)[:,60:,:]).data
#era5temp = xr.DataArray(era5file.t2m[:,60:,:])
#era5prec = xr.DataArray(era5file.tp)
#era5somo = xr.DataArray(era5file.swvl1[:,60:,:]).data
era5lon = xr.DataArray(era5file.lon).data
era5lat = xr.DataArray(era5file.lat[60:]).data
era5time = xr.DataArray(era5file.time).data
era5file.close()

#eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobsdata = xr.DataArray(getattr(eobsfile,var2)[:-26,60:,:]).data
##eobstemp = xr.DataArray(eobsfile.tg[:,60:,:]).data
#eobstime = xr.DataArray(eobsfile.time[:-26]).data
#eobslat = xr.DataArray(eobsfile.lat[60:]).data
#eobslon = xr.DataArray(eobsfile.lon).data
#eobsfile.close()

#eb_tmpfile = xr.open_dataset(eobsdir+'tg_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobstemp = xr.DataArray(eb_tmpfile.tg)
#eb_tmpfile.close()

#eb_prcfile = xr.open_dataset(eobsdir+'rr_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobsprec = xr.DataArray(eb_prcfile.rr)
#eb_prcfile.close()

nc_gmst = xr.open_dataset(era5dir+'era5_t2m_seasmean_1950-2014_gm.nc')
era5_gmst = xr.DataArray(nc_gmst.t2m).data
e5_longtime = xr.DataArray(nc_gmst.time).data
nc_gmst.close()

nc_gmst2 = xr.open_dataset(ncasdir+'hadcrut5_seasmean.nc')
eobs_gmst = xr.DataArray(nc_gmst2.tas_mean).data
gmst_time = xr.DataArray(nc_gmst2.time).data
nc_gmst2.close()

#if season == 'JJA':
##    data_em_tm = pl.mean(cmipmean[1:-1:4,:,:],axis=0)#/(0.1*1000)#*86400#
##    data_es_tm = pl.mean(cmipsprd[1::4,:,:],axis=0)#/(0.1*1000)#*86400
##    eobs_mn = pl.mean(eobsdata[1:-26:4],axis=0)# + 273.15
##    era5_mn = pl.mean(era5data[2:-1:4],axis=0)*1000
#    gmst_ts = gmst[:,1:-1:4]
#    e5_gmst = era5_gmst[2:-1:4]
##    eb_gmst = eobs_gmst[401:-26:4] + 273.15
#    data_ssn = modeldata[:,1:-1:4]#/(0.1*1000)#*86400#
##    c6tas_ssn = modeltemp[:,1:-1:4]
##    c6pr_ssn = modelprec[:,1:-1:4]*86400
#    c6sm_ssn = modelsomo[:,1:-1:4]/(0.1*1000)
#    era5_ssn = era5data[2:-1:4]#*1000
##    eobs_ssn = eobsdata[1:-26:4,:,:] + 273.15
##    e5tmp_ssn = era5temp[2:-14:4]
##    ebtg_ssn = eobstemp[1:-1:4,:,:] + 273.15
##    e5tp_ssn = era5prec[2:-1:4]*1000
##    ebrr_ssn = eobsprec[1:-26:4]
#    e5sm_ssn = era5somo[2:-1:4]
#elif season == 'DJF':
##    data_em_tm = pl.mean(cmipmean[3:-1:4,:,:],axis=0)/(0.1*1000)#*86400#
##    data_es_tm = pl.mean(cmipsprd[3::4,:,:],axis=0)#/(0.1*1000)#*86400
##    eobs_mn = pl.mean(eobsdata[3:-26:4],axis=0)# + 273.15
##    era5_mn = pl.mean(era5data[4:-14:4,:,:],axis=0)*1000
#    gmst_ts = gmst[:,3:-1:4]
#    e5_gmst = era5_gmst[4:-1:4]
##    eb_gmst = eobs_gmst[403:-26:4] + 273.15
#    data_ssn = modeldata[:,3:-1:4]#/(0.1*1000)#*86400#
# #   c6tas_ssn = modeltemp[:,3:-1:4]
##    c6pr_ssn = modelprec[:,3:-1:4]*86400
#    c6sm_ssn = modelsomo[:,3:-1:4]/(0.1*1000)
#    era5_ssn = era5data[4:-1:4]#*1000
##    eobs_ssn = eobsdata[3:-26:4,:,:] + 273.15
##    e5tmp_ssn = era5temp[4:-14:4]
##    ebtg_ssn = eobstemp[3:-1:4,:,:] + 273.15
##    e5tp_ssn = era5prec[4:-1:4]*1000
##    ebrr_ssn = eobsprec[3:-26:4]
#    e5sm_ssn = era5somo[4:-1:4]

#BIAS_mn = data_em_tm - era5_mn

#BIAS_mods = pl.zeros([len(allfiles),lat.size,lon.size])
#for i in range(len(allfiles)):
#    BIAS_mods[i] = pl.mean(data_ssn[i],axis=0) - era5_ssn[i]

###############################################################################
#SIGNALS = pl.zeros([len(allfiles),lat.size,lon.size])
#pvals = SN_ratios.copy()

#for i in range(len(allfiles)):
#    # smooth GMST with 15 years low-pass filter
#    # use Pandas rolling average
#    df = pd.DataFrame(gmst_ts[i])
#    gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
#    gmst_smth = pl.squeeze(pl.asarray(gmst_smth))
#    
#    SIGNALS[i] = S_N_ratio(c6sm_ssn[i],lat,lon,gmst_smth)
#
#signals_mean = pl.nanmean(SIGNALS,axis=0)
###############################################################################

###############################################################################
#df = pd.DataFrame(pl.squeeze(e5_gmst))
#e5gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
#e5gmst_smth = pl.squeeze(pl.asarray(e5gmst_smth))

#df = pd.DataFrame(pl.squeeze(eb_gmst))
#ebgmst_smth = df.rolling(15,min_periods=1,center=True).mean()
#ebgmst_smth = pl.squeeze(pl.asarray(ebgmst_smth))

#E5SIG = S_N_ratio(e5sm_ssn,era5lat,era5lon,e5gmst_smth)
#EBSIG = S_N_ratio(ebrr_ssn,eobslat.data,eobslon.data,ebgmst_smth)

###############################################################################

# rearrange the arrays to make slicing and averaging easier:
#sm0 = pl.zeros_like(signals_mean) # temporary array
#sm0[:,:120] = signals_mean[:,120:]
#sm0[:,120:] = signals_mean[:,:120]
#signals_mean = sm0.copy()
#del sm0
#
#e50 = pl.zeros_like(E5SIG)
#e50[:,:120] = E5SIG[:,120:]
#e50[:,120:] = E5SIG[:,:120]
#E5SIG = e50.copy()
#del e50

newlon = lon - 180

#bm0 = pl.zeros_like(BIAS_mods)
#bm0[:,:,:120] = BIAS_mods[:,:,120:]
#bm0[:,:,120:] = BIAS_mods[:,:,:120]
#BIAS_mods = bm0.copy()
#del bm0
#
#sg0 = pl.zeros_like(SIGNALS)
#sg0[:,:,:120] = SIGNALS[:,:,120:]
#sg0[:,:,120:] = SIGNALS[:,:,:120]
#SIGNALS = sg0.copy()
#del sg0
###############################################################################

nrth_eur = [(-5,48),(-11,51.8),(-11,58),(16.7,71),(28.9,71),(28.9,48)]
sth_eur = [(-11,36),(-11,43.8),(-5,48),(28.9,48),(28.6,41),(25.1,41),(25.1,36),
           (15,36),(10.3,38.1),(2.7,38.1),(-5.6,36)]

# calculate area-averages of model signal, observed signal & temperature bias

areas = AreasCalc2(lon,lat)
#nrth_msk = RegionMask(nrth_eur,lon,lat)
#sth_msk = RegionMask(sth_eur,lon,lat)

meanbias = pl.zeros([len(season)])
corrs = meanbias.copy()
slopes = meanbias.copy()
meansig = meanbias.copy()
obssig = meanbias.copy()

fig, ax  = pl.subplots(1,2,figsize=(11,5)) # (12,9) for 6 panels

for i in range(2):
    #gmst_ts, e5_gmst, eb_gmst, data_ssn, era5_ssn, eobs_ssn = \
    #                                    season_outs(season[i],gmst,era5_gmst,\
    #                                eobs_gmst,modeldata,era5data,eobsdata)
    
    gmst_ts, e5_gmst, data_ssn, era5_ssn = \
                                    season_outs2(season[i],gmst,era5_gmst,\
                                    modeldata,era5data)
    
    msk_in, title_lab = region_outs(region[i],nrth_eur,sth_eur,newlon,lat)
    
    BIAS_mods = pl.zeros([len(allfiles),lat.size,lon.size])
    
    
    SIGNALS = pl.zeros([len(allfiles),lat.size,lon.size])
    
    for F in range(len(allfiles)):
        # smooth GMST with 15 years low-pass filter
        # use Pandas rolling average
        df = pd.DataFrame(gmst_ts[F])
        gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
        gmst_smth = pl.squeeze(pl.asarray(gmst_smth))
        
        SIGNALS[F] = S_N_ratio(data_ssn[F],lat,lon,gmst_smth)
    
    signals_mean = pl.nanmean(SIGNALS,axis=0)
    
#    if i in (0,1):
#        df = pd.DataFrame(pl.squeeze(eb_gmst))
#        ebgmst_smth = df.rolling(15,min_periods=1,center=True).mean()
#        ebgmst_smth = pl.squeeze(pl.asarray(ebgmst_smth))
#        EBSIG = S_N_ratio(eobs_ssn,eobslat,eobslon,ebgmst_smth)
#        
#        OBSSIG = EBSIG.copy()
#        
#        for F in range(len(allfiles)):
#            BIAS_mods[F] = pl.mean(data_ssn[F],axis=0) - eobs_ssn[F]
#    elif i in (2,3):
#        df = pd.DataFrame(pl.squeeze(e5_gmst))
#        e5gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
#        e5gmst_smth = pl.squeeze(pl.asarray(e5gmst_smth))
#        E5SIG = S_N_ratio(era5_ssn,era5lat,era5lon,e5gmst_smth)
#        
#        OBSSIG = E5SIG.copy()
#        
#        for F in range(len(allfiles)):
#            BIAS_mods[F] = pl.mean(data_ssn[F],axis=0) - era5_ssn[F]
    
    if var == 'mrsos':
        df = pd.DataFrame(pl.squeeze(e5_gmst))
        e5gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
        e5gmst_smth = pl.squeeze(pl.asarray(e5gmst_smth))
        E5SIG = S_N_ratio(era5_ssn,era5lat,era5lon,e5gmst_smth)
        
        OBSSIG = E5SIG.copy()
        
        for F in range(len(allfiles)):
            BIAS_mods[F] = pl.mean(data_ssn[F],axis=0) - era5_ssn[F]
    
    # rearrange the arrays to make slicing and averaging easier:
    sm0 = pl.zeros_like(signals_mean) # temporary array
    sm0[:,:120] = signals_mean[:,120:]
    sm0[:,120:] = signals_mean[:,:120]
    signals_mean = sm0.copy()
    del sm0
    
#    e50 = pl.zeros_like(E5SIG)
#    e50[:,:120] = E5SIG[:,120:]
#    e50[:,120:] = E5SIG[:,:120]
#    E5SIG = e50.copy()
#    del e50
    
    os0 = pl.zeros_like(OBSSIG)
    os0[:,:120] = OBSSIG[:,120:]
    os0[:,120:] = OBSSIG[:,:120]
    OBSSIG = os0.copy()
    del os0
    
    bm0 = pl.zeros_like(BIAS_mods)
    bm0[:,:,:120] = BIAS_mods[:,:,120:]
    bm0[:,:,120:] = BIAS_mods[:,:,:120]
    BIAS_mods = bm0.copy()
    del bm0
    
    sg0 = pl.zeros_like(SIGNALS)
    sg0[:,:,:120] = SIGNALS[:,:,120:]
    sg0[:,:,120:] = SIGNALS[:,:,:120]
    SIGNALS = sg0.copy()
    del sg0
    
    #E5_ave = RegionCalc(msk_in,lon,lat,E5SIG.T,cmipmask.T,areas)
    #EB_ave = RegionCalc(msk_in,lon,lat,EBSIG.T,cmipmask.T,areas)
    OBS_ave = RegionCalc(msk_in,lon,lat,OBSSIG.T,cmipmask.T,areas)
    obssig[i] = OBS_ave
    
    bm_ave = pl.zeros([len(allfiles)])
    sm_ave = pl.zeros([len(allfiles)])
    for F in range(len(allfiles)):
        bm_ave[F] = RegionCalc(msk_in,lon,lat,BIAS_mods[F].T,cmipmask.T,areas)
        sm_ave[F] = RegionCalc(msk_in,lon,lat,SIGNALS[F].T,cmipmask.T,areas)
    
    BIAS_ave = bm_ave.mean()
    meanbias[i] = BIAS_ave
    sig_ave = sm_ave.mean()
    meansig[i] = sig_ave
    
    regout = stats.linregress(bm_ave,sm_ave-OBS_ave)
    corrs[i] = regout.rvalue
    slopes[i] = regout.slope
#if region == 'north':
#    msk_in = RegionMask(nrth_eur,newlon,lat)
#    title_lab = 'North'
#    filelab = 'nrth'
#elif region == 'south':
#    msk_in = RegionMask(sth_eur,newlon,lat)
#    title_lab = 'South'
#    filelab = 'sth'

    if var == 'tas':
        units = 'K'
    elif var == 'pr':
        units = 'mm day$^{-1}$'
    elif var == 'mrsos':
        units = 'm$^3$ m$^{-3}$'

    axx = pl.subplot(1,2,i+1)

    star = axx.plot(BIAS_ave,sig_ave-OBS_ave,marker='*',ms=13,zorder=10,lw=0,
            label='CMIP6 ensemble mean')

    pos = pl.where(sm_ave>0)
    neg = pl.where(sm_ave<0)


    if  i == 0:
        handles = []
        for M in range(len(modnames)):
            dots = axx.plot(bm_ave[M],sm_ave[M],lw=0,ms=0,marker='.',alpha=0,
                            color=None)#,label=nums[M]+'. '+modnames[M])
            handles.append(dots[0])
    axx.plot(bm_ave[pos[0]],sm_ave[pos[0]]-OBS_ave,marker='o',
            zorder=10,lw=0,color='goldenrod',label='S(mod)>0')
    axx.plot(bm_ave[neg[0]],sm_ave[neg[0]]-OBS_ave,marker='o',
            zorder=10,lw=0,color='grey',label='S(mod)<0')
    regline = axx.plot(pl.linspace(-6,6,11),regout[0]*pl.linspace(-6,6,11)+regout[1],
                                    ls='-',lw=1.5,color='lightcoral',zorder=6,
                                    label='line of best fit')


    xaxis = axx.axhline(y=0,ls='--',color='lightgrey')
    yaxis = axx.axvline(x=0,ls='--',color='lightgrey')

    if i == 0:
        pl.xlim(-0.1,0.2); pl.ylim(-0.025,0.005)
    elif i == 1:
        pl.xlim(-0.1,0.02); pl.ylim(-0.02,0.02)
#   elif i == 2:
#        pl.xlim(-1.0,0.4); pl.ylim(-0.6,0.2)
#    elif i == 3:
#        pl.xlim(-0.5,1.2); pl.ylim(-0.5,0.5)
#    elif i == 4:
#        pl.xlim(-2.0,5.1); pl.ylim(-2.0,2.1)
#    elif i == 5:
#        pl.xlim(-2.0,5.5); pl.ylim(-2.0,2.0)


    crd_x = bm_ave
    crd_y = sm_ave-OBS_ave
    texts = [axx.text(crd_x[j],crd_y[j],nums[j],zorder=50,size=9) for j in range(len(modnames))]
    adjust_text(texts,avoid_text=True,avoid_self=False,avoid_points=False,
                force_objects=(0.01,0.03),add_objects=[star[0]])
            #arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

    axx.tick_params(axis='both', which='major', labelsize=9,direction='in')
    #if i in (4,5):
    #    pl.xlabel('temperature bias ('+units+')',fontsize=10,labelpad=-0.5)
    
    if i in (0,2,4):
        pl.ylabel('S(mod) - S(obs) ('+units+')',fontsize=10,labelpad=-0.5)

    if i == 0:
        axx.annotate('cor = '+str(round(regout[2],2)),(0.07,-0.022),fontsize=11)
        axx.annotate('S(EOBS) = '+str(round(OBS_ave,3))+' '+units,(0.01,0.003),
                                                                fontsize=11)
    elif i == 1:
        axx.annotate('cor = '+str(round(regout[2],2)),(-0.03,-0.012),fontsize=11)
        axx.annotate('S(EOBS) = '+str(round(OBS_ave,3))+' '+units,(-0.095,0.017),
                                                                fontsize=11)
#    elif i == 2:
#        axx.annotate('cor = '+str(round(regout[2],2)),(-0.35,0.05),fontsize=11)
#        axx.annotate('S(ERA5) = '+str(round(OBS_ave,2))+' '+units,(-0.55,0.16),
#                                                                fontsize=11)
#    elif i == 3:
#        axx.annotate('cor = '+str(round(regout[2],2)),(0.75,0.2),fontsize=11)
#        axx.annotate('S(ERA5) = '+str(round(OBS_ave,2))+' '+units,(0.45,0.43),
#                                                                fontsize=11)
#    elif i == 4:
#        axx.annotate('cor = '+str(round(regout[2],2)),(3.5,1.1),fontsize=11)
#        axx.annotate('S(ERA5) = '+str(round(OBS_ave,2))+' '+units,(0.05,1.8),
#                                                                fontsize=11)
#    elif i == 5:
#        axx.annotate('cor = '+str(round(regout[2],2)),(3.5,1.1),fontsize=11)
#        axx.annotate('S(EOBS) = '+str(round(OBS_ave,2))+' '+units,(0.05,1.7),
#                                                                fontsize=11)
    pl.title(title_lab+' Europe '+season[i]+' '+var+' bias, '+var+' signal',
             fontsize=10,y=0.985)
    
    if i == 1:
        axx.legend(loc=3,fontsize=8)

    #pl.subplots_adjust(top=0.95,bottom=0.09)
leglabs = [nums[X]+'. '+modnames[X] for X in range(len(modnames))]
fig.legend(handles,labels=leglabs,loc=(0.06,0.01),ncol=9,columnspacing=0.7,
           handlelength=0,borderpad=0.5,fontsize=9)

pl.tight_layout()
#pl.subplots_adjust(top=0.975,bottom=0.11,hspace=0.16)
pl.subplots_adjust(top=0.96,bottom=0.14,hspace=0.16)

fig.text(0.005,0.96,'(a)',size=10)
fig.text(0.50,0.96,'(b)',size=10)
#fig.text(0.005,0.50,'(c)',size=10)
#fig.text(0.50,0.50,'(d)',size=10)
#fig.text(0.005,0.35,'(e)',size=10)
#fig.text(0.50,0.35,'(f)',size=10)

#pl.savefig(indecis+'figures/'+var+'_scatters_panels.png',dpi=400)
#pl.savefig(indecis+'figures/'+var+'_scatters_panels.pdf',dpi=400)
#
#print 'mean bias = ', BIAS_ave, ' mm/day'
#print 'correlation = ', regout.rvalue
#print 'slope = ', regout.slope, ' (mm/day)/(mm/day)'
#print 'mean signal = ', sig_ave, ' mm/day' 
#print 'EOBS signal = ', E5_ave, ' mm/day'