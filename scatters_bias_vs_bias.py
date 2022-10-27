# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:28:13 2021

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

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
cmipdir = ncasdir + 'CMIP6/'
era5dir = ncasdir + 'ERA5/'
eobsdir = ncasdir + 'EOBS/'
indecis = ncasdir + 'INDECIS/'
maskdir = cmipdir + '/masks/'

var = ['pr','mrsos']
#var2 = ['tg','rr']
#var3 = 'tp'
season = 'DJF'
region = 'north'
obs = 'ERA5'

if obs == 'ERA5':
    var2 = ['tp','swvl1']
elif obs == 'EOBS':
    var2 = ['tg','rr']

meanfile = xr.open_dataset(cmipdir + 'standard_grid/tas_cmip6_ensmean_new.nc')
##lat = xr.DataArray(meanfile.lat)
##lon = xr.DataArray(meanfile.lon)
time = xr.DataArray(meanfile.time)
#meantemp = xr.DataArray(getattr(meanfile,var[0]))
meanfile.close()
#
#meanfile = xr.open_dataset(cmipdir + 'standard_grid/test_prmean_nco.nc')#'+var[1]+'_cmip6_ensmean.nc')
#meanprec = xr.DataArray(getattr(meanfile,var[1]))
#time = xr.DataArray(meanfile.time)
#meanfile.close()

maskfile = xr.open_dataset(maskdir+'lsmask_cmip6_s1.5.nc')
cmipmask = xr.DataArray(maskfile.topo)
lat = xr.DataArray(maskfile.lat)
lon = xr.DataArray(maskfile.lon)
maskfile.close()

cm0 = pl.zeros_like(cmipmask.data) # temporary array
cm0[:,:120] = cmipmask.data[:,120:]
cm0[:,120:] = cmipmask.data[:,:120]
cmipmask = cm0.copy()
del cm0

tempfiles = glob.glob(cmipdir+'standard_grid/'+var[0]+'_1950-2014_new/*')
#tasfiles = glob.glob(cmipdir+'global_mean_tas/*')
precfiles = glob.glob(cmipdir+'standard_grid/'+var[1]+'_1950-2014_new/*')

tempdata = pl.zeros([len(tempfiles),len(time),lat.size,lon.size])
precdata = tempdata.copy()
modnames = []
#gmst = pl.zeros([len(tasfiles),len(time)])

for nci in range(len(tempfiles)):
    tempfile = xr.open_dataset(tempfiles[nci])
    tempdata[nci,:,:,:] = xr.DataArray(getattr(tempfile,var[0]))
    tempfile.close()
    
    precfile = xr.open_dataset(precfiles[nci])
    precdata[nci,:,:,:] = xr.DataArray(getattr(precfile,var[1]))
    precfile.close()
    
    split1 = tempfiles[nci].split('/')
    split2 = split1[-1].split('_')[2]
    modnames.append(split2)
    
#    tasfile = xr.open_dataset(tasfiles[nci])
#    gmst[nci,:] = pl.squeeze(xr.DataArray(tasfile.tas))
#    tasfile.close()

modnames[1] = modnames[1][:-2]
modnames[2] = modnames[2][:7]
modnames[6] = modnames[6][:-2]
modnames[9] = modnames[9][:-4]
modnames[12] = modnames[12][:-5]
modnames[13] = modnames[13][:-2]
modnames[14] = modnames[14][:-3]
modnames[16] = modnames[16][:6]

#eobsfile = xr.open_dataset(eobsdir+var2[0]+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobstemp = xr.DataArray(getattr(eobsfile,var2[0]))
#eobstime = xr.DataArray(eobsfile.time)
#eobslat = xr.DataArray(eobsfile.lat)
#eobslon = xr.DataArray(eobsfile.lon)
#eobsfile.close()

#eobsfile = xr.open_dataset(eobsdir+var2[1]+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobsprec = xr.DataArray(getattr(eobsfile,var2[1]))
#eobsfile.close()

era5file = xr.open_dataset(era5dir+'era5_land_seasmean_s1.5.nc')
era5temp = xr.DataArray(getattr(era5file,var2[0]))
era5prec = xr.DataArray(getattr(era5file,var2[1]))
#era5somo = xr.DataArray(era5file.swvl1)
era5lon = xr.DataArray(era5file.lon)
era5lat = xr.DataArray(era5file.lat)
era5time = xr.DataArray(era5file.time)
era5file.close()

if season == 'JJA':
    #temp_em_mn = pl.nanmean(meantemp[1:-1:4,:,:],axis=0)
    #prec_em_mn = pl.nanmean(meanprec[1:-1:4],axis=0)*86400
    
#    eobs_tmp_mn = pl.nanmean(eobstemp[1:-1:4,:,:],axis=0) + 273.15
#    eobs_prc_mn = pl.nanmean(eobsprec[1:-26:4,:,:],axis=0)
    
    era5_tmp_mn = pl.nanmean(era5temp[2:-14:4],axis=0)*1000
    era5_prc_mn = pl.nanmean(era5prec[2:-14:4])#/(0.1*1000)
    
    c6tas_ssn = tempdata[:,1:-1:4]*86400
    c6pr_ssn = precdata[:,1:-1:4]/(0.1*1000)#*86400
    
    e5t2m_ssn = era5temp[2:-14:4,:,:]*1000
    e5tp_ssn = era5prec[2:-14:4,:,:]#*1000
    
#    ebtg_ssn = eobstemp[1:-1:4,:,:] + 273.15
#    ebrr_ssn = eobsprec[1:-26:4]
elif season == 'DJF':
    #temp_em_mn = pl.nanmean(meantemp[3:-1:4,:,:],axis=0)
    #prec_em_mn = pl.nanmean(meanprec[3:-1:4],axis=0)*86400
    
#    eobs_tmp_mn = pl.nanmean(eobstemp[3:-1:4,:,:],axis=0) + 273.15
#    eobs_prc_mn = pl.nanmean(eobsprec[3:-26:4,:,:],axis=0)
    
    era5_tmp_mn = pl.nanmean(era5temp[4:-14:4],axis=0)*1000
    era5_prc_mn = pl.nanmean(era5prec[4:-14:4])#/(0.1*1000)#*1000
    
    c6tas_ssn = tempdata[:,3:-1:4]*86400
    c6pr_ssn = precdata[:,3:-1:4]/(0.1*1000)#*86400
    
    e5t2m_ssn = era5temp[4:-14:4,:,:]*1000
    e5tp_ssn = era5prec[4:-14:4,:,:]#*1000
    
#    ebtg_ssn = eobstemp[3:-1:4,:,:] + 273.15
#    ebrr_ssn = eobsprec[3:-26:4]


#TEMPBIAS = temp_em_mn - era5_tmp_mn
#PRECBIAS = prec_em_mn - era5_prc_mn

BIAS_tmp_mods = pl.zeros([len(tempfiles),lat.size,lon.size])
BIAS_prc_mods = BIAS_tmp_mods.copy()
for i in range(len(tempfiles)):
    BIAS_tmp_mods[i] = pl.mean(c6tas_ssn[i],axis=0) - e5t2m_ssn[i]
    BIAS_prc_mods[i] = pl.mean(c6pr_ssn[i],axis=0) - e5tp_ssn[i]


###############################################################################

newlon = lon.data - 180

#b0 = pl.zeros_like(TEMPBIAS)
#b0[:,:120] = TEMPBIAS[:,120:]
#b0[:,120:] = TEMPBIAS[:,:120]
#TEMPBIAS = b0.copy()
#del b0
#
#b0 = pl.zeros_like(PRECBIAS)
#b0[:,:120] = PRECBIAS[:,120:]
#b0[:,120:] = PRECBIAS[:,:120]
#PRECBIAS = b0.copy()
#del b0

bm0 = pl.zeros_like(BIAS_tmp_mods)
bm0[:,:,:120] = BIAS_tmp_mods[:,:,120:]
bm0[:,:,120:] = BIAS_tmp_mods[:,:,:120]
BIAS_tmp_mods = bm0.copy()
del bm0

bm0 = pl.zeros_like(BIAS_prc_mods)
bm0[:,:,:120] = BIAS_prc_mods[:,:,120:]
bm0[:,:,120:] = BIAS_prc_mods[:,:,:120]
BIAS_prc_mods = bm0.copy()
del bm0

#tm0 = pl.zeros_like(temp_em_mn)
#tm0[:,:120] = temp_em_mn[:,120:]
#tm0[:,120:] = temp_em_mn[:,:120]
#temp_em_mn = tm0.copy()
#del tm0
#
#ts0 = pl.zeros_like(pl.mean(c6tas_ssn,axis=1))
#ts0[:,:,:120] = pl.mean(c6tas_ssn,axis=1)[:,:,120:]
#ts0[:,:,120:] = pl.mean(c6tas_ssn,axis=1)[:,:,:120]
#c6tas_ssn_mn = ts0.copy()
#del ts0

et0 = pl.zeros_like(era5_tmp_mn)
et0[:,:120] = era5_tmp_mn[:,120:]
et0[:,120:] = era5_tmp_mn[:,:120]
era5_tmp_mn = et0.copy()
del et0

###############################################################################

nrth_eur = [(-5,48),(-11,51.8),(-11,58),(16.7,71),(28.9,71),(28.9,48)]
sth_eur = [(-11,36),(-11,43.8),(-5,48),(28.9,48),(28.6,41),(25.1,41),(25.1,36),
           (15,36),(10.3,38.1),(2.7,38.1),(-5.6,36)]

areas = AreasCalc2(lon.data,lat.data)

if region == 'north':
    msk_in = RegionMask(nrth_eur,newlon,lat)
    title_lab = 'North'
    filelab = 'nrth'
elif region == 'south':
    msk_in = RegionMask(sth_eur,newlon,lat)
    title_lab = 'South'
    filelab = 'sth'

if var[0] == 'tas' and var[1] == 'pr':
    units1 = 'K'
    units2 = 'mm day$^{-1}$'
elif var[0] == 'pr' and var[1] == 'tas':
    units1 = 'mm day$^{-1}$'
    units2 = 'K'
elif var[0] == 'tas' and var[1] == 'mrsos':
    units1 = 'K'
    units2 = 'm$^3$ m$^{-3}$'
elif var[0] == 'pr' and var[1] == 'mrsos':
    units1 = 'mm day$^{-1}$'
    units2 = 'm$^3$ m$^{-3}$'

bm_tmp_ave = pl.zeros([len(tempfiles)])
bm_prc_ave = pl.zeros([len(tempfiles)])
for i in range(len(tempfiles)):
    bm_tmp_ave[i] = RegionCalc(msk_in,lon,lat,BIAS_tmp_mods[i].T,cmipmask.T,areas)
    bm_prc_ave[i] = RegionCalc(msk_in,lon,lat,BIAS_prc_mods[i].T,cmipmask.T,areas)

TEMPBIAS_ave = bm_tmp_ave.mean()#RegionCalc(msk_in,lon,lat,TEMPBIAS.T,cmipmask.T,areas)
PRECBIAS_ave = bm_prc_ave.mean()#RegionCalc(msk_in,lon,lat,PRECBIAS.T,cmipmask.T,areas)

#tmp_ave = RegionCalc(msk_in,lon,lat,temp_em_mn.T,cmipmask.T,areas)
e5_tmp_ave = RegionCalc(msk_in,lon,lat,era5_tmp_mn.T,cmipmask.T,areas)

###############################################################################

regout = stats.linregress(bm_tmp_ave,bm_prc_ave)

fig, ax  = pl.subplots(figsize=(9,6))

star = ax.plot(TEMPBIAS_ave,PRECBIAS_ave,marker='*',ms=13,zorder=10,lw=0,
        label='CMIP6 ensemble mean')

scat = ax.scatter(bm_tmp_ave,bm_prc_ave,marker='o',zorder=10)

regline = ax.plot(pl.linspace(-6,6,11),regout[0]*pl.linspace(-6,6,11)+regout[1],
                                ls='-',lw=1.5,color='r',zorder=6,
                                label='line of best fit')

xaxis = ax.axhline(y=0,ls='--',color='lightgrey')
yaxis = ax.axvline(x=0,ls='--',color='lightgrey')

pl.xlim(-3.0,5.5); pl.ylim(-0.5,0.5)


crd_x = bm_tmp_ave
crd_y = bm_prc_ave
texts = [ax.text(crd_x[i],crd_y[i],modnames[i],zorder=50) for i in range(len(modnames))]
adjust_text(texts,avoid_text=False,avoid_self=False,avoid_points=True)#,force_objects=(0.0001,0.0006),
            #add_objects=[regline[0]])
            #arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

ax.tick_params(axis='both', which='major', labelsize=12)
pl.xlabel('temperature bias ('+units1+')',fontsize=12)
pl.ylabel('precipitation bias ('+units2+')',fontsize=12)

ax.legend(loc=3)

pl.tight_layout()

ax.annotate('cor = '+str(round(regout[2],2)),(2.8,1.35),fontsize=16)
#ax.annotate('S(ERA5) = '+str(round(EB_ave,2))+' '+units,(0.05,1.8),fontsize=14)
pl.title(title_lab+' Europe '+season+' '+var[0]+' bias vs '+var[1]+' bias (ERA5)',fontsize=12)

pl.subplots_adjust(top=0.95,bottom=0.09)

#pl.savefig(indecis+'figures/'+var[0]+'_'+var[1]+'_'+'biasvbias_era5_'+filelab+'_'+season.lower()+'.png',
#           dpi=360)

print 'mean temp bias = ', TEMPBIAS_ave, ' K'
print 'mean prec bias = ', PRECBIAS_ave, ' mm/day'
print 'correlation = ', regout.rvalue
print 'slope = ', regout.slope, ' (mm/day)/K'
#print 'mean signal = ', sig_ave, ' mm/day' 
#print 'EOBS signal = ', EB_ave, ' mm/day'