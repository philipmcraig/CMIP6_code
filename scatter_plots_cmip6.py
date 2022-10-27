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

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
cmipdir = ncasdir + 'CMIP6/'
era5dir = ncasdir + 'ERA5/'
eobsdir = ncasdir + 'EOBS/'
indecis = ncasdir + 'INDECIS/'
maskdir = cmipdir + '/masks/'

var = 'tas'
var2 = 'tg'
var3 = 't2m'
season = ['DJF','JJA']
region = ['south','north']

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
modelsomo = modeldata.copy()
modnames = []
gmst = pl.zeros([len(tasfiles),len(time)])

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var)[:,60:,:]).data
    ncfile.close
    
#    precfile = xr.open_dataset(precfiles[nci])
#    modelprec[nci,:,:,:] = xr.DataArray(precfile.pr)
#    precfile.close()
    
#    tempfile = xr.open_dataset(tempfiles[nci])
#    modeltemp[nci,:,:,:] = xr.DataArray(tempfile.tas)
#    tempfile.close()
    
    somofile = xr.open_dataset(somofiles[nci])
    modelsomo[nci,:,:,:] = xr.DataArray(somofile.mrsos[:,60:,:]).data
    somofile.close()
    
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
#era5temp = xr.DataArray(era5file.t2m)
#era5prec = xr.DataArray(era5file.tp)
era5somo = xr.DataArray(era5file.swvl1[:,60:,:]).data
era5lon = xr.DataArray(era5file.lon).data
era5lat = xr.DataArray(era5file.lat[60:]).data
era5time = xr.DataArray(era5file.time).data
era5file.close()

#eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobsdata = xr.DataArray(getattr(eobsfile,var2))
#eobstime = xr.DataArray(eobsfile.time)
#eobslat = xr.DataArray(eobsfile.lat)
#eobslon = xr.DataArray(eobsfile.lon)
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

#nc_gmst2 = xr.open_dataset(ncasdir+'hadcrut5_seasmean.nc')
#eobs_gmst = xr.DataArray(nc_gmst2.tas_mean)
#gmst_time = xr.DataArray(nc_gmst2.time)
#nc_gmst2.close()

if season == 'JJA':
#    data_em_tm = pl.mean(cmipmean[1:-1:4,:,:],axis=0)#/(0.1*1000)#*86400#
#    data_es_tm = pl.mean(cmipsprd[1::4,:,:],axis=0)#/(0.1*1000)#*86400
#    eobs_mn = pl.mean(eobsdata[1:-26:4],axis=0)# + 273.15
#    era5_mn = pl.mean(era5data[2:-1:4],axis=0)*1000
    gmst_ts = gmst[:,1:-1:4]
    e5_gmst = era5_gmst[2:-1:4]
#    eb_gmst = eobs_gmst[401:-26:4] + 273.15
    data_ssn = modeldata[:,1:-1:4]#/(0.1*1000)#*86400#
#    c6tas_ssn = modeltemp[:,1:-1:4]
#    c6pr_ssn = modelprec[:,1:-1:4]*86400
    c6sm_ssn = modelsomo[:,1:-1:4]/(0.1*1000)
    era5_ssn = era5data[2:-1:4]#*1000
#    eobs_ssn = eobsdata[1:-26:4,:,:] + 273.15
#    e5tmp_ssn = era5temp[2:-14:4]
#    ebtg_ssn = eobstemp[1:-1:4,:,:] + 273.15
#    e5tp_ssn = era5prec[2:-1:4]*1000
#    ebrr_ssn = eobsprec[1:-26:4]
    e5sm_ssn = era5somo[2:-1:4]
elif season == 'DJF':
#    data_em_tm = pl.mean(cmipmean[3:-1:4,:,:],axis=0)/(0.1*1000)#*86400#
#    data_es_tm = pl.mean(cmipsprd[3::4,:,:],axis=0)#/(0.1*1000)#*86400
#    eobs_mn = pl.mean(eobsdata[3:-26:4],axis=0)# + 273.15
#    era5_mn = pl.mean(era5data[4:-14:4,:,:],axis=0)*1000
    gmst_ts = gmst[:,3:-1:4]
    e5_gmst = era5_gmst[4:-1:4]
#    eb_gmst = eobs_gmst[403:-26:4] + 273.15
    data_ssn = modeldata[:,3:-1:4]#/(0.1*1000)#*86400#
 #   c6tas_ssn = modeltemp[:,3:-1:4]
#    c6pr_ssn = modelprec[:,3:-1:4]*86400
    c6sm_ssn = modelsomo[:,3:-1:4]/(0.1*1000)
    era5_ssn = era5data[4:-1:4]#*1000
#    eobs_ssn = eobsdata[3:-26:4,:,:] + 273.15
#    e5tmp_ssn = era5temp[4:-14:4]
#    ebtg_ssn = eobstemp[3:-1:4,:,:] + 273.15
#    e5tp_ssn = era5prec[4:-1:4]*1000
#    ebrr_ssn = eobsprec[3:-26:4]
    e5sm_ssn = era5somo[4:-1:4]

#BIAS_mn = data_em_tm - era5_mn

BIAS_mods = pl.zeros([len(allfiles),lat.size,lon.size])
for i in range(len(allfiles)):
    BIAS_mods[i] = pl.mean(data_ssn[i],axis=0) - era5_ssn[i]

###############################################################################
SIGNALS = pl.zeros([len(allfiles),lat.size,lon.size])
#pvals = SN_ratios.copy()

for i in range(len(allfiles)):
    # smooth GMST with 15 years low-pass filter
    # use Pandas rolling average
    df = pd.DataFrame(gmst_ts[i])
    gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
    gmst_smth = pl.squeeze(pl.asarray(gmst_smth))
    
    SIGNALS[i] = S_N_ratio(c6sm_ssn[i],lat,lon,gmst_smth)

signals_mean = pl.nanmean(SIGNALS,axis=0)
###############################################################################

###############################################################################
df = pd.DataFrame(pl.squeeze(e5_gmst))
e5gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
e5gmst_smth = pl.squeeze(pl.asarray(e5gmst_smth))

#df = pd.DataFrame(pl.squeeze(eb_gmst))
#ebgmst_smth = df.rolling(15,min_periods=1,center=True).mean()
#ebgmst_smth = pl.squeeze(pl.asarray(ebgmst_smth))

E5SIG = S_N_ratio(e5sm_ssn,era5lat,era5lon,e5gmst_smth)
#EBSIG = S_N_ratio(ebrr_ssn,eobslat.data,eobslon.data,ebgmst_smth)

###############################################################################

# rearrange the arrays to make slicing and averaging easier:
sm0 = pl.zeros_like(signals_mean) # temporary array
sm0[:,:120] = signals_mean[:,120:]
sm0[:,120:] = signals_mean[:,:120]
signals_mean = sm0.copy()
del sm0

e50 = pl.zeros_like(E5SIG)
e50[:,:120] = E5SIG[:,120:]
e50[:,120:] = E5SIG[:,:120]
E5SIG = e50.copy()
del e50

#eb0 = pl.zeros_like(EBSIG)
#eb0[:,:120] = EBSIG[:,120:]
#eb0[:,120:] = EBSIG[:,:120]
#EBSIG = eb0.copy()
#del eb0

#b0 = pl.zeros_like(BIAS_mn)
#b0[:,:120] = BIAS_mn[:,120:]
#b0[:,120:] = BIAS_mn[:,:120]
#BIAS_mn = b0.copy()
#del b0

#newlon = pl.zeros_like(lon.data)
#newlon[:120] = lon.data[120:] - 360.
#newlon[120:] = lon.data[:120]
newlon = lon - 180

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
###############################################################################

# split at 48N, remember that CMIP6 files are global
#split_ind = pc.NearestIndex(lat.data,48)
#north_ind = pc.NearestIndex(lat.data,71)
#south_ind = pc.NearestIndex(lat.data,36)
#west_ind = pc.NearestIndex(newlon,-11)
#east_ind = pc.NearestIndex(newlon,40)

nrth_eur = [(-5,48),(-11,51.8),(-11,58),(16.7,71),(28.9,71),(28.9,48)]
sth_eur = [(-11,36),(-11,43.8),(-5,48),(28.9,48),(28.6,41),(25.1,41),(25.1,36),
           (15,36),(10.3,38.1),(2.7,38.1),(-5.6,36)]

#proj = ccrs.PlateCarree()
#ext = [-15,42,35,70]
#
#axx = pl.axes(projection=proj,extent=ext)
#axx.coastlines(linewidth=0.5,resolution='50m')
#borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
#                                           '50m',edgecolor='grey',
#                                        facecolor='none')
#axx.add_feature(borders_50m,linewidth=0.5,zorder=5)
#
#poly_nrth = Polygon(nrth_eur,fc='none')
#pn = PatchCollection([poly_nrth])
#pn.set_edgecolor('b')
#pn.set_facecolor('none')
#pn.set_linewidth(2)
#axx.add_collection(pn)
#
#poly_sth = Polygon(sth_eur,fc='none')
#ps = PatchCollection([poly_sth])
#ps.set_edgecolor('b')
#ps.set_facecolor('none')
#ps.set_linewidth(2)
#axx.add_collection(ps)
#
#pl.tight_layout()
#pl.savefig(indecis+'/figures/north_south_europe_regions.png')
###############################################################################
# calculate area-averages of model signal, observed signal & temperature bias

areas = AreasCalc2(lon,lat)
#nrth_msk = RegionMask(nrth_eur,lon,lat)
#sth_msk = RegionMask(sth_eur,lon,lat)

if region == 'north':
    msk_in = RegionMask(nrth_eur,newlon,lat)
    title_lab = 'North'
    filelab = 'nrth'
elif region == 'south':
    msk_in = RegionMask(sth_eur,newlon,lat)
    title_lab = 'South'
    filelab = 'sth'

if var == 'tas':
    units = 'K'
elif var == 'pr':
    units = 'mm day$^{-1}$'
elif var == 'mrsos':
    units = 'm$^3$ m$^{-3}$'

E5_ave = RegionCalc(msk_in,lon,lat,E5SIG.T,cmipmask.T,areas)
#EB_ave = RegionCalc(msk_in,lon,lat,EBSIG.T,cmipmask.T,areas)

bm_ave = pl.zeros([len(allfiles)])
sm_ave = pl.zeros([len(allfiles)])
for i in range(len(allfiles)):
    bm_ave[i] = RegionCalc(msk_in,lon,lat,BIAS_mods[i].T,cmipmask.T,areas)
    sm_ave[i] = RegionCalc(msk_in,lon,lat,SIGNALS[i].T,cmipmask.T,areas)

BIAS_ave = bm_ave.mean()#RegionCalc(msk_in,lon,lat,BIAS_mn.T,cmipmask.T,areas)
sig_ave = sm_ave.mean()#RegionCalc(msk_in,lon,lat,signals_mean.T,cmipmask.T,areas)

regout = stats.linregress(bm_ave,sm_ave-E5_ave)

###############################################################################
# plot model signal minus observed signal as function of temperature bias

fig, ax  = pl.subplots(figsize=(9,6))

star = ax.plot(BIAS_ave,sig_ave-E5_ave,marker='*',ms=13,zorder=10,lw=0,
        label='CMIP6 ensemble mean')

pos = pl.where(sm_ave>0)
neg = pl.where(sm_ave<0)

#scat = ax.scatter(bm_sth_ave,sm_sth_ave-E5_sth_ave,marker='o',zorder=10,
#                  c=scatcol)
ax.plot(bm_ave[pos[0]],sm_ave[pos[0]]-E5_ave,marker='o',
        zorder=10,lw=0,color='darkgoldenrod',label='S(mod)>0')
ax.plot(bm_ave[neg[0]],sm_ave[neg[0]]-E5_ave,marker='o',
        zorder=10,lw=0,color='k',label='S(mod)<0')
regline = ax.plot(pl.linspace(-6,6,11),regout[0]*pl.linspace(-6,6,11)+regout[1],
                                ls='-',lw=1.5,color='r',zorder=6,
                                label='line of best fit')

#for i, txt in enumerate(modnames):
#    pl.annotate(txt,(bm_nrth_ave[i]+0.02,sm_nrth_ave[i]-E5_nrth_ave+0.02),
#                size=9,zorder=10)

xaxis = ax.axhline(y=0,ls='--',color='lightgrey')
yaxis = ax.axvline(x=0,ls='--',color='lightgrey')

pl.xlim(-2.5,5.00); pl.ylim(-0.02,0.02)


crd_x = bm_ave
crd_y = sm_ave-E5_ave
texts = [ax.text(crd_x[i],crd_y[i],modnames[i],zorder=50) for i in range(len(modnames))]
adjust_text(texts,avoid_text=False,avoid_self=False,avoid_points=True,force_objects=(0.0,0.001),
            add_objects=[regline[0]])
            #arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

ax.tick_params(axis='both', which='major', labelsize=12)
pl.xlabel('temperature bias ('+units+')',fontsize=12)
pl.ylabel('S(mod) - S(ERA5) (m$^{3}$ m$^{-3}$)',fontsize=12)


#axins = inset_axes(ax, width=1.3, height=0.9,projection=proj,extent=ext)
#ax1 = pl.subplot(ig[0],projection=proj,extent=ext)

#axins = fig.add_axes([0.65,0.15,0.3,0.3],projection=proj,extent=ext,alpha=0.1)#[0.65,0.65,0.3,0.3]
#axins.coastlines(linewidth=0.5,resolution='50m')
#polygon = Polygon(sth_eur,fc='none')
#p = PatchCollection([polygon])
#p.set_edgecolor('b')
#p.set_facecolor('none')
#p.set_linewidth(2)
#axins.add_collection(p)

ax.legend(loc=3)

pl.tight_layout()

#ax.annotate(season,(0.02,0.95),xycoords='axes fraction',
#            fontsize=16,bbox={'facecolor':'w'})
ax.annotate('cor = '+str(round(regout[2],2)),(2.8,-0.013),fontsize=16)
ax.annotate('S(ERA5) = '+str(round(E5_ave,3))+' '+units,(0.05,0.0185),fontsize=14)
pl.title(title_lab+' Europe '+season+' '+var+' bias, mrsos signal',fontsize=12)

pl.subplots_adjust(top=0.95,bottom=0.09)

#pl.savefig(indecis+'figures/'+var+'_mrsos_scatter_era5_'+filelab+'_'+season.lower()+'.png',
#           dpi=360)

print 'mean bias = ', BIAS_ave, ' mm/day'
print 'correlation = ', regout.rvalue
print 'slope = ', regout.slope, ' (mm/day)/(mm/day)'
print 'mean signal = ', sig_ave, ' mm/day' 
print 'EOBS signal = ', E5_ave, ' mm/day'