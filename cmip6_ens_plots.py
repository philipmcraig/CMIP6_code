# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 15:55:05 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import pandas as pd
import glob
from scipy import stats
import cartopy.crs as ccrs
import cartopy.util as util
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
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
    noise = pl.std(R,axis=0)
    #noise_djf = pl.std(R[1],axis=0)
    
    # signal-to-noise ratio
    # S/N
    sig_noi_ratio = signal/noise
    #sig_noi_ratio_djf = signal_djf/noise_djf
    
    #SN_ratios = pl.array([sig_noi_ratio_jja,sig_noi_ratio_djf])
    
    return sig_noi_ratio, pvalue

def CheckBiasSign(cmipdir,BIAS2,era5_seasmean,season,var,var3):
    """
    """
    allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014/*')
    
    #BIAS_mean = cmip_tm - eobs_mn
    BIAS_mods2 = pl.zeros([len(allfiles),BIAS2.shape[0],BIAS2.shape[1]])
    #BIAS_mods2 = BIAS_mods1.copy()
    #runtot = pl.zeros_like(BIAS_mean) # running total
    CHECK2 = pl.zeros([BIAS2.shape[0],BIAS2.shape[1]])
    CHECK2[:,:] = pl.float32('nan')
    #CHECK2 = CHECK1.copy()
    
    for m in range(len(allfiles)):
        filepath = allfiles[m]
        ncfile = xr.open_dataset(filepath)
        cmipdata = xr.DataArray(getattr(ncfile,var))
        ncfile.close()
        
        if season == 'JJA':
            model_tm = pl.nanmean(cmipdata[1::4,:,:],axis=0)
        elif season == 'DJF':
            model_tm = pl.nanmean(cmipdata[3::4,:,:],axis=0)
        
        if var == 'pr':
            model_tm = model_tm*86400
        elif var == 'mrsos':
            model_tm = model_tm/(0.1*1000)
        
        #BIAS_mods1[m] = model_tm - eobs_seasmean
        BIAS_mods2[m] = model_tm - era5_seasmean
    
    for i in range(BIAS_mods2.shape[1]):
        for j in range(BIAS_mods2.shape[2]):
#            point1 = BIAS_mods1[:,i,j]
#            if pl.all(pl.isnan(point1)) == True:
#                continue
#            else:
#                meansign1 = pl.sign(BIAS1[i,j])
#                meansign1 = int(meansign1.data)
#                SIGNS1 = pl.sign(point1)
#                samesign1 = pl.where(SIGNS1==meansign1)
#                
#                #pos_count1 = len(list(filter(lambda x: (x > 0), point1)))
#                #neg_count1 = len(list(filter(lambda x: (x < 0), point1)))
#                #zer_count1 = len(list(filter(lambda x: (x == 0), point1)))
#                
#                #A1 = pl.array([pos_count1,neg_count1,zer_count1])
#                
#                if samesign1[0].size > 17.*(2./3.):
#                    CHECK1[i,j] = 1.
            
            ###################################################################
            point2 = BIAS_mods2[:,i,j]
            if pl.all(pl.isnan(point2)) == True:
                continue
            else:
                meansign2 = pl.sign(BIAS2[i,j])
                meansign2 = int(meansign2.data)
                SIGNS2 = pl.sign(point2)
                samesign2 = pl.where(SIGNS2==meansign2)
#                pos_count2 = len(list(filter(lambda x: (x > 0), point2)))
#                neg_count2 = len(list(filter(lambda x: (x < 0), point2)))
#                zer_count2 = len(list(filter(lambda x: (x == 0), point2)))
#                
#                A2 = pl.array([pos_count2,neg_count2,zer_count2])
                
                if samesign2[0].size > 17.*(2./3.):
                    CHECK2[i,j] = 1.
    
    return CHECK2
    
def GridLines(ax,top,bottom,left,right):
    """
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
                  zorder=6)
    gl.xlabels_top = top; gl.xlabels_bottom = bottom
    gl.ylabels_left = left; gl.ylabels_right = right
    gl.xlocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20,30,40,50])
    gl.ylocator = mticker.FixedLocator([20,30,40,50,60,70,80])
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'k','size':10}
    gl.ylabel_style = {'color': 'k','size':10}
    
    return None

class MidpointNormalize(Normalize):
    """Found this on the internet. Use to centre colour bar at zero.
    """
    def __init__(self, vmin=-0.3, vmax=0.6, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        a, b = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return pl.ma.masked_array(pl.interp(value, a, b))

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
eobsdir = ncasdir + 'EOBS/'

var = 'mrsos'
#var2 = 'tg'
var3 = 'swvl1'
season = 'DJF'

meanfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_ensmean.nc')
lat = xr.DataArray(meanfile.lat)
lon = xr.DataArray(meanfile.lon)
time = xr.DataArray(meanfile.time)
cmipmean = xr.DataArray(getattr(meanfile,var))
meanfile.close()

allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014_new/*')
tasfiles = glob.glob(cmipdir+'global_means/*')

modeldata = pl.zeros([len(allfiles),len(time),lat.size,lon.size])

#if var == 'tas':
#    for nci in range(len(allfiles)):
#        ncfile = xr.open_dataset(allfiles[nci])
#        modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var))
#        ncfile.close
#    
#    gmst = modeldata.copy()
#else:
gmst = pl.zeros([len(tasfiles),len(time)])
#modnames = []

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var))
    ncfile.close
    
    tasfile = xr.open_dataset(tasfiles[nci])
    gmst[nci,:] = pl.squeeze(xr.DataArray(tasfile.tas))
    tasfile.close()
    
    split1 = allfiles[nci].split('/')
    split2 = split1[-1].split('_')[2]
#    modnames.append(split2)

sprdfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_enssprd.nc')
cmipsprd = xr.DataArray(getattr(sprdfile,var))
sprdfile.close()

#eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
#eobsdata = xr.DataArray(getattr(eobsfile,var2))
#eobstime = xr.DataArray(eobsfile.time)
#eobsfile.close()

#modnames[1] = modnames[1][:-2]
#modnames[2] = modnames[2][:7]
#modnames[6] = modnames[6][:-2]
#modnames[9] = modnames[9][:-4]
#modnames[12] = modnames[12][:-5]
#modnames[13] = modnames[13][:-2]
#modnames[14] = modnames[14][:-3]
#modnames[16] = modnames[16][:6]

era5file = xr.open_dataset(era5dir+'era5_land_seasmean_s1.5.nc')
era5data = xr.DataArray(getattr(era5file,var3))
era5time = xr.DataArray(era5file.time)
era5file.close()

#nc_gmst = xr.open_dataset(cmipdir+'standard_grid/tas_cmip6_ensmean_gm.nc')
#gmst = xr.DataArray(nc_gmst.tas)
#gmst = pl.squeeze(gmst)
#nc_gmst.close()

#depth = pl.zeros([len(modnames)])
#depth[:2] = 0.1; depth[2:] = 0.05

if season == 'JJA':
    data_em_tm = pl.mean(cmipmean[1::4,:,:],axis=0)/(0.1*1000)#*86400
    data_es_tm = pl.mean(cmipsprd[1::4,:,:],axis=0)/(0.1*1000)#*86400
    #eobs_mn = pl.mean(eobsdata[1::4],axis=0) + 273.15
    era5_mn = pl.mean(era5data[2:-14:4],axis=0)#*1000
    gmst_ts = gmst[:,1::4]
    data_ssn = modeldata[:,1::4]#*86400
elif season == 'DJF':
    data_em_tm = pl.mean(cmipmean[3::4,:,:],axis=0)/(0.1*1000)#*86400
    data_es_tm = pl.mean(cmipsprd[3::4,:,:],axis=0)/(0.1*1000)#*86400
    #eobs_mn = pl.mean(eobsdata[3::4],axis=0) + 273.15
    era5_mn = pl.mean(era5data[4:-14:4],axis=0)#*1000
    gmst_ts = gmst[:,3::4]
    data_ssn = modeldata[:,3::4]#*86400

#BIAS1 = data_em_tm - eobs_mn
BIAS2 = data_em_tm - era5_mn

CHECK2 = CheckBiasSign(cmipdir,BIAS2,era5_mn,season,var,var3)


SN_ratios = pl.zeros([len(allfiles),lat.size,lon.size])
pvals = SN_ratios.copy()

for i in range(len(allfiles)):
    # smooth GMST with 15 years low-pass filter
    # use Pandas rolling average
    df = pd.DataFrame(gmst_ts[i])
    gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
    gmst_smth = pl.squeeze(pl.asarray(gmst_smth))
    
    SN_ratios[i], pvals[i] = S_N_ratio(data_ssn[i],lat,lon,gmst_smth)

SN_rat_mean = pl.nanmean(SN_ratios,axis=0)

SIG = pl.zeros([lat.size,lon.size]); SIG[:,:] = pl.float32('nan')

for i in range(lat.size):
    for j in range(lon.size):
        point = pvals[:,i,j]
        if pl.all(pl.isnan(point)) == True:
            continue
        else:
            sigcount = pl.where(point<0.05)
            if sigcount[0].size > 17.*(2./3.):
                    SIG[i,j] = 1

###############################################################################

if var == 'tas':
    meanlevs = pl.linspace(260,305,10); cmap_a = 'RdYlBu_r'
    sprdlevs = pl.linspace(0.5,5,10)
    biaslevs = pl.linspace(-5,5,11); cmap_cd = 'seismic'
    snlevs = pl.linspace(0.6,2,8); cmap_e = 'YlOrRd'
    units = 'K'
elif var == 'pr':
    #data_es_tm = (data_es_tm/data_em_tm)*100 # coefficient of variation
    meanlevs = pl.linspace(0,5,11); cmap_a = 'YlGnBu'
    units = 'mm day$^{-1}$'
    sprdlevs = pl.linspace(0,2,11)#pl.linspace(0.0,2.0,9)
    biaslevs = pl.linspace(-3,3,13); cmap_cd = 'seismic_r'
    snlevs = pl.linspace(-0.5,0.5,11); cmap_e = 'RdBu'
elif var == 'mrsos':
    units = 'm m$^{-3}$'
    meanlevs = pl.linspace(0,0.4,9); cmap_a = 'YlGnBu'
    sprdlevs = pl.linspace(0,0.2,9)
    biaslevs = pl.linspace(-0.3,0.3,13); cmap_c = 'seismic_r'
    snlevs = pl.linspace(-0.6,0.6,13); cmap_d = 'RdBu'

data_em_cyc, lon_cyc = util.add_cyclic_point(data_em_tm.data, coord=lon)

data_es_cyc = util.add_cyclic_point(data_es_tm.data)

#BIAS1_cyc = util.add_cyclic_point(BIAS1.data)
BIAS2_cyc = util.add_cyclic_point(BIAS2.data)

#CHECK1_cyc = util.add_cyclic_point(CHECK1)
CHECK2_cyc = util.add_cyclic_point(CHECK2)

SN_ratio_cyc = util.add_cyclic_point(SN_rat_mean)
SIG_cyc = util.add_cyclic_point(SIG)


###############################################################################
proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                                           '50m',edgecolor='grey',
                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')

#meanlevs = pl.linspace(0,5,11)#pl.linspace(0,0.4,9)#pl.linspace(260,305,10)#
#sprdlevs = pl.linspace(0.0,2.0,9)#pl.linspace(0,0.15,16)#pl.linspace(0.5,5,10)#
#biaslevs = pl.linspace(-3,3,13)
#snlevs = pl.linspace(-2,2,9)

###############################################################################
#fig, ax = pl.subplots(3,2,figsize=(10,9))
"""fig = pl.figure(figsize=(10,10))
gs = gridspec.GridSpec(3,4)
ig = [gs[0,:2],gs[0,2:],gs[1,:2],gs[1,2:],gs[2,1:3]]

ax1 = pl.subplot(ig[0],projection=proj,extent=ext) # cmip6 model
ax1.coastlines(linewidth=0.5,resolution='50m')
ax1.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax1.contourf(lon_cyc,lat.data,data_em_cyc,cmap=cmap_a,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=meanlevs,
                extend='both')

ax1.annotate(season+' CMIP6 '+var+' ensemble mean',(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)

GridLines(ax1,True,False,True,False)

f = pl.gcf()
colax = f.add_axes([0.07,0.67,0.42,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=12,labelpad=1)
cb.set_ticks(meanlevs)
cb.set_ticklabels(meanlevs.astype(int))
cb.ax.tick_params(labelsize=10,pad=2,length=0)

###############################################################################
ax2 = pl.subplot(ig[1],projection=proj,extent=ext) # 0.25 eobs
ax2.coastlines(linewidth=0.5,resolution='50m')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax2.contourf(lon_cyc,lat.data,data_es_cyc,cmap='plasma',
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=sprdlevs,
                extend='max')

ax2.annotate(season+' CMIP6 '+var+' ensemble spread',(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax2,True,False,False,True)

f = pl.gcf()
colax = f.add_axes([0.535,0.67,0.42,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=12,labelpad=1)
cb.set_ticks(sprdlevs)
cb.set_ticklabels(sprdlevs)
cb.ax.tick_params(labelsize=10,pad=2,length=0)
###############################################################################
ax3 = pl.subplot(ig[2],projection=proj,extent=ext) # 0.25 eobs
ax3.coastlines(linewidth=0.5,resolution='50m')
ax3.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax3.contourf(lon_cyc,lat.data,BIAS1_cyc,cmap=cmap_cd,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=biaslevs,
                extend='both')
ax3.contourf(lon_cyc,lat.data,CHECK1_cyc,colors='none',hatches='...')

ax3.annotate(season+' CMIP6 '+var+' minus E-OBS '+var2,(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax3,False,False,True,False)
###############################################################################
ax4 = pl.subplot(ig[3],projection=proj,extent=ext) # 0.25 eobs
ax4.coastlines(linewidth=0.5,resolution='50m')
ax4.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax4.contourf(lon_cyc,lat.data,BIAS2_cyc,cmap=cmap_cd,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=biaslevs,
                extend='both')
ax4.contourf(lon_cyc,lat.data,CHECK2_cyc,colors='none',hatches='...')

ax4.annotate(season+' CMIP6 '+var+' minus ERA5 '+var3,(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax4,False,False,False,True)

f = pl.gcf()
colax = f.add_axes([0.29,0.3365,0.42,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=12,labelpad=1)
cb.set_ticks(biaslevs)
cb.set_ticklabels(biaslevs.astype(int))
cb.ax.tick_params(labelsize=10,pad=2,direction='in',length=0)
###############################################################################
ax5 = pl.subplot(ig[4],projection=proj,extent=ext)
ax5.coastlines(linewidth=0.5,resolution='50m')
ax5.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax5.contourf(lon_cyc,lat.data,SN_ratio_cyc,cmap=cmap_e,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=snlevs,
                    extend='both')
ax5.contourf(lon_cyc,lat.data,SIG_cyc,colors='none',hatches='...')#levels=[1,2],
cn = ax5.contour(lon_cyc,lat.data,SN_ratio_cyc,levels=[1,2],colors=['k','k'])
ax5.clabel(cn,manual=True,fmt="%0.0f")

ax5.annotate(season+' CMIP6 ensemble mean '+var+' S/N',(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax5,False,True,True,False)

f = pl.gcf()
colax = f.add_axes([0.73,0.025,0.015,0.27])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)
cb.set_label('S/N',fontsize=12,labelpad=3)
cb.set_ticks(snlevs)
cb.set_ticklabels(snlevs)
cb.ax.tick_params(labelsize=10,pad=5,direction='in',length=0)
###############################################################################
#ax6 = pl.subplot(326,projection=proj,extent=ext)
#ax6.coastlines(linewidth=0.5,resolution='50m')
#ax6.add_feature(ocean_50m,alpha=1,zorder=5)
#
#cs = ax6.contourf(lon_cyc,lat.data,BIAS_djf_cyc,cmap='seismic',
#                 transform=ccrs.PlateCarree(),alpha=0.7,levels=biaslevs,
#                extend='both')
#
#ax6.annotate('DJF CMIP6 '+var+' minus E-OBS '+var2,(-13.8,68.8),
#             bbox={'facecolor':'w'},zorder=7)
#GridLines(ax6,False,True,False,False)

#f = pl.gcf()
#colax = f.add_axes([0.92,0.04,0.02,0.285])
#cb = pl.colorbar(cs,orientation='vertical',cax=colax)
#cb.set_label('K',fontsize=12,labelpad=1)
#cb.set_ticks(biaslevs)
#cb.set_ticklabels(biaslevs.astype(int))
#cb.ax.tick_params(labelsize=10)
###############################################################################

pl.subplots_adjust(top=0.99,bottom=0.00,left=0.07,right=0.95,hspace=0.04,
                   wspace=0.29)

fig.text(0.005,0.95,'(a)',size=11)
fig.text(0.505,0.95,'(b)',size=11)
fig.text(0.005,0.615,'(c)',size=11)
fig.text(0.505,0.615,'(d)',size=11)
fig.text(0.235,0.285,'(e)',size=11)
#fig.text(0.48,0.31,'(f)',size=11)

#pl.savefig(indecis+'figures/'+var+'_ens_stats_bias_'+season.lower()+'_1950-2014.png',dpi=360)"""

###############################################################################

fig, ax = pl.subplots(2,2,figsize=(10,7))

ax1 = pl.subplot(221,projection=proj,extent=ext) # cmip6 model
ax1.coastlines(linewidth=0.5,resolution='50m')
ax1.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax1.contourf(lon_cyc,lat.data,data_em_cyc,cmap=cmap_a,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=meanlevs,
                extend='max')

ax1.annotate(season+' CMIP6 '+var+' ensemble mean',(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)

GridLines(ax1,True,False,True,False)

f = pl.gcf()
colax = f.add_axes([0.08,0.54,0.41,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=11,labelpad=1)
cb.set_ticks(meanlevs)
cb.set_ticklabels(meanlevs)
cb.ax.tick_params(labelsize=10,pad=2,length=0)
###############################################################################

ax2 = pl.subplot(222,projection=proj,extent=ext) # 0.25 eobs
ax2.coastlines(linewidth=0.5,resolution='50m')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax2.contourf(lon_cyc,lat.data,data_es_cyc,cmap='plasma',
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=sprdlevs,
                extend='max')

ax2.annotate(season+' CMIP6 '+var+' ensemble spread',(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax2,True,False,False,True)

f = pl.gcf()
colax = f.add_axes([0.535,0.54,0.41,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=11,labelpad=1)
cb.set_ticks(sprdlevs[::2])
cb.set_ticklabels(sprdlevs[::2])
cb.ax.tick_params(labelsize=10,pad=2,length=0)
###############################################################################

ax3 = pl.subplot(223,projection=proj,extent=ext) # 0.25 eobs
ax3.coastlines(linewidth=0.5,resolution='50m')
ax3.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax3.contourf(lon_cyc,lat.data,BIAS2_cyc,cmap=cmap_c,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=biaslevs,
                extend='both')
ax3.contourf(lon_cyc,lat.data,CHECK2_cyc,colors='none',hatches='...')

ax3.annotate(season+' CMIP6 '+var+' minus ERA5 '+var3,(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax3,False,False,True,False)

f = pl.gcf()
colax = f.add_axes([0.08,0.06,0.41,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label(units,fontsize=11,labelpad=1)
cb.set_ticks(biaslevs[::2])
cb.set_ticklabels(biaslevs[::2])
cb.ax.tick_params(labelsize=10,pad=2,direction='in',length=0)
###############################################################################

ax4 = pl.subplot(224,projection=proj,extent=ext)
ax4.coastlines(linewidth=0.5,resolution='50m')
ax4.add_feature(ocean_50m,alpha=1,zorder=5)

cs = ax4.contourf(lon_cyc,lat.data,SN_ratio_cyc,cmap=cmap_d,
                 transform=ccrs.PlateCarree(),alpha=0.7,levels=snlevs,
                    extend='both')
ax4.contourf(lon_cyc,lat.data,SIG_cyc,colors='none',hatches='...')#,levels=[1,2])
cn = ax4.contour(lon_cyc,lat.data,SN_ratio_cyc,levels=[-1,1],colors=['k','k'])
#ax4.clabel(cn,manual=True,fmt="%0.1f")

ax4.annotate(season+' CMIP6 ensemble mean '+var+' S/N',(-13.8,68.8),
             bbox={'facecolor':'w'},zorder=7)
GridLines(ax4,False,False,False,True)

f = pl.gcf()
colax = f.add_axes([0.53,0.06,0.41,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
#cb.set_label('K',fontsize=12,labelpad=1)
cb.set_ticks(snlevs[::2])
cb.set_ticklabels(snlevs[::2])
cb.ax.tick_params(labelsize=10)

###############################################################################

pl.subplots_adjust(top=0.98,bottom=0.07,left=0.07,right=0.95,hspace=0.10,
                   wspace=0.11)

fig.text(0.005,0.95,'(a)',size=11)
fig.text(0.505,0.95,'(b)',size=11)
fig.text(0.005,0.46,'(c)',size=11)
fig.text(0.505,0.46,'(d)',size=11)

#pl.savefig(indecis+'figures/'+var+'_ens_stats_bias_'+season.lower()+'_1950-2014.png',dpi=360)