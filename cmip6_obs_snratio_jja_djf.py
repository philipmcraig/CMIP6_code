# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:37:31 2022

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
    gl.xlabel_style = {'color': 'k','size':9}#10
    gl.ylabel_style = {'color': 'k','size':9}#10
    
    return None

def get_cmapBoundary(x, base=2):
	""" Function to round up to the nearest 2 (or other defined base) if the number is > 0 or down to 
		the nearest 2 (or other defined base) if the number is < 0.

		Inputs:
		 x - float or integer, containing the value to be rounded.
		 base - integer, containing the base to perform the rounding to.

		Outputs:
	"""
	if x>0: return int(base * pl.ceil(float(x)/base))
	else: return int(base * pl.floor(float(x)/base))


def get_eofColormap(eof):
    """ Function to get a diverging matplotlib colormap that is centred at zero, and has upper and lower 
    limits defined by the EOF to be plotted. Maximum upper limit is 10 and maximum lower limit is 
    -10.

    Inputs:
    eof - a 1-D NumPy array, containing the normalised EOF to be plotted.

    Outputs:
    cmap - a matplotlib color map, containing the user-defined color map to use when plotting.
    clevs - a 1-D NumPy array, containing the values associated with the colormap.
    """
    # 20 colour blue-white-red colour scale.
    colors=[(0/255.,0/255.,96/255.),\
        (0/255.,0/255.,128/255.),\
        (0/255.,0/255.,191.5/255.),\
        (0/255.,0/255.,255/255.),\
        (0/255.,64/255.,255/255.),\
        (0/255.,128/255.,255/255.),\
        (0/255.,191.5/255.,255/255.),\
        #(0/255.,225/255.,255/255.),\
        (191.5/255.,255/255.,255/255.),\
        #(255/255.,255/255.,255/255.),\
        #(255/255.,255/255.,255/255.),\
        (255/255.,255/255.,191.5/255.),\
        #(255/255.,225/255.,0/255.),\
        (255/255.,191.5/255.,0/255.),\
        (255/255.,128/255.,0/255.),\
        (255/255.,64/255.,0/255.),\
        (255/255.,0/255.,0/255.), \
        (192/255.,0/255.,0/255.),\
        (160/255.,0/255.,0/255.),\
        (128/255.,0/255.,0/255.)]

    # add end colors for extension of colour map
    #colors.insert(0,'midnightblue')
    #colors.append('maroon')

    # get levels and limits of colour scale
    bounds=pl.linspace(-2,2,17)
    #bounds = [-100000, -85000, -70000, -55000, -40000, -25000, -10000, -2000, 2000,
    #      10000, 25000, 40000, 55000, 70000, 85000, 100000]
    upperlim=pc.NearestIndex(bounds,pl.nanmax(eof))#pl.where(bounds==get_cmapBoundary(pl.nanmax(eof)))[0]
    lowerlim=pc.NearestIndex(bounds,pl.nanmin(eof))#pl.where(bounds==get_cmapBoundary(pl.nanmin(eof)))[0]
    clevs=bounds#[lowerlim:upperlim+1]

    # create colour map
    cmap=pl.matplotlib.colors.ListedColormap(colors)
    cmap.set_over(color=(96/255.,0/255.,0/255.))
    cmap.set_under(color=(0/255.,0/255.,64/255.))

    return cmap, clevs

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
era5dir = ncasdir + 'ERA5/'
indecis = ncasdir + 'INDECIS/'
cmipdir = ncasdir + 'CMIP6/'
eobsdir = ncasdir + 'EOBS/'

var = 'tas'
var2 = 'tg'
var3 = 't2m'
seasons = ['JJA','DJF']

meanfile = xr.open_dataset(cmipdir + 'standard_grid/'+var+'_cmip6_ensmean.nc')
lat = xr.DataArray(meanfile.lat)
lon = xr.DataArray(meanfile.lon)
time = xr.DataArray(meanfile.time)
cmipmean = xr.DataArray(getattr(meanfile,var))
meanfile.close()

allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014_new/*')
tasfiles = glob.glob(cmipdir+'global_means/*')

modeldata = pl.zeros([len(allfiles),len(time),lat.size,lon.size])

gmst = pl.zeros([len(tasfiles),len(time)])

modnames = pl.zeros([len(allfiles)],dtype='object')

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var))
    ncfile.close
    
    tasfile = xr.open_dataset(tasfiles[nci])
    gmst[nci,:] = pl.squeeze(xr.DataArray(tasfile.tas))
    tasfile.close()
    
    split1 = allfiles[nci].split('/')
    split2 = split1[-1].split('_')[2]
    
    modnames[nci] = split2

era5file = xr.open_dataset(era5dir+'era5_land_seasmean_1950-2014.nc')
era5data = xr.DataArray(getattr(era5file,var3))
era5time = xr.DataArray(era5file.time)
era5lat = xr.DataArray(era5file.latitude)
era5lon = xr.DataArray(era5file.longitude)
era5file.close()

nce5_gmst = xr.open_dataset(era5dir+'era5_t2m_seasmean_1950-2014_gm.nc')
gmst_e5 = xr.DataArray(nce5_gmst.t2m)
nce5_gmst.close()

eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_0.25deg_reg_v23.0e.nc')
eobsdata = xr.DataArray(getattr(eobsfile,var2))
eobstime = xr.DataArray(eobsfile.time)
eobslat = xr.DataArray(eobsfile.latitude)
eobslon = xr.DataArray(eobsfile.longitude)
eobsfile.close()

nceo_gmst = xr.open_dataset(eobsdir+'tg_seasmean_remapbil1.0_v23.0e_gm.nc')
gmst_eo = xr.DataArray(nceo_gmst.tg)
nceo_gmst.close()

SN_ratios_cmip = pl.zeros([len(seasons),len(allfiles),lat.size,lon.size])
SN_rat_mean = pl.zeros([len(seasons),lat.size,lon.size])
pvals_cmip = SN_ratios_cmip.copy()
SIG_cmip = pl.zeros([2,lat.size,lon.size]); SIG_cmip[:,:,:] = pl.float16('nan')
#
SN_ratios_era5 = pl.zeros([len(seasons),era5lat.size,era5lon.size])
pvals_era5 = SN_ratios_era5.copy()
SIG_era5 = pl.zeros([2,era5lat.size,era5lon.size]); SIG_era5[:,:,:] = pl.float16('nan')
#
SN_ratios_eobs = pl.zeros([len(seasons),eobslat.size,eobslon.size])
pvals_eobs = SN_ratios_eobs.copy()
SIG_eobs = pl.zeros([2,eobslat.size,eobslon.size]); SIG_eobs[:,:,:] = pl.float16('nan')
#
for i in range(len(seasons)):
    if seasons[i] == 'JJA':
        data_em_tm = pl.mean(cmipmean[1::4,:,:],axis=0)#*86400####/(0.1*1000)
        eobs_mn = pl.mean(eobsdata[1:-1:4],axis=0) + 273.15
        era5_mn = pl.mean(era5data[2:-1:4],axis=0)#*1000
        gmst_ts = gmst[:,1::4]
        gmst_e5_ssn = gmst_e5.data[2:-1:4]
        gmst_eo_ssn = gmst_eo[1:-1:4]
        data_ssn = modeldata[:,1::4]#*86400####/(0.1*1000)
        data_e5 = era5data[2:-1:4]#*1000
        data_eo = eobsdata[1:-1:4] + 273.15
    elif seasons[i] == 'DJF':
        data_em_tm = pl.mean(cmipmean[3::4,:,:],axis=0)#*86400####/(0.1*1000)
        eobs_mn = pl.mean(eobsdata[3:-1:4],axis=0) + 273.15
        era5_mn = pl.mean(era5data[4:-1:4],axis=0)#*1000
        gmst_ts = gmst[:,3::4]
        gmst_e5_ssn = gmst_e5.data[4:-1:4]
        gmst_eo_ssn = gmst_eo[3:-1:4]
        data_ssn = modeldata[:,3::4]#*86400##/(0.1*1000)
        data_e5 = era5data[4:-1:4]#*1000
        data_eo = eobsdata[3:-1:4] + 273.15
    
    for j in range(len(allfiles)):
        # smooth GMST with 15 years low-pass filter
        # use Pandas rolling average
        df = pd.DataFrame(gmst_ts[j])
        gmst_smth = df.rolling(15,min_periods=1,center=True).mean()
        gmst_smth = pl.squeeze(pl.asarray(gmst_smth))
        
        SN_ratios_cmip[i,j], pvals_cmip[i,j] = S_N_ratio(data_ssn[j],lat,lon,
                                                                    gmst_smth)

    SN_rat_mean[i] = pl.nanmean(SN_ratios_cmip[i],axis=0)
    
#
    df_e5 = pd.DataFrame(pl.squeeze(gmst_e5_ssn))
    gmst_e5_smth = df_e5.rolling(15,min_periods=1,center=True).mean()
    gmst_e5_smth = pl.squeeze(pl.asarray(gmst_e5_smth))
    
    df_eo = pd.DataFrame(pl.squeeze(gmst_eo_ssn))
    gmst_eo_smth = df_eo.rolling(15,min_periods=1,center=True).mean()
    gmst_eo_smth = pl.squeeze(pl.asarray(gmst_eo_smth))
#    
    SN_ratios_era5[i], pvals_era5[i] = S_N_ratio(data_e5,era5lat,era5lon,
                                                                gmst_e5_smth)
#    
    SN_ratios_eobs[i], pvals_eobs[i] = S_N_ratio(data_eo,eobslat,eobslon,
                                                                gmst_eo_smth)
#    
#    
    for phi in range(lat.size):
        for lam in range(lon.size):
            point1 = pvals_cmip[i,:,phi,lam]
            if pl.all(pl.isnan(point1)) == True:
                continue
            else:
                sigcount1 = pl.where(point1<0.05)
                if sigcount1[0].size > 17.*(2./3.):
                        SIG_cmip[i,phi,lam] = 1
            
    for phi in range(era5lat.size):
        for lam in range(era5lon.size):
            point2 = pvals_era5[i,phi,lam]
            if pl.isnan(point2) == True:
                continue
            else:
                if point2 < 0.05:
                    SIG_era5[i,phi,lam] = 1

    for phi in range(eobslat.size):
        for lam in range(eobslon.size):
            point3 = pvals_eobs[i,phi,lam]
            if pl.isnan(point3) == True:
                continue
            else:
                if point3 < 0.05:
                    SIG_eobs[i,phi,lam] = 1

spread = pl.nanstd(SN_ratios_cmip,axis=1)

###############################################################################
###############################################################################

SN_cmip_cyc, lon_cyc = util.add_cyclic_point(SN_rat_mean, coord=lon)
SIG_cmip_cyc = util.add_cyclic_point(SIG_cmip)
SN_ratios_cyc = util.add_cyclic_point(SN_ratios_cmip)
spread_cyc = util.add_cyclic_point(spread)

SN_era5_cyc, e5lon_cyc = util.add_cyclic_point(SN_ratios_era5,coord=era5lon)
SIG_era5_cyc = util.add_cyclic_point(SIG_era5)

SN_eobs_cyc, eolon_cyc = util.add_cyclic_point(SN_ratios_eobs,coord=eobslon)
SIG_eobs_cyc = util.add_cyclic_point(SIG_eobs)

snlevs = pl.linspace(-2,2,17)
cmap, clevs = get_eofColormap(SN_cmip_cyc)#'seismic_r'#'RdYlBu_r'#'BrBG'#
alpha = 0.6
###############################################################################
###############################################################################


fig, ax = pl.subplots(3,2,figsize=(10,9)) # (10,9) for 6 panels, (8,5) for 4 panels
proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
#borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
#                                           '50m',edgecolor='grey',
#                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')

###############################################################################

ax1 = pl.subplot(321,projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m')
ax1.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax1.contourf(lon_cyc,lat.data,SN_cmip_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
                #norm=pl.Normalize(-2,2))
ax1.contourf(lon_cyc,lat.data,SIG_cmip_cyc[0],colors='none',hatches='....')
#ax1.contour(lon_cyc,lat.data,SN_cmip_cyc[0],colors='k',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
#ax1.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax1,True,False,True,False)
ax1.text(-26,69,'a',size=14)
ax1.annotate('JJA CMIP6 '+var+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax2 = pl.subplot(322,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)
ax2.contourf(lon_cyc,lat.data,SN_cmip_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
                #norm=pl.Normalize(-2,2))
ax2.contourf(lon_cyc,lat.data,SIG_cmip_cyc[1],colors='none',hatches='....')
#ax2.contour(lon_cyc,lat.data,SN_cmip_cyc[1],colors='k',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
#ax2.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax2,True,False,False,True)
ax2.text(-20,69,'b',size=14)
ax2.annotate('DJF CMIP6 '+var+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax3 = pl.subplot(323,projection=proj,extent=ext)
ax3.coastlines(linewidth=0.5,resolution='50m')
ax3.add_feature(ocean_50m,alpha=1,zorder=5)
ax3.contourf(e5lon_cyc,era5lat.data,SN_era5_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax3.contourf(e5lon_cyc,era5lat.data,SIG_era5_cyc[0],colors='none',hatches='....')
#ax3.contour(e5lon_cyc,era5lat.data,SN_era5_cyc[0],colors='k',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
cn = ax3.contour(e5lon_cyc,era5lat.data,SN_era5_cyc[0],colors='k',
            levels=[-3,3],linewidths=0.75)
#ax3.clabel(cn,manual=True,fmt="%0.0f")

GridLines(ax3,False,False,True,False)
ax3.text(-26,69,'c',size=14)
ax3.annotate('JJA ERA5 '+var3+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax4 = pl.subplot(324,projection=proj,extent=ext)
ax4.coastlines(linewidth=0.5,resolution='50m')
ax4.add_feature(ocean_50m,alpha=1,zorder=5)
ax4.contourf(e5lon_cyc,era5lat.data,SN_era5_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax4.contourf(e5lon_cyc,era5lat.data,SIG_era5_cyc[1],colors='none',hatches='....')
#ax4.contour(e5lon_cyc,era5lat.data,SN_era5_cyc[1],colors='k',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
#ax4.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax4,False,False,False,True)
ax4.text(-20,69,'d',size=14)
ax4.annotate('DJF ERA5 '+var3+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax5 = pl.subplot(325,projection=proj,extent=ext)
ax5.coastlines(linewidth=0.5,resolution='50m')
ax5.add_feature(ocean_50m,alpha=1,zorder=5)
ax5.contourf(eolon_cyc,eobslat.data,SN_eobs_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax5.contourf(eolon_cyc,eobslat.data,SIG_eobs_cyc[0],colors='none',hatches='....')
#ax5.contour(eolon_cyc,eobslat.data,SN_eobs_cyc[0],colors='k',
#            levels=[-3,-2,-1,0,1,2,3],linewidths=0.75)
#ax5.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax5,False,False,True,False)
ax5.text(-26,69,'e',size=14)
ax5.annotate('JJA E-OBS '+var2+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

##############################################################################

ax6 = pl.subplot(326,projection=proj,extent=ext)
ax6.coastlines(linewidth=0.5,resolution='50m')
ax6.add_feature(ocean_50m,alpha=1,zorder=5)
ax6.contourf(eolon_cyc,eobslat.data,SN_eobs_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax6.contourf(eolon_cyc,eobslat.data,SIG_eobs_cyc[1],colors='none',hatches='....')
#ax6.contour(eolon_cyc,eobslat.data,SN_eobs_cyc[1],colors='k',
#            levels=[-3,-2,-1,0,1,2,3],linewidths=0.75)
ax6.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax6,False,False,False,True)
ax6.text(-20,69,'f',size=14)
ax6.annotate('DJF E-OBS '+var2+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

f = pl.gcf()
colax = f.add_axes([0.1,0.06,0.8,0.02])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
#cb.set_label('K',fontsize=12,labelpad=1)
cb.set_ticks(snlevs[::2])
cb.set_ticklabels(snlevs[::2])
cb.ax.tick_params(labelsize=12)

pl.tight_layout()
pl.subplots_adjust(top=0.96,bottom=0.10,hspace=0.04,left=0.0,right=1.0,
                   wspace=-0.15)

#pl.savefig(indecis+'figures/'+var+'_cmip6_era5_eobs_sn_jja_djf.png',dpi=400)
#pl.savefig(indecis+'figures/'+var+'_cmip6_era5_eobs_sn_jja_djf.pdf',dpi=400)

"""fig, ax = pl.subplots(2,2,figsize=(8,5))

ax1 = pl.subplot(221,projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m')
ax1.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax1.contourf(lon_cyc,lat.data,SN_cmip_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax1.contourf(lon_cyc,lat.data,SIG_cmip_cyc[0],colors='none',hatches='....')
#ax1.contour(lon_cyc,lat.data,SN_cmip_cyc[0],colors='k',
#           levels=[-2,-1,0,1,2],linewidths=0.75)

GridLines(ax1,True,False,True,False)
ax1.text(-27,69,'a',size=14)
ax1.annotate('JJA CMIP6 '+var+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax2 = pl.subplot(222,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)
ax2.contourf(lon_cyc,lat.data,SN_cmip_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax2.contourf(lon_cyc,lat.data,SIG_cmip_cyc[1],colors='none',hatches='....')
#ax2.contour(lon_cyc,lat.data,SN_cmip_cyc[1],colors='k',linestyles='-',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
#ax2.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax2,True,False,False,True)
ax2.text(-20,69,'b',size=14)
ax2.annotate('DJF CMIP6 '+var+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax3 = pl.subplot(223,projection=proj,extent=ext)
ax3.coastlines(linewidth=0.5,resolution='50m')
ax3.add_feature(ocean_50m,alpha=1,zorder=5)
ax3.contourf(e5lon_cyc,era5lat.data,SN_era5_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax3.contourf(e5lon_cyc,era5lat.data,SIG_era5_cyc[0],colors='none',hatches='....')
#ax3.contour(e5lon_cyc,era5lat.data,SN_era5_cyc[0],colors='k',linestyles='-',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
#cn = ax3.contour(e5lon_cyc,era5lat.data,SN_era5_cyc[0],colors='k',linestyles='-',
#            levels=[-3,3],linewidths=0.75)
#ax3.clabel(cn,manual=True,fmt="%0.0f")

GridLines(ax3,False,False,True,False)
ax3.text(-27,69,'c',size=14)
ax3.annotate('JJA ERA5 '+var3+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax4 = pl.subplot(224,projection=proj,extent=ext)
ax4.coastlines(linewidth=0.5,resolution='50m')
ax4.add_feature(ocean_50m,alpha=1,zorder=5)
ax4.contourf(e5lon_cyc,era5lat.data,SN_era5_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)
ax4.contourf(e5lon_cyc,era5lat.data,SIG_era5_cyc[1],colors='none',hatches='....')
#ax4.contour(e5lon_cyc,era5lat.data,SN_era5_cyc[1],colors='k',linestyles='-',
#            levels=[-2,-1,0,1,2],linewidths=0.75)
#ax4.clabel(cn,manual=False,fmt="%0.0f")

GridLines(ax4,False,False,False,True)
ax4.text(-20,69,'d',size=14)
ax4.annotate('DJF ERA5 '+var3+' S/N',(-13.8,69),bbox={'facecolor':'w'},zorder=20)

###############################################################################

f = pl.gcf()
colax = f.add_axes([0.1,0.06,0.8,0.02])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_ticks(snlevs[::2])
cb.set_ticklabels(snlevs[::2])
cb.ax.tick_params(labelsize=12)

pl.subplots_adjust(top=0.94,left=0.08,right=0.96,hspace=0.07,wspace=0.04)"""

#pl.savefig(indecis+'figures/'+var+'_cmip6_era5_sn_jja_djf.png',dpi=400)
#pl.savefig(indecis+'figures/'+var+'_cmip6_era5_sn_jja_djf.pdf',dpi=400)

###############################################################################

"""fig, ax = pl.subplots(6,3,figsize=(8.0,9.5))

for i in range(SN_ratios_cmip.shape[1]):
    axx = pl.subplot(6,3,i+1,projection=proj,extent=ext)
    axx.coastlines(linewidth=0.5,resolution='50m')
    axx.add_feature(ocean_50m,alpha=1,zorder=5)
    cs = axx.contourf(lon_cyc,lat.data,SN_ratios_cyc[0,i],cmap=cmap,extend='both',
                 alpha=alpha,transform=ccrs.PlateCarree(),levels=snlevs)

    axx.text(-13.8,68.8,modnames[i],{'fontsize':8,'bbox':dict(facecolor='w')},
             zorder=7)

    if i == 0:
        GridLines(axx,True,False,True,False)
    elif i == 1:
        GridLines(axx,True,False,False,False)
    elif i == 2:
        GridLines(axx,True,False,False,True)
    elif i in (3,6,9,12):
        GridLines(axx,False,False,True,False)
    elif i in (4,7,10,13):
        GridLines(axx,False,False,False,False)
    elif i in (5,8,11,14):
        GridLines(axx,False,False,False,True)
    elif i == 15:
        GridLines(axx,False,False,True,False)
    elif i == 16:
        GridLines(axx,False,False,False,False)

f = pl.gcf()
colax = f.add_axes([0.06,0.020,0.57,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_ticks(snlevs[::2])
cb.set_ticklabels(snlevs[::2])
cb.ax.tick_params(labelsize=9,pad=2,direction='in',length=3)

ax18 = pl.subplot(6,3,18,projection=proj,extent=ext)
ax18.coastlines(linewidth=0.5,resolution='50m')
ax18.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax18.contourf(lon_cyc,lat.data,spread_cyc[0],levels=pl.linspace(0.1,0.6,6),
                   cmap='plasma',alpha=alpha,transform=ccrs.PlateCarree(),
                    extend='both')

ax18.text(-13.8,69.2,var+' ensemble spread',
          {'fontsize':9,'bbox':dict(facecolor='w',edgecolor='none',pad=0.9)},
             zorder=7)
GridLines(ax18,False,False,False,True)

f = pl.gcf()
colax = f.add_axes([0.68,0.020,0.26,0.015])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
#cb.set_ticks(snlevs[::2])
#cb.set_ticklabels(snlevs[::2])
cb.ax.tick_params(labelsize=8.5,pad=2,direction='in',length=3)
#ax18.text(36,30,'mm day$^{-1}$',{'fontsize':10})

pl.tight_layout()
pl.subplots_adjust(bottom=0.04,left=0.02,right=0.98,hspace=0.1,wspace=-0.1,
                   top=0.97)"""

#pl.savefig(indecis+'figures/cmip6_sn_models_sprd_'+var+'_jja.png',dpi=400)
#pl.savefig(indecis+'figures/cmip6_sn_models_sprd_'+var+'_jja.pdf',dpi=400)"""