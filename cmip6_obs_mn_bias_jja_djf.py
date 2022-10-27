# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 12:26:05 2022

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
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mplPath
import pcraig_funcs as pc

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
        #(0/255.,0/255.,128/255.),\
        (0/255.,0/255.,192/255.),\
        #(0/255.,0/255.,255/255.),\
        (0/255.,64/255.,255/255.),\
        #(0/255.,128/255.,255/255.),\
        (0/255.,191.5/255.,255/255.),\
        #(0/255.,225/255.,255/255.),\
        (191.5/255.,255/255.,255/255.),\
        #(255/255.,255/255.,255/255.),\
        #(255/255.,255/255.,255/255.),\
        (255/255.,255/255.,191.5/255.),\
        #(255/255.,225/255.,0/255.),\
        (255/255.,191.5/255.,0/255.),\
        #(255/255.,128/255.,0/255.),\
        (255/255.,64/255.,0/255.),\
        #(255/255.,0/255.,0/255.), \
        (192/255.,0/255.,0/255.),\
        #(160/255.,0/255.,0/255.),\
        (128/255.,0/255.,0/255.)]

    # add end colors for extension of colour map
    #colors.insert(0,'midnightblue')
    #colors.append('maroon')

    # get levels and limits of colour scale
    bounds=pl.linspace(-5,5,11)
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

def CheckBiasSign(cmipdir,BIAS1,BIAS2,eobs_seasmean,era5_seasmean,season,
                                                          var,var2,var3):
    """
    """
    allfiles = glob.glob(cmipdir+'standard_grid/'+var+'_1950-2014/*')
    
    #BIAS_mean = cmip_tm - eobs_mn
    BIAS_mods1 = pl.zeros([len(allfiles),BIAS2.shape[0],BIAS2.shape[1]])
    BIAS_mods2 = BIAS_mods1.copy()
    #runtot = pl.zeros_like(BIAS_mean) # running total
    CHECK1 = pl.zeros([BIAS2.shape[0],BIAS2.shape[1]])
    CHECK1[:,:] = pl.float32('nan')
    CHECK2 = CHECK1.copy()
    
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
        
        BIAS_mods1[m] = model_tm - eobs_seasmean
        BIAS_mods2[m] = model_tm - era5_seasmean
    
    for i in range(BIAS_mods2.shape[1]):
        for j in range(BIAS_mods2.shape[2]):
            point1 = BIAS_mods1[:,i,j]
            if pl.all(pl.isnan(point1)) == True:
                continue
            else:
                meansign1 = pl.sign(BIAS1[i,j])
                #print i,j
                meansign1 = int(meansign1)
                SIGNS1 = pl.sign(point1)
                samesign1 = pl.where(SIGNS1==meansign1)
                
                #pos_count1 = len(list(filter(lambda x: (x > 0), point1)))
                #neg_count1 = len(list(filter(lambda x: (x < 0), point1)))
                #zer_count1 = len(list(filter(lambda x: (x == 0), point1)))
                
                #A1 = pl.array([pos_count1,neg_count1,zer_count1])
                
                if samesign1[0].size > 17.*(2./3.):
                    CHECK1[i,j] = 1.
            
            ###################################################################
            point2 = BIAS_mods2[:,i,j]
            if pl.all(pl.isnan(point2)) == True:
                continue
            else:
                meansign2 = pl.sign(BIAS2[i,j])
                meansign2 = int(meansign2)
                SIGNS2 = pl.sign(point2)
                samesign2 = pl.where(SIGNS2==meansign2)
#                pos_count2 = len(list(filter(lambda x: (x > 0), point2)))
#                neg_count2 = len(list(filter(lambda x: (x < 0), point2)))
#                zer_count2 = len(list(filter(lambda x: (x == 0), point2)))
#                
#                A2 = pl.array([pos_count2,neg_count2,zer_count2])
                
                if samesign2[0].size > 17.*(2./3.):
                    CHECK2[i,j] = 1.
    
    return CHECK1, CHECK2

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

modeldata = pl.zeros([len(allfiles),len(time),lat.size,lon.size])

for nci in range(len(allfiles)):
    ncfile = xr.open_dataset(allfiles[nci])
    modeldata[nci,:,:,:] = xr.DataArray(getattr(ncfile,var))
    ncfile.close

    split1 = allfiles[nci].split('/')
    split2 = split1[-1].split('_')[2]

era5file = xr.open_dataset(era5dir+'era5_surf_seasmean_s1.5.nc')
era5data = xr.DataArray(getattr(era5file,var3))
era5time = xr.DataArray(era5file.time)
era5lat = xr.DataArray(era5file.lat)
era5lon = xr.DataArray(era5file.lon)
era5file.close()

eobsfile = xr.open_dataset(eobsdir+var2+'_seasmean_remapbil1.0_v23.0e_s1.5.nc')
eobsdata = xr.DataArray(getattr(eobsfile,var2))
eobstime = xr.DataArray(eobsfile.time)
eobslat = xr.DataArray(eobsfile.lat)
eobslon = xr.DataArray(eobsfile.lon)
eobsfile.close()

CMIP_ssn = pl.zeros([len(seasons),lat.size,lon.size])

BIAS_era5 = pl.zeros([len(seasons),lat.size,lon.size])
CHECK_e5 = BIAS_era5.copy()
#BIAS_era5_mn = pl.zeros_like([len(seasons),lat.size,lon.size])

BIAS_eobs = pl.zeros([len(seasons),lat.size,lon.size])
CHECK_eo = BIAS_eobs.copy()
#BIAS_eobs_mn = pl.zeros([len(seasons),lat.size,lon.size])

for i in range(len(seasons)):
    if seasons[i] == 'JJA':
        data_em_tm = pl.mean(cmipmean[1::4,:,:],axis=0)#/(0.1*1000)#*86400#
        eobs_mn = pl.mean(eobsdata[1::4],axis=0) + 273.15
        era5_mn = pl.mean(era5data[2:-13:4],axis=0)#*1000
        #gmst_ts = gmst[:,1::4]
        #gmst_e5_ssn = gmst_e5.data[2:-1:4]
        #gmst_eo_ssn = gmst_eo[1:-1:4]
        data_ssn = modeldata[:,1::4]/(0.1*1000)#*86400#
        data_e5 = era5data[2:-13:4]
        data_eo = eobsdata[1:-1:4] + 273.15#*1000
    elif seasons[i] == 'DJF':
        data_em_tm = pl.mean(cmipmean[3::4,:,:],axis=0)#/(0.1*1000)#*86400#
        eobs_mn = pl.mean(eobsdata[3::4],axis=0) + 273.15
        era5_mn = pl.mean(era5data[4:-13:4],axis=0)#*1000
        #gmst_ts = gmst[:,3::4]
        #gmst_e5_ssn = gmst_e5.data[4:-1:4]
        #gmst_eo_ssn = gmst_eo[3:-1:4]
        data_ssn = modeldata[:,3::4]#/(0.1*1000)#*86400#
        data_e5 = era5data[4:-13:4]
        data_eo = eobsdata[3:-1:4] + 273.15

    CMIP_ssn[i] = data_em_tm

    #for j in range(len(allfiles)):
        #model_tm = pl.mean(data_ssn[j],axis=0)
        #e5_tm = pl.mean(data_e5,axis=0)
        
    BIAS_era5[i] = data_em_tm - era5_mn
    BIAS_eobs[i] = data_em_tm - eobs_mn
        
    
    CHECK_eo[i],CHECK_e5[i] = CheckBiasSign(cmipdir,BIAS_eobs[i],BIAS_era5[i],
                                    eobs_mn,era5_mn,seasons[i],var,var2,var3)
    #CHECK_e5[i] = CheckBiasSign(cmipdir,BIAS_era5[i],era5_mn,seasons[i],var,var3)

#BIAS_era5_mn = pl.mean(BIAS_era5,axis=1)
#BIAS_eobs_mn = pl.mean(BIAS_eobs,axis=1)

###############################################################################

nrth_eur = [(-5,48),(-11,51.8),(-11,58),(16.7,71),(28.9,71),(28.9,48)]
sth_eur = [(-11,36),(-11,43.8),(-5,48),(28.9,48),(28.6,41),(25.1,41),(25.1,36),
           (15,36),(10.3,38.1),(2.7,38.1),(-5.6,36)]


###############################################################################
CMIP_ssn_cyc, lon_cyc = util.add_cyclic_point(CMIP_ssn,lon)
BIAS_e5_cyc = util.add_cyclic_point(BIAS_era5)
BIAS_eo_cyc = util.add_cyclic_point(BIAS_eobs)

CHECK_e5_cyc = util.add_cyclic_point(CHECK_e5)
CHECK_eo_cyc = util.add_cyclic_point(CHECK_eo)

biaslevs = pl.linspace(-5,5,11)
meanlevs = pl.linspace(260,305,10)
cmap, clevs = get_eofColormap(BIAS_e5_cyc)
alpha=0.6

###############################################################################

fig, ax = pl.subplots(2,2,figsize=(10,9)) # (10,9) for 6 panels, (8,5) for 4 panels
proj = ccrs.PlateCarree()
ext = [-15,42,35,70]
borders_50m = cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                                           '50m',edgecolor='grey',
                                        facecolor='none')
ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='none',
                                        facecolor='w')

###############################################################################

ax1 = pl.subplot(321,projection=proj,extent=ext)
ax1.coastlines(linewidth=0.5,resolution='50m')
ax1.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax1.contourf(lon_cyc,lat.data,CMIP_ssn_cyc[0],cmap='RdYlBu_r',extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=meanlevs)

GridLines(ax1,True,False,True,False)
ax1.text(-24,69,'a',size=14)
ax1.annotate('JJA CMIP6 '+var+' ensemble mean',(-13.8,69),
                                             bbox={'facecolor':'w'},zorder=20)


poly_nrth = Polygon(nrth_eur,fc='none')
pn = PatchCollection([poly_nrth],zorder=6)
pn.set_edgecolor('k')
pn.set_facecolor('none')
pn.set_linewidth(2)
ax1.add_collection(pn)

poly_sth = Polygon(sth_eur,fc='none')
ps = PatchCollection([poly_sth],zorder=6)
ps.set_edgecolor('k')
ps.set_facecolor('none')
ps.set_linewidth(2)
ax1.add_collection(ps)
###############################################################################

ax2 = pl.subplot(322,projection=proj,extent=ext)
ax2.coastlines(linewidth=0.5,resolution='50m')
ax2.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax2.contourf(lon_cyc,lat.data,CMIP_ssn_cyc[1],cmap='RdYlBu_r',extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=meanlevs)

GridLines(ax2,True,False,False,True)
ax2.text(-20,69,'b',size=14)
ax2.annotate('DJF CMIP6 '+var+' ensemble mean',(-13.8,69),
                                             bbox={'facecolor':'w'},zorder=20)

f = pl.gcf()
colax = f.add_axes([0.93,0.685,0.0172,0.27])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)
#cb.set_label('K',fontsize=12,labelpad=1)
cb.set_ticks(meanlevs.astype(int))
cb.set_ticklabels(meanlevs.astype(int))
cb.ax.tick_params(labelsize=9)
cb.set_label('K',fontsize=12,labelpad=1)

###############################################################################

ax3 = pl.subplot(323,projection=proj,extent=ext)
ax3.coastlines(linewidth=0.5,resolution='50m')
ax3.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax3.contourf(lon_cyc,lat.data,BIAS_e5_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=biaslevs)
ax3.contourf(lon_cyc,lat.data,CHECK_e5_cyc[0],colors='none',hatches='....')

GridLines(ax3,False,False,True,False)
ax3.text(-24,69,'c',size=14)
ax3.annotate('JJA CMIP6 '+var+' minus ERA5 '+var3,(-13.8,69),
                                             bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax4 = pl.subplot(324,projection=proj,extent=ext)
ax4.coastlines(linewidth=0.5,resolution='50m')
ax4.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax4.contourf(lon_cyc,lat.data,BIAS_e5_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=biaslevs)
ax4.contourf(lon_cyc,lat.data,CHECK_e5_cyc[1],colors='none',hatches='....')

GridLines(ax4,False,False,False,True)
ax4.text(-20,69,'d',size=14)
ax4.annotate('DJF CMIP6 '+var+' minus ERA5 '+var3,(-13.8,69),
                                             bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax5 = pl.subplot(325,projection=proj,extent=ext)
ax5.coastlines(linewidth=0.5,resolution='50m')
ax5.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax5.contourf(lon_cyc,lat.data,BIAS_eo_cyc[0],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=biaslevs)
ax5.contourf(lon_cyc,lat.data,CHECK_eo_cyc[0],colors='none',hatches='....')

GridLines(ax5,False,True,True,False)
ax5.text(-24,69,'e',size=14)
ax5.annotate('JJA CMIP6 '+var+' minus E-OBS '+var2,(-13.8,69),
                                             bbox={'facecolor':'w'},zorder=20)

###############################################################################

ax6 = pl.subplot(326,projection=proj,extent=ext)
ax6.coastlines(linewidth=0.5,resolution='50m')
ax6.add_feature(ocean_50m,alpha=1,zorder=5)
cs = ax6.contourf(lon_cyc,lat.data,BIAS_eo_cyc[1],cmap=cmap,extend='both',
             alpha=alpha,transform=ccrs.PlateCarree(),levels=biaslevs)
ax6.contourf(lon_cyc,lat.data,CHECK_eo_cyc[1],colors='none',hatches='....')

GridLines(ax6,False,True,False,True)
ax6.text(-20,69,'f',size=14)
ax6.annotate('DJF CMIP6 '+var+' minus E-OBS '+var2,(-13.8,69),
                                             bbox={'facecolor':'w'},zorder=20)

###############################################################################

pl.tight_layout()
pl.subplots_adjust(top=0.96,bottom=0.04,hspace=0.18,left=0.04,right=0.91)

f = pl.gcf()
colax = f.add_axes([0.93,0.035,0.02,0.6])
cb = pl.colorbar(cs,orientation='vertical',cax=colax)
#cb.set_label('K',fontsize=12,labelpad=1)
cb.set_ticks(biaslevs.astype(int))
cb.set_ticklabels(biaslevs.astype(int))
cb.ax.tick_params(labelsize=12)
cb.set_label('K',fontsize=12,labelpad=1)

#pl.savefig(indecis+'figures/'+var+'_cmip6_era5_eobs_mn_bias_jja_djf.png',dpi=400)


#ax1 = pl.subplot(221,projection=proj,extent=ext)
#ax1.coastlines(linewidth=0.5,resolution='50m')
#ax1.add_feature(ocean_50m,alpha=1,zorder=5)
#cs = ax1.contourf(lon_cyc,lat.data,CMIP_ssn_cyc[0],cmap='Greens',extend='max',
#             alpha=alpha,transform=ccrs.PlateCarree(),levels=meanlevs)
#
#GridLines(ax1,False,False,True,False)
#ax1.text(-25,69,'a',size=14)
#pl.title('JJA CMIP6 '+var+' ensemble mean',size=10)
##ax1.annotate('JJA CMIP6 '+var+' ensemble mean',(-13.8,69),bbox={'facecolor':'w'},zorder=20)
#
################################################################################
#
#ax2 = pl.subplot(222,projection=proj,extent=ext)
#ax2.coastlines(linewidth=0.5,resolution='50m')
#ax2.add_feature(ocean_50m,alpha=1,zorder=5)
#ax2.contourf(lon_cyc,lat.data,CMIP_ssn_cyc[1],cmap='Greens',extend='max',
#             alpha=alpha,transform=ccrs.PlateCarree(),levels=meanlevs)
#
#GridLines(ax2,False,False,False,False)
#ax2.text(-19,69,'b',size=14)
#pl.title('DJF CMIP6 '+var+' ensemble mean',size=10)
##ax2.annotate('DJF CMIP6 '+var+' ensemble mean',(-13.8,69),bbox={'facecolor':'w'},zorder=20)
#
#f = pl.gcf()
#colax = f.add_axes([0.905,0.54,0.015,0.41])
#cb = pl.colorbar(cs,orientation='vertical',cax=colax)
#cb.set_ticks(meanlevs[::2])
#cb.set_ticklabels(meanlevs[::2])
#cb.ax.tick_params(labelsize=10)
#cb.set_label('m$^3$ m$^{-1}$',fontsize=11,labelpad=5)
#
#
################################################################################
#
#ax3 = pl.subplot(223,projection=proj,extent=ext)
#ax3.coastlines(linewidth=0.5,resolution='50m')
#ax3.add_feature(ocean_50m,alpha=1,zorder=5)
#ax3.contourf(lon_cyc,lat.data,BIAS_e5_cyc[0],cmap=cmap,extend='both',
#             alpha=alpha,transform=ccrs.PlateCarree(),levels=biaslevs)
#ax3.contourf(lon_cyc,lat.data,CHECK_e5_cyc[0],colors='none',hatches='....')
#
#GridLines(ax3,False,True,True,False)
#ax3.text(-25,69,'c',size=14)
#pl.title('JJA CMIP6 '+var3+' minus ERA5 '+var3,size=10)
##ax3.annotate('JJA CMIP6 '+var3+' minus ERA5 '+var3,(-13.8,69),
##                                             bbox={'facecolor':'w'},zorder=20)
#
################################################################################
#
#ax4 = pl.subplot(224,projection=proj,extent=ext)
#ax4.coastlines(linewidth=0.5,resolution='50m')
#ax4.add_feature(ocean_50m,alpha=1,zorder=5)
#cs = ax4.contourf(lon_cyc,lat.data,BIAS_e5_cyc[1],cmap=cmap,extend='both',
#             alpha=alpha,transform=ccrs.PlateCarree(),levels=biaslevs)
#ax4.contourf(lon_cyc,lat.data,CHECK_e5_cyc[1],colors='none',hatches='....')
#
#GridLines(ax4,False,True,False,False)
#ax4.text(-19,69,'d',size=14)
#pl.title('DJF CMIP6 '+var3+' minus ERA5 '+var3,size=10)
##ax4.annotate('DJF CMIP6 '+var3+' minus ERA5 '+var3,(-13.8,69),
##                                             bbox={'facecolor':'w'},zorder=20)
#
#f = pl.gcf()
#colax = f.add_axes([0.905,0.055,0.015,0.41])
#cb = pl.colorbar(cs,orientation='vertical',cax=colax)
#cb.set_ticks(biaslevs[::2])
#cb.set_ticklabels(biaslevs[::2])
#cb.ax.tick_params(labelsize=10)
#cb.set_label('m$^3$ m$^{-1}$',fontsize=11,labelpad=-1)
#
#pl.subplots_adjust(top=0.97,left=0.08,right=0.90,hspace=0.07,wspace=0.1,bottom=0.04)
#
#pl.savefig(indecis+'figures/'+var+'_cmip6_era5_mn_bias_jja_djf.png',dpi=400)