# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 16:13:23 2019

@author: pmcraig
"""

from __future__ import division
import pylab as pl
import glob
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import timeit
import matplotlib.ticker as mticker
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
from adjustText import adjust_text
import pcraig_funcs as pc

def SixHourSteps(DATESEL,HOURSEL,inittime):
    """
    """
    ind = int((DATESEL-1)*4 + HOURSEL/6)
    frc = 0
    time = inittime.values[ind]
    time = list(time)
    time[0], time[3] = time[3], time[0]
    time[1], time[4] = time[4], time[1]
    time = ''.join(time)
    print time
    
    return ind, frc, time

def ForecastSteps(DATESEL,HOURSEL,inititme,foretime):
    """
    """
    F = pl.array([foretime.values[0],foretime.values[1]],dtype='timedelta64[h]')
    ind = int((DATESEL-1)*4 + int(HOURSEL/6))
    frc = 1
    time = inittime.values[ind]
    time = list(time)
    time[0], time[3] = time[3], time[0]
    time[1], time[4] = time[4], time[1]
    HR = str("{0:0=2d}".format(int(time[-6:-4])+F[1].astype('int')))
    time[-6] = HR[0]; time[-5] = HR[1]
    time = ''.join(time)
    print time
    
    return ind, frc, time

def BetweenSteps(DATESEL,HOURSEL,inittime):
    """
    """
    H = int(HOURSEL/3)
    ind1 = int((DATESEL-1)*4 + int(HOURSEL/6))#; ind2 = (DATESEL-1)*4 + H+1
    if H % 2 == 0:
        ind2 = ind1
        frc1 = 0; frc2 = 1
    else:
        ind2 = ind1 + 1
        frc1 = 1; frc2 = 0
    time = inittime.values[ind1]
    time = list(time)
    time[0], time[3] = time[3], time[0]
    time[1], time[4] = time[4], time[1]
    HR = str("{0:0=2d}".format(HOURSEL))
    time[-6] = HR[0]; time[-5] = HR[1]
    time = ''.join(time)
    print time
    
    return ind1, ind2, frc1, frc2, time

def MakeAverage(data_in,ind1,ind2,frc1,frc2,HOURSEL):
    """
    """
    N = HOURSEL-int(HOURSEL/3)*3 # number of hours to muliply dp/dt below
    M = data_in[ind1,frc1]
    M2 = data_in[ind2,frc2]
    dp_dt = (M2-M)/(3*60*60)
    MI = M + dp_dt*(N*60*60)
    
    return MI

def GetDWRobs(ncasdir,HOURSEL,time,year,statinds):
    """
    """
    dwrfiles = glob.glob(ncasdir+'DWRcsv/'+year+'/*')
    dwrfiles = pl.sort(dwrfiles)
    
    datestrings = pl.zeros([dwrfiles.shape[0],10],dtype='object')
    for i in range(dwrfiles.shape[0]):
        datestrings[i] = pl.asarray(list(dwrfiles[i][-14:-4]))
    
    T = pl.asarray(list(time))
    
    dwrind = pl.where((datestrings[:,5]==T[3]) & (datestrings[:,6]==T[4])
                       & (datestrings[:,-2]==T[0]) & (datestrings[:,-1]==T[1]))
    
    if HOURSEL < 12:
        dwrind = dwrind[0][0]
        obsfile = pl.genfromtxt(dwrfiles[dwrind],delimiter=',')
        minus99999 = pl.where(obsfile==-99999.0)
        obsfile[minus99999[0],minus99999[1]] = pl.float32('nan')
        obsdata = obsfile[statinds[0][:-1],3]
    elif HOURSEL > 12:
        dwrind = dwrind[0][0] + 1
        obsfile = pl.genfromtxt(dwrfiles[dwrind],delimiter=',')
        minus99999 = pl.where(obsfile==-99999.0)
        obsfile[minus99999[0],minus99999[1]] = pl.float32('nan')
        obsdata = obsfile[statinds[0],1]

    return obsdata

def MSLP_contours(axx,X,Y,mslp,mslp_levs,land):
    """
    """
    axx.set_extent([-20.001,15.001,37.999,65],ccrs.PlateCarree())
    axx.coastlines(color='grey',resolution='50m',linewidth=0.5)
    cn = axx.contour(X,Y,mslp/100,colors='grey',levels=mslp_levs,
                    alpha=0.5,transform=ccrs.PlateCarree(),linewidths=0.75)
    
    axx.add_feature(land_50m,alpha=0.5)
    pl.clabel(cn,inline=True,fmt="%.0f",zorder=3,inline_spacing=5,manual=True,
              fontsize=9)
    
    return None

def SPRD_contours(axx,X,Y,sprd,sprd_levs):
    """
    """
    axx.set_extent([-20.001,15.001,37.999,65],ccrs.PlateCarree())
    axx.coastlines(color='grey',resolution='50m',linewidth=0.5)
    cs = axx.contourf(X,Y,sprd/100,norm=pl.Normalize(0,12),#colors='grey',
            levels=sprd_levs,alpha=0.4,cmap='OrRd',
            extend='max',transform=ccrs.PlateCarree())
    
    return cs

def StationMarkers(axx,coords):
    """
    """
    axx.plot(coords[:,1],coords[:,0],marker='.',color='k',linewidth=0,
        transform=ccrs.PlateCarree(),alpha=0.5,ms=5)
    
    return None

def PresLabels(presobs,coords):
    """
    Args:
        presobs (array): pressure observations
        coords (array): latitude, longitude co-ordinates of stations
    """
    p2 = pl.around(presobs,0); p2 = p2[~pl.isnan(presobs)].astype(int)
    crd_lon = coords[:,1][~pl.isnan(presobs)]
    crd_lat = coords[:,0][~pl.isnan(presobs)]
    texts = [pl.text(crd_lon[i],crd_lat[i],p2[i].astype(str),zorder=10,size=9) for i in range(len(p2))]
    adjust_text(texts,avoid_text=True,avoid_self=False)

    return None

def ErrorCalc(varobs,ensmean,spread,loc,lon,lat):
    """
    varmean (array): ensemble mean
    varsprd (array): ensemble spread
    locs (array): longitude, latitude co-ordinates of stations
    lon (array): longitude array from 20CR
    lat (array): latitude array from 20CR
    """
    errs = pl.zeros_like(varobs)
    for i in range(len(varobs)):
        #print i, varobs[i]
        if pl.isnan(varobs[i]) == True:
            errs[i] = pl.float32('nan')
        else:
            statmean = pc.BilinInterp((loc[i,1],loc[i,0]),lon,lat,ensmean)
            #print statmean.values
            statsprd = pc.BilinInterp((loc[i,1],loc[i,0]),lon,lat,spread)
            #print statsprd.values
            errs[i] = (statmean - varobs[i])/statsprd
    
    return errs

def ErrLabels(errors,coords,lon,lat):
    """
    """

    e2 = pl.around(errors,2); e2 = e2[~pl.isnan(errors)]
    crd_lon = coords[:,1][~pl.isnan(errors)]
    crd_lat = coords[:,0][~pl.isnan(errors)]
    texts = [pl.text(crd_lon[i],crd_lat[i],e2[i].astype(str),size=9) for i in range(len(e2))]
    adjust_text(texts,avoid_text=True,avoid_self=False)
    
    return None

def GridLines(ax,top,left,right,bottom):
    """
    """
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
                  zorder=3)
    gl.xlabels_top = top; gl.xlabels_bottom = bottom
    gl.ylabels_left = left; gl.ylabels_right = right
    gl.xlocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20,26,30])
    gl.ylocator = mticker.FixedLocator([30,40,50,60,70])
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'k','size':10}
    gl.ylabel_style = {'color': 'k','size':10}
    
    return None

def ISPDstations(axx,ispddir,year,month,DATESEL,HOURSEL):
    """
    """
    ispdfile = ispddir + 'ISPD47_'+year+'.txt'
    with open(ispdfile,'r') as f:
        ispddata = f.readlines()
    
    ispd_lons = []; ispd_lats = []
    for i in range(len(ispddata)):
        if ispddata[i][4:6] == month and ispddata[i][6:8] == str(DATESEL) and\
                                HOURSEL-3<=int(ispddata[i][8:10])<=HOURSEL+3\
                            and 40<=float(ispddata[i][35:40])<=65:
            if int(ispddata[i][8:10]) == HOURSEL+3 and int(ispddata[i][10:12])>0:
                pass
            else:
                ispd_lons.append(float(ispddata[i][27:34]))
                ispd_lats.append(float(ispddata[i][35:40]))
    
    ispd_locs = pl.array([pl.asarray(ispd_lons),pl.asarray(ispd_lats)]).T
    
    axx.plot(ispd_locs[:,0],ispd_locs[:,1],marker='x',color='k',linewidth=0,
        transform=ccrs.PlateCarree(),alpha=0.5,ms=6)
    
    return None

def GravityCorrection(year,presobs):
    """
    """
    if int(year) == 1902:#i == 0:
        presobs[9] = presobs[9] + 0.04 # Sumburgh Head
        Nx = [10,11,24,25,26,27,28] # north of Malin Head & Shields
        presobs[Nx] = presobs[Nx] + 0.03
        # south of Malin Head & Shields, north of Jersey & Scilly
        Sx = [12,13,14,15,16,17,18,19,22,23,29,30,31,32,33]
        presobs[Sx] = presobs[Sx] + 0.02
        presobs[[20,21]] + 0.01 # Scilly & Jersey
    elif int(year) == 1903:
        presobs[9] = presobs[9] + 0.04 # Sumburgh Head
        Nx = [10,11,24,25,26,27,28] # north of Malin Head & Shields
        presobs[Nx] = presobs[Nx] + 0.03
        # south of Malin Head & Shields, north of Jersey & Scilly
        Sx = [12,13,14,15,16,17,18,19,22,23,29,30,31,32,33,-2,-1]
        presobs[Sx] = presobs[Sx] + 0.02
        presobs[[20,21]] + 0.01 # Scilly & Jersey
    
    return presobs

def AddColourBar(cs,levels):
    """
    """
    f = pl.gcf()
    colax = f.add_axes([0.2,0.04,0.6,0.02])                   
    cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
    cb.set_ticks(levels)
    cb.set_label('hPa',fontsize=12,labelpad=0)
    cb.ax.tick_params(labelsize=10,direction='in',pad=2)
    
    return None

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
CR20dir = '/gws/nopw/j04/glosat/development/data/raw/20CRv3/451/subdaily/'
#CR20dir = ncasdir + '20CR/'
ispddir = '/gws/nopw/j04/ncas_climate_vol1/users/ehawkins/ISPD/'

years = ['1863','1872']
months = ['01','01']
dates = [20,18]
hours = [8,8]

fig, ax = pl.subplots(2,2,figsize=(11,8.5))
mslp_levs = pl.linspace(952,1060,28)
sprd_levs = pl.linspace(0,12,13)#[1,1.5,2,2.5,3,3.5,4.5,5,6,7,8]
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

pl.tight_layout()
pl.subplots_adjust(bottom=0.06,top=0.98,hspace=0.02,wspace=0.08,right=0.96)

# loop over the two events
for i in range(len(years)):
    YR = years[i]
    MN = months[i]
    MONSEL = int(months[i])
    DATESEL = dates[i]
    HOURSEL = hours[i]
    
    data = pl.genfromtxt(ncasdir+'latlon_files/latlon'+YR+'_pc.txt')
    if HOURSEL < 12:
        #stats_uki = pl.where(data[:,-4]==HOURSEL)
        allstats = pl.where((HOURSEL-3<data[:,-4]) & (data[:,-4]<HOURSEL+3))
    elif HOURSEL > 12:
        #stats_uki = pl.where(data[:,-2]==HOURSEL)
        allstats = pl.where((HOURSEL-3<data[:,-2]) & (data[:,-2]<HOURSEL+3))
    
    coords = data[allstats[0]][:,1:3]
    
    allfiles = glob.glob(CR20dir+YR+'/PRMSL*')
    enslen = len(allfiles)
     # 256 lat, 512 lon
    #extract range lon=(-45,45), lat=(20,80)
    presdata = pl.zeros([enslen,128,512])
    
    for nc in range(enslen):
        ncfile = xr.open_dataset(allfiles[nc])
        if nc == 0:
            lat = xr.DataArray(ncfile.lat[:128]).data
            lon = xr.DataArray(ncfile.lon).data
            time = xr.DataArray(ncfile.time).data
                
            timelist = [str(pd.to_datetime(j)) for j in time]
            month = pl.asarray([j[5:7] for j in timelist]).astype(int)
            day = pl.asarray([j[8:10] for j in timelist]).astype(int)
            hour = pl.asarray([j[11:13] for j in timelist]).astype(int)
        
        if HOURSEL % 3 == 0:
            if nc == 0:
    
                A = pl.where((month==MONSEL) & (day==DATESEL) & (hour==HOURSEL))
                ind = A[0][0]
        
            presdata[nc,:,:] = xr.DataArray(ncfile.PRMSL[ind,:128,:]).data
            ncfile.close()
        else:
            # need to interpolate between adjacent timesteps
            if nc == 0:
                H = hour[:8]
                NI = pc.NearestIndex(H,HOURSEL)
                if NI < HOURSEL:
                    lower = H[NI]; upper = H[NI+1]
                else:
                    lower = H[NI-1]; upper = H[NI]
                
                A = pl.where((month==MONSEL) & (day==DATESEL) & (hour==lower))
                B = pl.where((month==MONSEL) & (day==DATESEL) & (hour==upper))
                ind1 = A[0][0]; ind2 = B[0][0]
            
            # extract pres at timestep before HOURSEL
            P1 = xr.DataArray(ncfile.PRMSL[ind1,:128,:]).data
            # extract pres at timestep after HOURSEL
            P2 = xr.DataArray(ncfile.PRMSL[ind2,:128,:]).data
            ncfile.close()
            
            # interpolate
            N = HOURSEL - int(HOURSEL/3)*3
            dp_dt = (P2-P1)/(3*60*60)
            presdata[nc,:,:] = P1 + dp_dt*(N*60*60)

    pres_mn = pl.mean(presdata,axis=0)
    pres_sd = pl.std(presdata,axis=0)
    
    lontemp = pl.zeros_like(lon)
    lontemp[:int(lon.size/2)] = lon[int(lon.size/2):] - 360
    lontemp[int(lon.size/2):] = lon[:int(lon.size/2)]
    lon = lontemp.copy()
    del lontemp
    
    ptemp = pl.zeros_like(pres_mn)
    ptemp[:,:int(lon.size/2)] = pres_mn[:,int(lon.size/2):]
    ptemp[:,int(lon.size/2):] = pres_mn[:,:int(lon.size/2)]
    pres_mn = ptemp.copy()
    del ptemp
    
    stemp = pl.zeros_like(pres_sd)
    stemp[:,:int(lon.size/2)] = pres_sd[:,int(lon.size/2):]
    stemp[:,int(lon.size/2):] = pres_sd[:,:int(lon.size/2)]
    pres_sd = stemp.copy()
    del stemp

   # obsdata = GetDWRobs(ncasdir,HOURSEL,time,YR,allstats)
    # do gravity correction here
    #obsdata = GravityCorrection(YR,obsdata)
    
    ax1 = pl.subplot(2,2,i+1,projection=ccrs.PlateCarree())
    ISPDstations(ax1,ispddir,YR,MN,DATESEL,HOURSEL)
    MSLP_contours(ax1,lon,lat,pres_mn,mslp_levs,land_50m)
    
    df1 = pd.read_csv(ncasdir+'latlon_files/latlon'+YR+'_pc.txt',
                                                  header=None,delimiter=' ')
    statdata = pl.array(df1)
    df2 = pd.read_csv(ncasdir+'DWRcsv/'+YR+'/DWR_'+YR+'_'+"{0:0=2d}".format(MONSEL)+\
                    '_'+"{0:0=2d}".format(DATESEL)+'.csv',header=None,skiprows=1)
    logs = pl.array(df2)
    
    pvals = pl.where(logs[:,1]!=999.00)[0]
    coords = pl.zeros([pvals.size,2])
    presobs = logs[pvals,1].astype(float)
    
    for k in range(pvals.size):
        statind = pl.where(statdata[:,0]==logs[pvals[k],0])[0]
        coords[k,0] = statdata[statind,1]
        coords[k,1] = statdata[statind,2]
    
    if i == 0:
        #correct for gravity and change to hPa:
        shieldsind = pl.where(logs[pvals,0]=='SHIELDS')[0][0]
        N = pl.where(coords[:,0]>=coords[shieldsind,0])[0]# Shields and north (not Sumburgh Head)
        presobs[N] = presobs[N] + 0.03
        
        jerseyind = pl.where(logs[pvals,0]=='JERSEY')[0][0]
        heldind = pl.where(logs[pvals,0]=='HELDER')[0][0]
        S = pl.where((coords[:,0]>coords[jerseyind,0]) & (coords[:,0]<coords[shieldsind,0]))[0]
        S = pl.delete(S,pl.where(S==heldind)[0][0])
        presobs[S] = presobs[S] + 0.02
        presobs[jerseyind] = presobs[jerseyind] + 0.01
        
        presobs = presobs*33.8639
        
        #obs_in = pl.delete(obsdata,[0,6,obsdata.size-1])
        crd_in = pl.delete(coords,[0,6,coords.shape[0]-1],axis=0)
        StationMarkers(ax1,coords)
        PresLabels(presobs,coords)
    elif i == 1:
        #correct for gravity and change to hPa:
        sumind = pl.where(logs[pvals,0]=='SUMBURGHHEAD')[0][0]
        presobs[sumind] = presobs[sumind] + 0.04
        #shieldsind = pl.where(logs[pvals,0]=='SHIELDS')[0][0]
        norstats = ['THURSO','WICK','NAIRN','ABERDEEN','LEITH','ARDROSSAN',
                    'GREENCASTLE','SHIELDS']
        N = [pl.where(logs[pvals,0]==k)[0][0] for k in norstats]
        #N = pl.where(coords[:,0]>=coords[shieldsind,0])[0]# Shields and north (not Sumburgh Head)
        presobs[N] = presobs[N] + 0.03
        
        soustats = ['SCARBOROUGH','LIVERPOOL','HOLYHEAD','YARMOUTH','VALENTIA',
                    'ROCHESPOINT','LONDON','DOVER','PORTSMOUTH','PLYMOUTH']
        S = [pl.where(logs[pvals,0]==k)[0][0] for k in soustats]
        
        scillyind = pl.where(logs[pvals,0]=='SCILLY')[0][0]
        
        presobs[S] = presobs[S] + 0.02
        presobs[scillyind] = presobs[scillyind] + 0.01
        
        presobs = presobs*33.8639
        
        #obs_in = pl.delete(obsdata,[0,6,obsdata.size-5,obsdata.size-4,
         #                                                   obsdata.size-3])
        crd_in = pl.delete(coords,[0,6,coords.shape[0]-6,coords.shape[0]-5,
                                   coords.shape[0]-4,coords.shape[0]-1],axis=0)
        StationMarkers(ax1,coords[:-1])
        PresLabels(presobs,coords)
    
    ax2 = pl.subplot(2,2,i+3,projection=ccrs.PlateCarree())
    cs = SPRD_contours(ax2,lon,lat,pres_sd,sprd_levs)

    errs = ErrorCalc(presobs,pres_mn/100,pres_sd/100,coords,lon,lat)
    StationMarkers(ax2,coords)
    ErrLabels(errs,coords,lon,lat)
    

    if i == 0:
        GridLines(ax1,True,True,False,False)
        GridLines(ax2,False,True,False,False)
        label = str(dates[0]) + '/' + months[0] + '/' + years[0] + ' ' + \
                    '(' + "{0:0=2d}".format(hours[0]) + ':00)'
        ax1.annotate('20CRv3 ensemble mean mslp\n& DWR observations '+label,#+time,
                     (-19.5,62.5),bbox={'facecolor':'w'},size=9)
        ax2.annotate('20CRv3 ensemble spread & z-scores\n'+label,#+time,
                     (-19.5,62.5),bbox={'facecolor':'w'},size=9)
    elif i == 1:
        GridLines(ax1,True,False,True,False)
        GridLines(ax2,False,False,True,False)
        label = str(dates[1]) + '/' + months[1] + '/' + years[1] + ' ' + \
                    '(' + "{0:0=2d}".format(hours[1]) + ':00)'
        ax1.annotate('20CRv3 ensemble mean mslp\n& DWR observations '+label,#+time,
                     (-19.5,62.5),bbox={'facecolor':'w'},size=9)
        ax2.annotate('20CRv3 ensemble spread & z-scores\n'+label,#+time,
                     (-19.5,62.5),bbox={'facecolor':'w'},size=9)



AddColourBar(cs,sprd_levs)
fig.text(0.01,0.95,'(a)',size=12)
fig.text(0.49,0.95,'(b)',size=12)
fig.text(0.01,0.49,'(c)',size=12)
fig.text(0.49,0.49,'(d)',size=12)

pl.savefig(jashome+'/mslp_obs_spread_panels_paper_2022.png',dpi=400)
pl.savefig(jashome+'/mslp_obs_spread_panels_paper_2022.pdf',dpi=400)