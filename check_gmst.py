# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 13:58:54 2021

@author: pmcraig
"""

import pylab as pl
import xarray as xr
import glob

pl.close('all')

jashome = '/home/users/pmcraig/'
ncasdir = '/gws/nopw/j04/ncas_climate_vol1/users/pmcraig/'
cmipdir = ncasdir + 'CMIP6/'

meanfile = xr.open_dataset(cmipdir+'standard_grid/tas_cmip6_ensmean_gm_new.nc')
tasmean = xr.DataArray(meanfile.tas)
meantime = xr.DataArray(meanfile.time)
meanfile.close()

tasmean = pl.squeeze(tasmean)

insts = ['AS-RCEC','AWI','BCC','CCCma','CMCC','CNRM-CERFACS','CSIRO','CSIRO-ARCCSS',
         'EC-Earth-Consortium','FIO-QLNM','MIROC','MOHC','MPI-M','MRI','NCAR',
         'NCC','NOAA-GFDL']

models = ['TaiESM1','AWI-ESM-1-1-LR','BCC-ESM1','CanESM5','CMCC-ESM2',
          'CNRM-ESM2-1','ACCESS-ESM1-5','ACCESS-CM2','EC-Earth3',
          'FIO-ESM-2-0','MIROC6','UKESM1-0-LL','MPI-ESM1-2-HR','MRI-ESM2-0',
          'CESM2','NorESM2-MM','GFDL-ESM4']

variants = ['r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f2',
            'r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f2',
            'r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1']

grid = ['gn','gn','gn','gn','gn','gr','gn','gn','gr','gn','gn','gn','gn','gn',
        'gn','gn','gr1']

GMST = pl.zeros([len(insts),259])
tlist = GMST.copy()

filenames = glob.glob(cmipdir+'global_means/*')

for i in range(len(insts)):
    #path = cmipdir+insts[i]+'/'+models[i]+'/'+variants[i]+'/'
    #filename = glob.glob(cmipdir+'global_means)
    ncfile = xr.open_dataset(filenames[i])
    tasmodel = xr.DataArray(ncfile.tas)
    time = xr.DataArray(ncfile.time)
    ncfile.close()
    
    tlist[i] = time.data
    
    #if tasmodel.shape[0] == 261:
    GMST[i] = pl.squeeze(tasmodel.data[:])
    #else:
    #    GMST[i] = pl.squeeze(tasmodel.data[401:-1])

GMST_mn = pl.mean(GMST,axis=0)

pl.plot(tasmean[:],lw=4,label='mean from cdo')
pl.plot(GMST_mn[:],lw=1.5,label='mean from Python',zorder=10)
pl.legend()
pl.xlim(0,260)