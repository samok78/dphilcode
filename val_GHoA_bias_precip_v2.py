# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:35:08 2015

@author: Dann
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 09:18:24 2015

@author: Dann
"""

# looking at biases in the temperature and precip for the GHoA long rains region.
# all data must be interpolated to the same grid (model grid in this case)
 
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncfile
from pylab import *
import os
import glob
import cartopy.crs as ccrs
import matplotlib as mpl
from scipy import ndimage

test_mode = 1 # if set to one, run through test diagnostics. Set to 0 for most

# Define what you want to look at here in terms of temporal and spaatial scales
# ==============================================================================
# Plot either tas or prp
var = 'prp'

# Either a monthly mean map, or a monthly mean time series, or annual cycle; map or ts or cycle
diag = 'cycle'

# temporal: look at 1 = monthly averages
temporal = 1

# spatial: look at 'Whole', 'GHoA', 'User', 'ME', 'Kericho', 'North Africa', 'West Sahel' where User must supple coords
spatial = 'GHoA'

# if user is set, define region here

lat_upper = 3.8
lat_lower = -1.8
lon_upper = 38.5
lon_lower = 33.5

# compare with: 1 = None, 2 = TRMM, 3 = CRU, 4 = None_TRMM, 5 = None_CRU
compare = 'TRMM'

# ==============================================================================

# Path where data is stored
path = '/ouce-home/staff/sedm4922/Validation/'

# factors to convert units to mm/month
fac = 3600*24*30 # for the model
fac2 = 24*30.5 # for the obs

# load landsea mask
nc=ncfile(path+'lsm_africa.nc')
lsm = nc.variables['lsm'][0,0,:,:]*-1.
lsm[lsm == 0] = np.nan

# Define months array
mon = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
dimo = [31,28,31,30,31,30,31,31,30,31,30,31]

if test_mode == 1:
    precip = np.zeros((12,25,138,200))
    precip_ens = np.zeros((12,25,100,138,200))
       
    for moncount in np.arange(len(mon)):
        
        if var == 'tas':
            nc=ncfile(path + mon[moncount]+'_tas_Africa_1986-2010_Dec-Nov.nc')
            temp = nc.variables['tas'][:]
            temp[temp > 350] = np.nan
            temp[temp < 200] = np.nan # if temperature in kelvin, this makes sense
            temp[temp == 0] = np.nan
            precip[moncount,:,:,:] = np.nanmean(temp,axis=1)
            precip_ens[moncount,:,:,:,:] = temp - 273.15
            precip_ens[7,:,44,:,:] = np.nan
            
        if var == 'prp':
            nc=ncfile(path + mon[moncount]+'_precip_Africa_1986-2010_Dec-Nov.nc')
            temp = nc.variables['precip'][:]*fac
            temp[temp > np.std(temp)*10] = np.nan
            temp[temp < 0] = np.nan # if temperature in kelvin, this makes sense
            temp[temp == 0] = np.nan
            precip[moncount,:,:,:] = np.nanmean(temp,axis=1)
            precip_ens[moncount,:,:,:,:] = temp #- 273.15
            precip_ens[:,:,40,:,:] = np.nan           
            precip_ens[:,:,11,:,:] = np.nan
        
    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]

    if var == 'tas':
        precip_mod = np.nanmean(precip,axis=1)-273.15   
        np.save(path+'tas_mod.npy',precip_mod)

    if var == 'prp':
        precip_mod = np.nanmean(precip,axis=1)   
        np.save(path+'prp_mod.npy',precip_mod)

# faster to do it this way, but less control
if test_mode == 0:
    lat = np.load(path+'lat_afr.npy')
    lon = np.load(path+'lon_afr.npy')


# define region to look over

if spatial == 'GHoA':

    lat1 = np.array(np.where((lat < 15.63) & (lat > -4.05)))[0,:]
    lon1 = np.array(np.where((lon < 51.25) & (lon > 30)))[0,:]

if spatial == 'Whole':

    #lat1 = np.array(np.where((lat < 40.) & (lat > -60)))[0,:]
    #lon1 = np.array(np.where((lon < 60) & (lon > -20)))[0,:]

    lat1 = np.arange(len(lat))
    lon1 = np.arange(len(lon))

if spatial == 'User':

    lat1 = np.array(np.where((lat < lat_upper) & (lat > lat_lower)))[0,:]
    lon1 = np.array(np.where((lon < lon_upper) & (lon > lon_lower)))[0,:]
    
if spatial == 'ME':

    lat1 = np.array(np.where((lat < 45.) & (lat > 10)))[0,:]
    lon1 = np.array(np.where((lon < 60) & (lon > 30)))[0,:] 
    
if spatial == 'Kericho':

    lat1 = np.array(np.where((lat < 2.8) & (lat > -1.8)))[0,:]
    lon1 = np.array(np.where((lon < 37.5) & (lon > 33.5)))[0,:] 

if spatial == 'North Africa':

    lat1 = np.array(np.where((lat < 36) & (lat > 20)))[0,:]
    lon1 = np.array(np.where((lon < 40) & (lon > -20)))[0,:] 
    
if spatial == 'West Sahel':

    lat1 = np.array(np.where((lat < 20) & (lat > 10)))[0,:]
    lon1 = np.array(np.where((lon < 0) & (lon > -20)))[0,:]

# Read in model data and 3 observational estimates
# Note that TAMSAT looks crap, so will not use it for now.

# set up dummy arrays
precip_mod = np.zeros((12,200,200))
t0 = np.zeros((12,200,200))
t1 = np.zeros((12,200,200))
t2 = np.zeros((12,200,200))

if var == 'prp':
    precip_mod = np.load(path+'precip_mod.npy')
    
    t0=np.load(path + 'prp_tam_mon_clim.npy')
    t1=np.load(path + 'prp_trmm_mon_clim.npy')
    t2=np.load(path + 'prp_cru_mon_clim.npy')

if var == 'tas':
    precip_mod = np.load(path+'tas_mod.npy')
    
    t2=np.load(path + 'tas_cru_mon_clim.npy')
    
# select only the region of interest



temp_precip = precip_mod[:,:,lon1]
precip_mod2 = temp_precip[:,lat1,:]

if test_mode == 1:
    temp_precip = np.nanmean(precip_ens[:,:,:,:,lon1],axis=4)
    precip_mod3 = np.nanmean(temp_precip[:,:,:,lat1],axis=3)
    
    if diag == 'cycle':
        precip_cycle = np.nanmean(precip_mod3[:,:,0:65],axis=1)
        
    if diag == 'ts':
        precip_ts = np.nanmean(precip_mod3[:,:,0:65],axis=0)    
    
        if var == 'tas':    
            s2=np.load(path + 'tas_cru_mon_ts.npy')
    
            temp_precip = s2[:,:,lon1]
            s22 = temp_precip[:,lat1,:]
            s22[s22 == 0] = np.nan
            
        if var == 'prp':
            s0=np.load(path + 'prp_tam_ts.npy')
            s1=np.load(path + 'prp_trmm_ts.npy')
            s2=np.load(path + 'prp_cru_ts.npy')
    
            temp_precip = s0[:,:,lon1]
            s02 = temp_precip[:,lat1,:]
        
            temp_precip = s1[:,:,lon1]
            s12 = temp_precip[:,lat1,:]
        
            temp_precip = s2[:,:,lon1]
            s22 = temp_precip[:,lat1,:]
            
            s22[s22 == 0] = np.nan
            s12[s12 == 0] = np.nan
            s02[s02 == 0] = np.nan

temp_precip = t0[:,:,lon1]
t02t = temp_precip[:,lat1,:]

temp_precip = t1[:,:,lon1]
t12t = temp_precip[:,lat1,:]

temp_precip = t2[:,:,lon1]
t22t = temp_precip[:,lat1,:]

# NOT SURE WHAT THESE DO, BUT IF PRECIP =0 IT SHOULD BE OK, SO I COMMENTED OUT ONE OF THEM FOR BIASES
t02t[t02t == 0] = np.nan
t12t[t12t == 0] = np.nan
#t22t[t22t == 0] = np.nan  

temp_precip = lsm[:,lon1]
lsm = temp_precip[lat1,:]

lon = lon[lon1]
lat = lat[lat1]

# Define plotting levels and color map
if var == 'prp':
    levels = np.arange(10)*20
    levels2 = -50 + np.arange(11)*10
    
if var == 'tas':
    levels = 0 + np.arange(10)*4
    levels2 = -5 + np.arange(11)*1   
    
my_cmap = matplotlib.cm.get_cmap('bwr')
#if var == 'tas': my_cmap = matplotlib.cm.get_cmap('Reds')
# plotting
# ==============================================================================

if diag == 'cycle':
    
    t22t[t22t > 1000] = np.nan
    
    for encount in np.arange(65):
        plt.plot(np.arange(12), precip_cycle[:,encount], color = 'lightblue')
    
    mo, = plt.plot(np.arange(12),np.nanmean(precip_cycle,axis=1),label='HadAM3P Clim.',color='blue')
    
    if var == 'tas':
        cr, = plt.plot(np.arange(12),np.average(np.nanmean(t22t,axis=2),axis=1),label='CRU Clim.',color='k',linewidth=3)
        plt.ylabel('Temperature (Degrees C)')
        plt.title('Temperature annual cycle, Region: ' + spatial)
        
    if var == 'prp':
        cr, = plt.plot(np.arange(12),np.average(np.nanmean(t22t,axis=2),axis=1),label='CRU Clim.',color='k',linewidth=3)    
        trm, = plt.plot(np.arange(12),np.average(np.nanmean(t12t,axis=2),axis=1),label='TRMM Clim.',color='k',linewidth=3,linestyle='-.')    
        tams, = plt.plot(np.arange(12),np.average(np.nanmean(t02t,axis=2),axis=1),label='TAMSAT Clim.',color='k',linewidth=3,linestyle='--')    


        plt.ylabel('Precip. (mm/month)')
        plt.title('Precip. annual cycle, Region: ' + spatial) 
        
    if var == 'tas': plt.legend(handles=[mo,cr],loc=2,fontsize=10,frameon=False)
    if var == 'prp': plt.legend(handles=[mo,cr,trm,tams],loc=2,fontsize=10,frameon=False)        
        
    plt.xlim(0,11)
    locs, labels = plt.xticks() 
    plt.xticks(np.arange(12),mon)

if diag == 'ts':
    
    t22t[t22t > 1000] = np.nan
    
    for encount in np.arange(65):
        plt.plot(np.arange(25)+1987, precip_ts[:,encount],color = 'lightblue')
    
    if var == 'tas':

        mo, = plt.plot(np.arange(25)+1987,np.nanmean(precip_ts,axis=1),label='HadAM3P',color='blue')
        cr, = plt.plot(np.arange(25)+1987,np.average(np.nanmean(s22[86:111,:,:],axis=2),axis=1),label='CRU',color='k',linewidth=3)    

        plt.ylabel('Temperature (Degrees C)')
        plt.title('Annual mean temperature, Region: ' + spatial)
        
    if var == 'prp':
        
        mo, = plt.plot(np.arange(25)+1987,np.nanmean(precip_ts,axis=1),label='HadAM3P',color='blue')
        cr, = plt.plot(np.arange(25)+1987,np.average(np.nanmean(s22[2:27,:,:],axis=2),axis=1),label='CRU',color='k',linewidth=3)    
        trm, = plt.plot(np.arange(17)+1998,np.average(np.nanmean(s12[:,:,:],axis=2),axis=1),label='TRMM',color='k',linewidth=3,linestyle='-.')    
        #tams, = plt.plot(np.arange(12)+1987,np.average(np.nanmean(s02[:,:,:],axis=2),axis=1),label='TAMSAT',color='k',linewidth=3,linestyle='--')    
 
               
        plt.ylabel('Precip. (mm/month)')
        plt.title('Annual mean precip., Region: ' + spatial)   
        #plt.ylim(0,100)
    #plt.legend(handles=[mo,cr],loc=2,fontsize=10,frameon=False)
    if var == 'tas': plt.legend(handles=[mo,cr],loc=2,fontsize=10,frameon=False)
    if var == 'prp': plt.legend(handles=[mo,cr,trm],loc=2,fontsize=10,frameon=False) 
    plt.xlim(1985,2012)
    #locs, labels = plt.xticks() 
    #plt.xticks(np.arange(),mon)

# calculate bias smoothed and output for bias correction
# NOTE COMMENTED OUT SOMETHING ABOVE TO MAKE THIS WORK

bias_t = np.zeros((t22t.shape))
bias_t[:,:,:] = np.nan

for tcount in np.arange(12):
    bias2 = precip_mod2[tcount,:,:] - t22t[tcount,:,:]*lsm
    bias_t[tcount,:,:] = ndimage.filters.median_filter(bias2,size=4,mode='nearest')

    for lcount in np.arange(len(t22t[0,:,0])):
        for llcount in np.arange(len(t22t[0,0,:])):
            if np.isnan(bias_t[tcount,lcount,llcount]) == True and np.isnan(bias2[lcount,llcount]) == False:
                bias_t[tcount,lcount,llcount] = bias2[lcount,llcount]

if var == 'tas': np.save(path + 'monthly_bias_smoothed_t.npy',bias_t)
if var == 'prp': np.save(path + 'monthly_bias_smoothed_prp.npy',bias_t)   

if diag == 'map':

    fig = plt.figure(figsize =[8,10])
    
    for moncount in np.arange(12):
        sp = plt.subplot(4,3,moncount+1,projection=ccrs.PlateCarree())
        
        if compare == 'None':
            im = plt.contourf(lon,lat,precip_mod2[moncount,:,:]*lsm,extend='both',levels=levels)
    
        if compare == 'TRMM':
            im = plt.contourf(lon,lat,precip_mod2[moncount,:,:] - t12t[moncount,:,:]*lsm,extend='both',levels=levels2,cmap= my_cmap)
    
        if compare == 'CRU':
            #im = plt.contourf(lon,lat,precip_mod2[moncount,:,:] - t22t[moncount,:,:]*lsm,extend='both',levels=levels2,cmap= my_cmap)
            im = plt.contourf(lon,lat,bias_t[moncount,:,:],extend='both',levels=levels2,cmap= my_cmap)
    
        if compare == 'None_TRMM':
            im = plt.contourf(lon,lat,t12t[moncount,:,:]*lsm,extend='both',levels=levels)
    
        if compare == 'None_CRU':
            im = plt.contourf(lon,lat,t22t[moncount,:,:]*lsm,extend='both',levels=levels)
    
        plt.title(mon[moncount])
        sp.coastlines()
    
    if var == 'prp': suptitle('Precip bias HadRM3P vs. '+compare+' (mm/month)', fontsize=15, y=0.99,x=0.5)
    if var == 'tas': suptitle('Temp bias HadRM3P vs. '+compare+' (K)', fontsize=15, y=0.99,x=0.5)
    
    
    cbaxes = fig.add_axes([0.05, 0.03, 0.9, 0.02]) 
    plt.colorbar(im,cax=cbaxes,orientation='horizontal')
    plt.tight_layout()
    plt.show()
    #savefig('/Users/Dann/CPDN/Africa/validation_precip_mon_maps.pdf')

