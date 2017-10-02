import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncfile
from pylab import *
import os
import glob
import cartopy.crs as ccrs
import matplotlib as mpl
from scipy import ndimage

#Choose variable to look at:
var = 'prp'
#Choose annual cycle
diag = 'ts'
#Choose region
spatial = 'GHoA'
#Choose path where data stored
path = '/ouce-home/staff/sedm4922/Validation/'
#Choose factor to convert to mm/month

fac = 3600*24*30 # for the model
fac2 = 24*30.5 # for the obs
#mon = ['01','02','03','04','05','06','07','08','09','10','11','12']
mon = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
dimo = [31,28,31,30,31,30,31,31,30,31,30,31]

precip = np.zeros((12,25,138,200))
precip_ens = np.zeros((12,25,100,138,200))

for moncount in np.arange(len(mon)):

        if var == 'prp':
            nc=ncfile(path + mon[moncount]+'_precip_Africa_1986-2010_Dec-Nov.nc')
            temp = nc.variables['precip'][:]*fac
            temp[temp > np.std(temp)*10] = np.nan
            temp[temp < 0] = np.nan
            temp[temp == 0] = np.nan 
            precip[moncount,:,:,:] = np.nanmean(temp,axis=1)
            precip_ens[moncount,:,:,:,:] = temp
            precip_ens[:,:,40,:,:] = np.nan           
            precip_ens[:,:,11,:,:] = np.nan 

lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]


if spatial == 'Whole':

    lat1 = np.arange(len(lat))
    lon1 = np.arange(len(lon))

if spatial == 'GHoA':

    lat1 = np.array(np.where((lat < 15.63) & (lat > -4.05)))[0,:]
    lon1 = np.array(np.where((lon < 51.25) & (lon > 30)))[0,:]

if spatial == 'ME':

    lat1 = np.array(np.where((lat < 45.) & (lat > 10)))[0,:]
    lon1 = np.array(np.where((lon < 60) & (lon > 30)))[0,:] 

if diag == 'ts':
	print precip_ens.shape
	temp_precip = np.nanmean(precip_ens[:,:,:,:,lon1],axis=4)
    	precip_mod3 = np.nanmean(temp_precip[:,:,:,lat1],axis=3)
	print precip_mod3.shape
	precip_ts = np.nanmean(precip_mod3[:,:,0:65],axis=0)
	print precip_ts.shape

if var == 'prp':
            s0=np.load(path + 'prp_chrps_mon_ts.npy') 
            s1=np.load(path + 'prp_trmm_ts.npy') 
            s2=np.load(path + 'prp_cru_ts.npy')
            print s2.shape

            temp_precip = s0[:,:,lon1]
            s02 = temp_precip[:,lat1,:]
        
            temp_precip = s1[:,:,lon1]
            s12 = temp_precip[:,lat1,:]
        
            temp_precip = s2[:,:,lon1]
            s22 = temp_precip[:,lat1,:]
            
            s22[s22 == 0] = np.nan
            s12[s12 == 0] = np.nan
            s02[s02 == 0] = np.nan
	    
#Sets the levels for the plotting
#if var == 'prp':
    #levels = np.arange(10)*20
    #levels2 = -50 + np.arange(11)*10

if diag == 'ts':
 
    
    for encount in np.arange(65):
        plt.plot(np.arange(25)+1987, precip_ts[:,encount],color = 'lightblue')  
        mo, = plt.plot(np.arange(25)+1987,np.nanmean(precip_ts,axis=1),label='HadAM3P',color='blue')

        cr, = plt.plot(np.arange(25)+1987,np.average(np.nanmean(s22[2:27,:,:],axis=2),axis=1),label='CRU',color='k',linewidth=3)    
        trm, = plt.plot(np.arange(17)+1998,np.average(np.nanmean(s12[:,:,:],axis=2),axis=1),label='TRMM',color='k',linewidth=3,linestyle='-.')    
        chirps, = plt.plot(np.arange(30)+1987,np.average(np.nanmean(s02[6:36,:,:],axis=2),axis=1),label='CHIRPS',color='k',linewidth=3,linestyle='--')    
 
               
        plt.ylabel('Precip. (mm/month)')
        plt.title('Annual mean precip., Region: GHoA Africa')   
        #plt.ylim(0,100)
  	plt.legend(handles=[mo,cr,trm,chirps],loc=2,fontsize=10,frameon=False) 
    	plt.xlim(1985,2016)
    #locs, labels = plt.xticks() 
    #plt.xticks(np.arange(),mon)
plt.show()
savefig(path+'Timeseries.pdf')
