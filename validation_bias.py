import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkids.basemap import Basemap
from netCDF4 import Dataset as ncfile
from pylab import *
import os
import glob
import cartopy.crs as ccrs
import matplotlib as mpl
from scipy import ndimage
from sys import exit

#Choose variable to look at:
var = 'prp'
#Choose path where data stored
path = '/ouce-home/staff/sedm4922/Validation/'
#Choose factor to convert to mm/month

fac = 3600*24*30 # for the model
fac2 = 24*30.5 # for the obs

#Choose mask
nc=ncfile(path+'ethiopiamask_12.nc')
eth = nc.variables['mask'][:,:,:]


mon = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
dimo = [31,28,31,30,31,30,31,31,30,31,30,31]

precip = np.zeros((12,25,138,200))
precip_ens = np.zeros((12,25,100,138,200))

# compare with: 1 = None, 2 = TRMM, 3 = CRU, 4 = CHIRPS
compare = 'None'

precip = np.zeros((12,25,138,200))
precip_ens = np.zeros((12,25,100,138,200))

for moncount in np.arange(len(mon)): 
	if var == 'prp':
		#if not os.path.exists(path+'precip_mod.npy'):
		
            		nc=ncfile(path + mon[moncount]+'_precip_Africa_1986-2010_Dec-Nov.nc')
            		temp = nc.variables['precip'][0,0,:,:]*fac
            		#temp[temp > np.std(temp)*10] = np.nan
            		#temp[temp < 0] = np.nan # if temperature in kelvin, this makes sense
            		#temp[temp == 0] = np.nan        
            		#precip[moncount,:,:,:] = np.nanmean(temp,axis=1)
	    		precip[moncount,0,0,0,] = moncount +1
			lat = nc.variables['latitude'][:]
        		lon = nc.variables['longitude'][:]
			print temp.shape
			print temp[:,0,0,0]
			print precip[:,0,0,0,]
			precip_mod = np.nanmean(precip,axis=1)
			np.save(path+'prp_mod.npy',precip_mod)
		
		else: 
			precip_mod = np.load(path+'precip_mod.npy')
			nc=ncfile(path + mon[moncount]+'_precip_Africa_1986-2010_Dec-Nov.nc')
			lat = nc.variables['latitude'][:]
        		lon = nc.variables['longitude'][:]
			precip_mod=np.ma.masked_where(eth.mask,precip_mod)

		t0=np.load(path + 'prp_chrps_mon_clim.npy')
    		t1=np.load(path + 'prp_trmm_mon_clim.npy')
    		t2=np.load(path + 'prp_cru_mon_clim.npy')
		
		t0=np.ma.masked_where(eth.mask,t0)
		t1=np.ma.masked_where(eth.mask,t1)
		t2=np.ma.masked_where(eth.mask,t2)

# Define plotting levels and color map

    		levels = np.arange(10)*20
    		levels2 = -50 + np.arange(11)*10

#Plotting


fig = plt.figure(figsize =[8,10])
    
for moncount in np.arange(12):
	sp = plt.subplot(4,3,moncount+1,projections=ccrs.PlateCaree())
        
	if compare == 'None':
		im = plt.contourf(lon,lat,precip_mod[moncount,:,:], extend='both',levels=levels)
	    
	if compare == 'TRMM':
		im = plt.contourf(lon,lat,precip_mod2[moncount,:,:] - t12[moncount,:,:],extend='both',levels=levels2)
	    
	if compare == 'CRU':
		field =  precip_mod2[moncount,:,:] - t02t[moncount,:,:]
		im = plt.contourf(lon,lat,field,extend='both',levels=levels2)
		   
	plt.title(mon[moncount])
	sp.coastlines()
    
if var == 'prp': suptitle('Precip bias HadRM3P vs. '+compare+' (mm/month)', fontsize=15, y=0.99,x=0.5)
 
cbaxes = fig.add_axes([0.05, 0.03, 0.9, 0.02]) 
plt.colorbar(im,cax=cbaxes,orientation='horizontal')
plt.tight_layout()
plt.show()
savefig(path+'Africa_Maps.pdf')
