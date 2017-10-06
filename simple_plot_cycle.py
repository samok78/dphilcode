import numpy as np
from netCDF4 import Dataset as ncfile
import matplotlib.pyplot as plt
from sys import exit
import os

cwd = os.getcwd()
subpath = r'NEWINPUT/WAH_MODEL/'
fname = r'NEWINPUT/WAH_MODEL/item5216_monthly_mean_h000_1985-12_1986-12.nc'
path = os.path.join(cwd,fname)
dir_path = os.path.join(cwd,subpath)
print dir_path
print path

onlyfiles = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
print onlyfiles
files_split = [[f,f[0:-3].split('_')] for f in onlyfiles]
for el in files_split:
  fname = el.pop(0)
  el[0].append(fname)
files_split = [item for sublist in files_split for item in sublist]
#files_split = [item for sublist in files_split for item in sublist]
print files_split
keys = ['variable_name','time_period','stat','ensemble_member','start_date','end_date','filename']
files_dict = [dict(zip(keys,f)) for f in files_split]
print files_dict
selection = [el['filename'] for el in files_dict if el['end_date'][0:4]=='1990']
print selection
#convert filenames to a dict
#write a loop that opens each file and adds to an empty np array



"""
dataset = ncfile(path)
print dataset.file_format
print dataset.dimensions.keys()
print dataset.variables.keys()
lats = dataset.variables['global_latitude0'][:]
lons = dataset.variables['global_longitude0'][:]
time = dataset.variables['global_longitude0'][:]
pr = dataset.variables['item5216_monthly_mean'][:]

print lats
print lats.shape
print lons
print lons.shape
print time
print time.shape
print pr 
print pr.shape
"""

"""
path = '/ouce-home/staff/sedm4922/Validation/'

# factors to convert units to mm/month
fac = 86400 # for the model
fac2 = 30.5 # for the obs
spatial = 'GHoA'
# Define months array
mon = ['01','02','03','04','05','06','07','08','09','10','11','12']

precip = np.zeros((12,25,146,209))
precip_ens = np.zeros((12,25,100,146,209))
       
for moncount in np.arange(len(mon)):
            nc=ncfile(path + mon[moncount]+'_precip_Asia_1986-2010_Dec-Nov.nc')
            temp = nc.variables['pr'][:]*fac
	    #print temp.shape
	    #exit()
            temp[temp > np.std(temp)*10] = np.nan
            temp[temp < 0] = np.nan # if temperature in kelvin, this makes sense
            temp[temp == 0] = np.nan
            precip[moncount,:,:,:] = np.nanmean(temp,axis=1)
            precip_ens[moncount,:,:,:,:] = temp 
            precip_ens[:,:,40,:,:] = np.nan           
            precip_ens[:,:,11,:,:] = np.nan
print precip_ens.shape
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

chirps = ncfile(path + 'chirps_clim_mn_mmd_jan_dec.nc')	    
chirps2 = chirps.variables['precip'][:]
lat2 = chirps.variables['latitude'][:]
lon2 = chirps.variables['longitude'][:]
#exit()
print chirps2.shape
#exit()
#print lat.shape
#print lon.shape
#exit()
if spatial == 'GHoA':
    lat1 = np.array(np.where((lat < 15.63) & (lat > -4.05)))[0,:]
    lon1 = np.array(np.where((lon < 51.25) & (lon > 30)))[0,:]
    lat2 = np.array(np.where((lat2 < 15.63) & (lat2 > -4.05)))[0,:]
    lon2 = np.array(np.where((lon2 < 51.25) & (lon2 > 30)))[0,:]


    #lat1 = np.array(np.where((lat < 13) & (lat > 8)))[0,:]
    #lon1 = np.array(np.where((lon < 43) & (lon > 38)))[0,:]
    #lat2 = np.array(np.where((lat2 < 13) & (lat2 > 8)))[0,:]
    #lon2 = np.array(np.where((lon2 < 43) & (lon2 > 38)))[0,:]

#lat1 = np.arange(len(lat))
#lon1 = np.arange(len(lon))
#lat2 = np.arange(len(lat2))
#lon2 = np.arange(len(lon2))
print lat1.shape
print lon1.shape
print lat2
print lon2
print lat1
print lon1
temp_precip = np.nanmean(precip_ens[:,:,:,:,lon1],axis=4)
precip_mod3 = np.nanmean(temp_precip[:,:,:,lat1],axis=3)
precip_cycle = np.nanmean(precip_mod3[:,:,0:65],axis=1)

temp_obs = np.nanmean(chirps2[:,:,lon2],axis=2)
chrps = np.nanmean(temp_obs[:,lat2],axis=1)

    
for encount in np.arange(65):
        plt.plot(np.arange(12), precip_cycle[:,encount], color = 'lightblue')
    	mo, = plt.plot(np.arange(12),np.nanmean(precip_cycle,axis=1),label='HadAM3P Clim.',color='blue') 
        tams, = plt.plot(np.arange(12),chrps,label='CHIRPS Clim.',color='k',linewidth=3,linestyle='--')   


plt.ylabel('Precip. (mm/day)')
plt.title('Precipitation Annual Cycle: GHoA - ASIA Region') 
plt.legend(handles=[mo,tams],loc=2,fontsize=10,frameon=False)        
        
plt.xlim(0,11)
locs, labels = plt.xticks() 
plt.xticks(np.arange(12),mon)

plt.show()
#savefig(path+'Africa_Maps.pdf')
"""

