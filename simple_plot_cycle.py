import numpy as np
from netCDF4 import Dataset as ncfile
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.basemap import Basemap, cm
from sys import exit
import os


def main():
    cwd = os.getcwd()
    subpath = r'NEWINPUT/WAH_MODEL/'
    fname = r'NEWINPUT/WAH_MODEL/item5216_monthly_mean_h000_1985-12_1986-12.nc'
    path = os.path.join(cwd,fname)
    dir_path = os.path.join(cwd,subpath)
    print dir_path
    print path
    
    
    #get lat/lon dims
    dataset = ncfile(path)
    lats = dataset.variables['global_latitude0'][:]
    lons = dataset.variables['global_longitude0'][:]
    lat_dim = lats.shape[0]
    lon_dim = lats.shape[1]
    
    
    #get all file names
    onlyfiles = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
    #print onlyfiles
    files_split = [[f,f[0:-3].split('_')] for f in onlyfiles]
    for el in files_split:
      fname = el.pop(0)
      el[0].append(fname)
    files_split = [item for sublist in files_split for item in sublist]
    #files_split = [item for sublist in files_split for item in sublist]
    #print files_split
    keys = ['variable_name','time_period','stat','ensemble_member','start_date','end_date','filename']
    files_dict = [dict(zip(keys,f)) for f in files_split]
    #print files_dict
    for el in files_dict:
      el['end_year'] = el['end_date'][0:4]
      el['start_year'] = el['start_date'][0:4]
    selection = [el['filename'] for el in files_dict if el['end_date'][0:4]=='1990']
    
    #get number of years
    years = list(set([el['end_year'] for el in files_dict]))
    print years
    year_dim = len(years)
    
    
    #get max number of runs
    n_runs = []
    for year in years:
      selection = [el for el in files_dict if el['end_year']==year]
      n_runs.append(len(selection))
    
    max_runs = max(n_runs)
    print max_runs
    
    #initialise the data arrays
    
    precip = np.zeros((13,year_dim,lat_dim,lon_dim))
    precip_ens = np.zeros((13,year_dim, max_runs, lat_dim, lon_dim))
    
    
    #path = '/ouce-home/staff/sedm4922/Validation/'
    
    # factors to convert units to mm/month
    fac = 86400 # for the model
    fac2 = 30.5 # for the obs
    spatial = 'GHoA'
    # Define months array
    mon = ['01','02','03','04','05','06','07','08','09','10','11','12']
    
    #precip = np.zeros((12,25,146,209))
    #precip_ens = np.zeros((12,25,100,146,209))
    
    # get all the model run data into shape
    for year_ind in range(year_dim):
      selection = [el['filename'] for el in files_dict if el['end_year']==years[year_ind]]
      for run_ind in range(max_runs):
        file_path = os.path.join(dir_path, selection[run_ind])
        #print file_path
        ds = ncfile(path)
        ds_precip = ds.variables['item5216_monthly_mean'][:]
        ds_precip = np.squeeze(ds_precip)
        precip_ens[:,year_ind,run_ind,:,:] = ds_precip
    print 'precip', precip_ens.shape
    
    #average it by axis
    precip_output = np.nanmean(precip_ens, axis=4) #lat
    precip_output = np.nanmean(precip_output, axis=3) #lon
    precip_output = np.nanmean(precip_output, axis=1) #year -> 13 mos x N runs
    
    #get the average of all the runs
    precip_mean = np.nanmean(precip_output, axis=1)
    
    
    #filter GHOA extents
    min_lat = min(lats.flatten())
    max_lat = max(lats.flatten())
    min_lon = min(lons.flatten())
    max_lon = max(lons.flatten())
    
    #print min_lat, max_lat, min_lon, max_lon
    
    
    #load in the historic data
    subpath = r'NEWINPUT/'
    fname = r'chirps_clim.nc'
    file_path = os.path.join(subpath,fname)
    chirps_ds = ncfile(file_path)
    #print chirps_ds.dimensions.keys()
    #print chirps_ds.variables.keys()
    chirps_precip = chirps_ds.variables['precip'][:]
    
    #get the lats and lons and times
    lats = chirps_ds.variables['latitude'][:]
    lons = chirps_ds.variables['longitude'][:]
    times = chirps_ds.variables['time'][:]
    print times
    
    #create masks using the GHOA extents
    lat_mask = ((lats>min_lat)&(lats<max_lat))
    lon_mask = ((lons>min_lon)&(lons<max_lon))
    
    #tile one dimension of those badboys out
    tile = np.tile(lat_mask,(len(lons),1))
    
    #create a mask of same global shape
    chirps_mask = np.ones((12,len(lats),len(lons)), dtype=bool)
    
    #multiply the boolean tile and lon mask, broadcast into the global mask
    chirps_mask[:,:,:] = np.multiply(tile.T, lon_mask)
    print chirps_mask.shape
    
    #get the masked data, (convert to 'not ~') np.ma.where(data, mask)
    #chirps_ghoa = chirps_precip[~chirps_mask]
    #print chirps_ghoa.shape
    chirps_ghoa = np.ma.masked_where(chirps_precip,~chirps_mask)
    chirps_mean = np.nanmean(chirps_ghoa, axis = (1,2))
    
    #years, months, basemap
    
    #gen months
    print chirps_mean
    precip_mean = np.mean(precip_ens, axis = (3,4))
    precip_mean = np.delete(precip_mean, (0), axis=0)
    print precip_mean.shape
    print precip_mean
    shape = precip_mean.shape
    dats = np.random.rand(shape[0],shape[1],shape[2])
    #gr_monthly(chirps_mean, dats)
    gr_map(chirps_ghoa, precip_ens, min_lat, max_lat, min_lon, max_lon)

def gr_monthly(chirps_mean, data):
    fig, ax = plt.subplots()
    plt_chirps, = ax.plot(chirps_mean, label = 'CHIRPS', color = (0,0,0), linewidth = 3)
    shape = data.shape
    
    for i in range(shape[1]):
        for j in range(shape[2]):
            ax.plot(data[:,i,j], color = '#75bbfd')
    
    model = mlines.Line2D([], [], color='#75bbfd', label='MODEL HXXX')
    plt.legend(handles=[plt_chirps, model])
    #ax.set_axis_off()
    ax.set_xlim(0,11)
    ax.xaxis.set_ticks(np.arange(0,12))
    plt.ylabel('Monthly Precipitation [mm/month]')
    mos = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    ax.set_xticklabels(mos)

    plt.show()

def gr_yearly(chirps_mean, data):
    fix, ax = plt.subplots()

def gr_map(chirps_mean, data, min_lat, max_lat, min_lon, max_lon):
    print 'hi'
    print 'extents', min_lat, max_lat, min_lon, max_lon
    data = np.mean(data, axis = (1,2))[1,:,:]
    print data.shape
    print np.max(data.flatten())
    print np.min(data.flatten())
    m = Basemap(projection='merc',llcrnrlat=min_lat,urcrnrlat=max_lat,llcrnrlon=min_lon,urcrnrlon=max_lat,lat_ts=20,resolution='c')
    m.drawcoastlines()
    #m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    ny = data.shape[0]; nx = data.shape[1]
    lons, lats = m.makegrid(nx, ny)
    x, y = m(lons, lats)
    clevs = [0.0001,0.0002,0.0003,0.0004]
    cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
    m.drawmapboundary(fill_color='aqua')
    plt.title("Mercator Projection")
    plt.show()
    

if __name__ == "__main__":
    main()

"""
#rcp26 = dataset.variables['p50'][:]
"""
"""
       
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

"""
"""
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

