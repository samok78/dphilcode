import numpy as np
from netCDF4 import Dataset,MFDataset
from netCDF4 import Dataset as s
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm	
import glob,os,sys
from pandas import DataFrame,concat
workdir = '/ouce-home/staff/sedm4922/batch204/field90_annual/ENSMEAN/'
workdir2 = '/ouce-home/staff/sedm4922/Validation/CHIRPS/'
months = ['02','03','04','05','06','07','08','09']


fig = plt.figure(figsize =[8,10])
precip_mod=workdir+'batch204_clim.nc'
t02t=workdir2+'chirps_clim.nc'

precip_mod2=s(precip_mod, mode='r')
obs_mod2=s(t02t, mode='r')
precip_model=precip_mod2.variables['item5216_monthly_mean'][:,0,:,:]*86400
lat = precip_mod2.variables['latitude0'][:]
lon = precip_mod2.variables['longitude0'][:]
for i,m in enumerate(months):
	# Set up data for analysis
	#xdata = trmm_data[loc]
	xdata2 = chirps_data[loc]
	#ydata = data[1997].flatten()
	ydata=np.concatenate(data[loc].values()).flatten()
	print xdata.shape,ydata.shape
	df1=DataFrame(xdata,columns=['TRMM'])
	df2=DataFrame(ydata,columns=['weather@home 1998-2011'])
	df3=DataFrame(xdata2,columns=['CHIRPS'])

	# Calculate quantiles
	#quantiles=np.arange(0.005,1,.01)
	#quantiles=np.arange(0.025,1,.05)
	quantiles=np.arange(0.5,1,.005)
	# Value of data at quantiles
	q1=df1.quantile(quantiles)
	q2=df2.quantile(quantiles)
	q3=df3.quantile(quantiles)

	qvals=concat([q1,q2,q3],axis=1)

	sp = plt.subplot(4,3,i+1)
	qvals.plot.scatter('CHIRPS','weather@home 1998-2011',xlim=[0,max],ylim=[0,max])
	plt.title(loc+ ': Hadrm3p vs CHIRPS: quantiles .5 to .995')
	plt.xlabel('CHIRPS daily Nov precip (mm/day)')
	plt.ylabel('Hadrm3p daily Nov precip (mm/day)')
	plt.plot([0,max],[0,max],'r')
#	plt.show()
	plt.savefig('figs/'+loc+'_wah2_vs_chirps.png')
		
		# TRMM vs CHRIPS
		#

		# qq plot
	#qvals.plot.scatter('TRMM','CHIRPS',xlim=[0,max],ylim=[0,max])
	#plt.title(loc+ ': CHIRPS vs TRMM: quantiles .5 to .995')
	#plt.xlabel('TRMM daily Nov precip (mm/day)')
	#plt.ylabel('CHIRPS daily Nov precip (mm/day)')
	#plt.plot([0,max],[0,max],'r')
#	plt.show()
	#plt.savefig('figs/'+loc+'_chirps_vs_trmm.png')
