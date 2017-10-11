#Sarah test for github

import os
from sys import exit
import numpy as np
from netCDF4 import Dataset as s
#import Scientific.IO.NetCDF as s
import matplotlib.pyplot as plt
import numpy as numpy
import scipy.stats as sci
import pylab as pl

os.chdir('/ouce-home/staff/sedm4922')
workdir1='/ouce-home/staff/sedm4922/batch371/field90_mod/CAT/KARSTEN/SEASMN2'
workdir2='/ouce-home/staff/sedm4922/batch372/field90_mod/CAT/KARSTEN/SEASMN2'
workdir3='/ouce-home/staff/sedm4922/batch204/field90_mod/CAT/KARSTEN/SEASMN2'

workdir4='/ouce-home/staff/sedm4922/'

workdir5='/ouce-home/staff/sedm4922/batch_370/batch_370/ga.pe/field90_mod/SEAS'
workdir6='/ouce-home/staff/sedm4922/batch_367/batch_367/ga.pe/field90_mod/SEAS'
workdir7='/ouce-home/staff/sedm4922/batch387/field90_mod/CAT/KARSTEN/SEASMN2'

month = ['JJA']
fig = plt.figure()
c=0
for mo in month:
    c=c+1
    ifile_nat=workdir2+'/TEST4_372_'+mo+'.nc'
    ifile_all=workdir1+'/TEST4_371_'+mo+'.nc'
    ifile_obs=workdir4+'/OBS_O/CHIRPS/chirps_mon_mmday_'+mo+'_NE_FLD_SEAS.nc'
    ifile_clim=workdir3+'/TEST4_204_'+mo+'.nc'

    ifile_nat2=workdir5+'/TEST4_370_'+mo+'2.nc'
    ifile_all2=workdir6+'/TEST4_367_'+mo+'2.nc'
    ifile_clim2=workdir7+'/TEST4_387_'+mo+'.nc'


    d_nat=s(ifile_nat, mode='r')
    pr_nat=d_nat.variables['pr'][:]#*86400
    ensmember=d_nat.variables['ensmember'][:]

    d_all=s(ifile_all, mode='r')
    pr_all=d_all.variables['pr'][:]#*86400
    ensmember=d_all.variables['ensmember'][:]

    d_clim=s(ifile_clim, mode='r')
    pr_clim=d_clim.variables['pr'][:]#*86400
    ensmember=d_clim.variables['ensmember'][:]
    ax = fig.add_subplot(2,1,c)

    d_nat2=s(ifile_nat2, mode='r')
    pr_nat2=d_nat2.variables['pr'][:]#*86400
    time4=d_nat2.variables['ensmember'][:]

    d_all2=s(ifile_all2, mode='r')
    pr_all2=d_all2.variables['pr'][:]#*86400
    time5=d_all2.variables['ensmember'][:]

    d_clim2=s(ifile_clim2, mode='r')
    pr_clim2=d_clim2.variables['pr'][:]#*86400
    time6=d_clim2.variables['ensmember'][:]
    ax = fig.add_subplot(2,1,c)


    d_obs=s(ifile_obs, mode='r')
    pr_obs=d_obs.variables['precip'][:,0,0]

    e_all=np.sort(pr_all)#[::-1]
    e_nat=np.sort(pr_nat)#[::-1]
    e_clim=np.sort(pr_clim)#[::-1]
    e_obs=np.sort(pr_obs)

    e_all2=np.sort(pr_all2)#[::-1]
    e_nat2=np.sort(pr_nat2)#[::-1]
    e_clim2=np.sort(pr_clim2)#[::-1]
    #e_obs=np.sort(pr_obs)

    fit = sci.norm.pdf(e_all, np.mean(e_all), np.std(e_all))
    plt.plot(e_all,fit,'-o', label='ALL -AFRICA', alpha=0.3, fillstyle="none")
    #if you want histogram add following
    #plt.hist(e_all,normed=True, label='ALL', alpha=0.3)
    print e_all
    print fit
    print np.mean(e_all)
    print np.std(e_all)

    fit2 = sci.norm.pdf(e_nat, np.mean(e_nat), np.std(e_nat))
    plt.plot(e_nat,fit2,'-o', label='NATURAL _AFRICA', alpha=0.3, fillstyle="none")
    #plt.hist(e_nat,normed=True, label='NATURAL', alpha=0.3)

    fit3 = sci.norm.pdf(e_clim, np.mean(e_clim), np.std(e_clim))
    plt.plot(e_clim,fit3,'-o', label='HADAM3P 1987-2011 -AFRICA', alpha=0.3, fillstyle="none")
    #plt.hist(e_clim,normed=True, label='HADAM3P 1987-2011', alpha=0.3)

    fit4 = sci.norm.pdf(e_obs, np.mean(e_obs), np.std(e_obs))
    plt.plot(e_obs,fit4,'-o', label='CHIRPS 1987-2011', alpha=0.3, fillstyle="none")
    #plt.hist(e_obs,normed=True, label='CHIRPS', alpha=0.3)
    #plt.ylim([0, 1.2])

    fit5 = sci.norm.pdf(e_all2, np.mean(e_all2), np.std(e_all2))
    plt.plot(e_all2,fit5,'-o', label='ALL - ASIA', alpha=0.3, fillstyle="none")
    #if you want histogram add following
    #plt.hist(e_all2,normed=True, label='ALL', alpha=0.3)

    fit6 = sci.norm.pdf(e_nat2, np.mean(e_nat2), np.std(e_nat2))
    plt.plot(e_nat2,fit6,'-o', label='NATURAL- ASIA', alpha=0.3, fillstyle="none")
    #plt.hist(e_nat2,normed=True, label='NATURAL', alpha=0.3)

    fit7 = sci.norm.pdf(e_clim2, np.mean(e_clim2), np.std(e_clim2))
    plt.plot(e_clim2,fit7,'-o', label='HADAM3P 1987-2011 - ASIA', alpha=0.3, fillstyle="none")
    #plt.hist(e_clim2,normed=True, label='HADAM3P 1987-2011', alpha=0.3)

    plt.xlabel('mm/day',labelpad=5,fontsize=12)
    plt.ylabel('freq',labelpad=5,fontsize=12)
    plt.title(mo, fontsize=12)
    plt.legend(loc='best', fontsize=10)

fig.suptitle('NORM-PDF SEAS MEAN PRECIP: ETHIOPIA AFR VS AFRICAN DOMAIN', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.savefig(workdir4+'/IMAGES_O/HIST_PR_ETH_AFRICA.jpg', orientation = 'portrait',format='jpg', )
plt.show()
