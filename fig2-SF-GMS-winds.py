import sys,pickle,time
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm,LogNorm
from netCDF4 import Dataset,num2date
import pandas as pd

# if you want to regenerate the El Nino, La Nina, and climatological mean profiles
# makedata should be set to true
basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/meteo_clim/zonal_winds/core/'
basedir2 = '/work/bb1018/b380873/MCS_clim/ausgabe/meteo_clim/merid_winds/core/'

ENDJFuu1 = np.load(basedir + 'u_ENDJF1.npy')
LNDJFuu1 = np.load(basedir + 'u_LNDJF1.npy')
ENDJFuu2 = np.load(basedir + 'u_ENDJF2.npy')
LNDJFuu2 = np.load(basedir + 'u_LNDJF2.npy')
ENDJFuu3 = np.load(basedir + 'u_ENDJF3.npy')
LNDJFuu3 = np.load(basedir + 'u_LNDJF3.npy')

ENDJFvv1 = np.load(basedir2 + 'v_ENDJF1.npy')
LNDJFvv1 = np.load(basedir2 + 'v_LNDJF1.npy')
ENDJFvv2 = np.load(basedir2 + 'v_ENDJF2.npy')
LNDJFvv2 = np.load(basedir2 + 'v_LNDJF2.npy')
ENDJFvv3 = np.load(basedir2 + 'v_ENDJF3.npy')
LNDJFvv3 = np.load(basedir2 + 'v_LNDJF3.npy')

ENwind1 = np.sqrt(ENDJFuu1**2 + ENDJFvv1**2)
LNwind1 = np.sqrt(LNDJFuu1**2 + LNDJFvv1**2)
ENwind2 = np.sqrt(ENDJFuu2**2 + ENDJFvv2**2)
LNwind2 = np.sqrt(LNDJFuu2**2 + LNDJFvv2**2)
ENwind3 = np.sqrt(ENDJFuu3**2 + ENDJFvv3**2)
LNwind3 = np.sqrt(LNDJFuu3**2 + LNDJFvv3**2)


# create the average profile and its spread for all sets
alt = np.asarray([15418.43,14415.44,13470.22,12579.42,11739.87,10948.81,10191.54,\
       9458.63,8749.51,8064.40,7403.97,6769.08,6160.62,5579.46,5026.39,4502.09,\
       4007.11,3541.93,3106.86,2702.14,2327.89,1984.11,1670.7,1387.43,1133.93,\
       909.7,714.05,546.11,404.74,288.57,195.85,124.48])
press = np.asarray([122.6137,142.9017,165.0886,189.1466,215.0251,242.6523,272.0593,\
         303.2174,336.0439,370.4072,406.1328,443.0086,480.7907,519.2093,557.9734,\
         596.7774,635.3060,673.2403,710.2627,746.0635,780.3455,812.8303,843.2634,\
         871.4203,897.1118,920.1893,940.5511,958.1477,972.9868,985.1399,994.7472,\
         1002.0236])

windavg = np.zeros((6,32))
windstd = np.zeros((6,32))
windhi = np.zeros((6,32))

# WIND AVERAGES
windavg[0] = np.nanmean(ENwind1,axis=0)
windavg[1] = np.nanmean(ENwind2,axis=0)
windavg[2] = np.nanmean(ENwind3,axis=0)
windavg[3] = np.nanmean(LNwind1,axis=0)
windavg[4] = np.nanmean(LNwind2,axis=0)
windavg[5] = np.nanmean(LNwind3,axis=0)

# WIND STANDARD DEVIATIONS
windstd[0] = np.nanstd(ENwind1,axis=0)
windstd[1] = np.nanstd(ENwind2,axis=0)
windstd[2] = np.nanstd(ENwind3,axis=0)
windstd[3] = np.nanstd(LNwind1,axis=0)
windstd[4] = np.nanstd(LNwind2,axis=0)
windstd[5] = np.nanstd(LNwind3,axis=0)

# ZONAL WIND EXTREMES
windhi[0] = np.percentile(ENwind1,99.,axis=0)   # easterly
windhi[1] = np.percentile(ENwind2,99.,axis=0)
windhi[2] = np.percentile(ENwind3,99.,axis=0)
windhi[3] = np.percentile(LNwind1,99.,axis=0)
windhi[4] = np.percentile(LNwind2,99.,axis=0)
windhi[5] = np.percentile(LNwind3,99.,axis=0)

fs = 13
fig = plt.figure(figsize=(9.3,5.76))

# MEAN WIND PROFILES
ax1 = plt.subplot2grid((1,3),(0,1))
ax1.plot((windavg[0]-windavg[3])/windavg[3]*100.,press,color='blue',linewidth=1.25,label=r'Least deep')
ax1.plot((windavg[1]-windavg[4])/windavg[4]*100.,press,color='green',linewidth=1.25,label=r'Intermediate')
ax1.plot((windavg[2]-windavg[5])/windavg[5]*100.,press,color='red',linewidth=1.25,label=r'Deepest')
ax1.plot([0.,0.],[0.,15.],color='black',linewidth=0.75,linestyle='--')
plt.text(0.05,0.92,'(c)',fontsize=fs,fontweight='bold',transform=ax1.transAxes)
plt.ylim([300,1000.])
ax1.invert_yaxis()
plt.xlabel(r'(EN-LN)/LN relative difference in '
           '\n'
           'mean wind speed [m s$^{-1}$]',fontsize=fs)
plt.ylabel('Pressure [hPa]',fontsize=fs)
#ax1.update(wspace=0.1,hspace=0.1)
plt.tight_layout()
ax1.tick_params(axis='both',labelsize=fs-1)
#ax1.legend(loc='center left',fontsize=fs-2)
plt.text(0.05,0.8,'Deepest',fontsize=fs-1,transform=ax1.transAxes,color='red')
plt.text(0.05,0.7,'Intermediate',fontsize=fs-1,transform=ax1.transAxes,color='green')
plt.text(0.05,0.6,'Least deep',fontsize=fs-1,transform=ax1.transAxes,color='blue')

# Load the GMS and surface flux values.
basedir = '/work/bb1018/b380873/MCS_clim/scripts/figs/GMS/'
ENd1 = np.load(basedir + 'GMS_ENDJF1.npy')   #_pmax_gt_1.npy')
LNd1 = np.load(basedir + 'GMS_LNDJF1.npy')
ENd3 = np.load(basedir + 'GMS_ENDJF3_pmax_gt_10.npy')
LNd3 = np.load(basedir + 'GMS_LNDJF3_pmax_gt_10.npy')

ax2 = plt.subplot2grid((2,3),(0,2))
d = 0.1; u = 1.3; n = 13;
h1,bar_edges = np.histogram(ENd1,bins=np.linspace(d,u,n),density=True)
h2,_ = np.histogram(LNd1,bins=np.linspace(d,u,n),density=True)
bar_center = np.asarray([(bar_edges[i] + bar_edges[i+1])/2 for i in np.arange(len(bar_edges)-1)])
ww = np.asarray([bar_edges[i]-bar_edges[i-1] for i in np.arange(1,len(bar_edges))])
l1 = ENd1.shape[0]
l2 = LNd1.shape[0]
plt.ylim([-40,35])
plt.bar(bar_center,(h1-h2)/h2*100.,width=ww,color='purple',edgecolor='black')
ax2.tick_params(axis='both',labelsize=fs-1)
plt.text(0.05,0.92,'(d) Least deep',fontsize=13,fontweight='bold',transform=ax2.transAxes)
plt.text(0.05,0.8,'KL = 0.0453 / 0.0471',fontsize=fs-1,transform=ax2.transAxes)
ax2.set_ylabel('$\Delta p$ (EN-LN)/LN [%]',fontsize=fs)

ax3 = plt.subplot2grid((2,3),(1,2))
h1,_ = np.histogram(ENd3,bins=np.linspace(d,u,n),density=True)
h2,_ = np.histogram(LNd3,bins=np.linspace(d,u,n),density=True)
l1 = ENd3.shape[0]
l2 = LNd3.shape[0]
plt.ylim([-40,35])
plt.bar(bar_center,(h1-h2)/h2*100.,width=ww,color='purple',edgecolor='black')
plt.text(0.05,0.4,'(e) Deepest',fontsize=13,fontweight='bold',transform=ax3.transAxes)
plt.text(0.05,0.2,'KL = 0.0537 / 0.0538',fontsize=fs-1,transform=ax3.transAxes)
ax3.tick_params(axis='both',labelsize=fs-1)
ax3.set_ylabel('$\Delta p$ (EN-LN)/LN [%]',fontsize=fs)
ax3.set_xlabel(r'NGMS ($\dot{P}_{99}$)',fontsize=fs)

basedir = '/work/bb1018/b380873/MCS_clim/scripts/figs/SF/'
ENd1 = (np.load(basedir + 'LH_ENDJF1_pmax_gt_1.npy') + np.load(basedir + 'SH_ENDJF1_pmax_gt_1.npy'))/(-3*3600)
LNd1 = (np.load(basedir + 'LH_LNDJF1_pmax_gt_1.npy') + np.load(basedir + 'SH_LNDJF1_pmax_gt_1.npy'))/(-3*3600)
ENd3 = (np.load(basedir + 'LH_ENDJF3.npy') + np.load(basedir + 'SH_ENDJF3.npy'))/(-3*3600) #_pmax_gt_10.npy
LNd3 = (np.load(basedir + 'LH_LNDJF3.npy') + np.load(basedir + 'SH_LNDJF3.npy'))/(-3*3600)

ax4 = plt.subplot2grid((2,3),(0,0))
d = 10; u = 1200; n = 13;
h1,bar_edges = np.histogram(ENd1,bins=np.linspace(d,u,n),density=True)
h2,_ = np.histogram(LNd1,bins=np.linspace(d,u,n),density=True)
bar_center = np.asarray([(bar_edges[i] + bar_edges[i+1])/2 for i in np.arange(len(bar_edges)-1)])
ww = np.asarray([bar_edges[i]-bar_edges[i-1] for i in np.arange(1,len(bar_edges))])
plt.bar(bar_center,(h1-h2)/h2*100.,width=ww,color='orange',edgecolor='black')
plt.ylim([-40,35])
plt.text(0.05,0.92,'(a) Least deep',fontsize=13,fontweight='bold',transform=ax4.transAxes)
plt.text(0.05,0.83,'KL = 0.0524 / 0.0593',fontsize=fs-1,transform=ax4.transAxes)
ax4.set_ylabel('$\Delta p$ (EN-LN)/LN [%]',fontsize=fs)
ax4.tick_params(axis='both',labelsize=fs-1)

ax5 = plt.subplot2grid((2,3),(1,0))
h1,_ = np.histogram(ENd3,bins=np.linspace(d,u,n),density=True)
h2,_ = np.histogram(LNd3,bins=np.linspace(d,u,n),density=True)
plt.bar(bar_center,(h1-h2)/h2*100.,width=ww,color='orange',edgecolor='black')
plt.text(0.05,0.8,'KL = 0.0278 / 0.0271',fontsize=fs-1,transform=ax5.transAxes)
plt.ylim([-40,35])
plt.text(0.05,0.92,'(b) Deepest',fontsize=13,fontweight='bold',transform=ax5.transAxes)
ax5.set_ylabel('$\Delta p$ (EN-LN)/LN [%]',fontsize=fs)
ax5.set_xlabel(r'Surface flux ($\dot{P}_{99}$) [W m$^{-2}$]',fontsize=fs)
ax5.tick_params(axis='both',labelsize=fs-1)

fig.savefig('./figures/fig2-SF-GMS-winds.pdf',bbox_inches='tight')
plt.show()
