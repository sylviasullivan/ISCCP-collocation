import sys,pickle,time
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm,LogNorm
from netCDF4 import Dataset,num2date
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.io import loadmat
import scipy.stats
import scipy.integrate as integrate

basedir3 = '/work/bb1018/b380873/MCS_clim/scripts/figs/drag/'

# if you want to regenerate the El Nino, La Nina, and climatological mean profiles
# makedata should be set to true
def mean_confidence_interval(data,confidence):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.median(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return h

makedata = False
filterdata = False
basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/meteo_clim/updrafts/other/'
basedir2 = '/work/bb1018/b380873/MCS_clim/ausgabe/meteo_clim/divergence/core/'

ENDJFww1 = np.load(basedir + 'w_ENDJF1.npy')
LNDJFww1 = np.load(basedir + 'w_LNDJF1.npy')
ENDJFww2 = np.load(basedir + 'w_ENDJF2.npy')
LNDJFww2 = np.load(basedir + 'w_LNDJF2.npy')
ENDJFww3 = np.load(basedir + 'w_ENDJF3.npy')
LNDJFww3 = np.load(basedir + 'w_LNDJF3.npy')
fs = 11

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
wavg = np.zeros((6,32)); whi = np.zeros((6,32)); wlo = np.zeros((6,32))

# UPDRAFT AVERAGES
wavg[0] = np.nanmean(ENDJFww1,axis=0)
wavg[1] = np.nanmean(ENDJFww2,axis=0)
wavg[2] = np.nanmean(ENDJFww3,axis=0)
wavg[3] = np.nanmean(LNDJFww1,axis=0)
wavg[4] = np.nanmean(LNDJFww2,axis=0)
wavg[5] = np.nanmean(LNDJFww3,axis=0)

# UPDRAFT EXTREMES
whi[0] = np.percentile(ENDJFww1,0.001,axis=0)
whi[1] = np.percentile(ENDJFww2,0.001,axis=0)
whi[2] = np.percentile(ENDJFww3,0.001,axis=0)
whi[3] = np.percentile(LNDJFww1,0.001,axis=0)
whi[4] = np.percentile(LNDJFww2,0.001,axis=0)
whi[5] = np.percentile(LNDJFww3,0.001,axis=0)

wlo[0] = np.percentile(ENDJFww1,99.99,axis=0)
wlo[1] = np.percentile(ENDJFww2,99.99,axis=0)
wlo[2] = np.percentile(ENDJFww3,99.99,axis=0)
wlo[3] = np.percentile(LNDJFww1,99.99,axis=0)
wlo[4] = np.percentile(LNDJFww2,99.99,axis=0)
wlo[5] = np.percentile(LNDJFww3,99.99,axis=0)

fig = plt.figure(figsize=(9.25,10))

# VELOCITY MEAN PROFILES
y1 = 1002; y2 = 122
ax1 = plt.subplot2grid((2,2),(0,0))
ax1.plot(wavg[0],press,color='blue',linewidth=1.25)
ax1.plot(wavg[1],press,color='green',linewidth=1.25,label=r'EN - 65 K $<$ depth $\leq$ 85 K')
ax1.plot(wavg[2],press,color='red',linewidth=1.25,label=r'EN - depth $>$ 85 K')
ax1.plot(wavg[3],press,color='blue',linestyle='--',linewidth=1.25,label=r'LN - depth $\leq$ 65 K')
ax1.plot(wavg[4],press,color='green',linestyle='--',linewidth=1.25,label=r'LN - 65 K $<$ depth $\leq$ 85 K')
ax1.plot(wavg[5],press,color='red',linestyle='--',linewidth=1.25,label=r'LN - depth $>$ 85 K')
ax1.plot([0.,0.],[y2,y1],color='black',linewidth=0.75,linestyle='--')
plt.text(0.05,0.92,'(a)',fontsize=fs+3,fontweight='bold',transform=ax1.transAxes)
plt.xlim([-0.02,0.02])
plt.ylim([y2,y1])
plt.xlabel(r'Mean pressure velocity [Pa s$^{-1}$]',fontsize=fs+3)
plt.ylabel('Pressure [hPa]',fontsize=fs+3)
ax1.tick_params(axis='both',labelsize=fs+1)
plt.gca().invert_yaxis()
plt.xticks([-0.02,-0.01,0,0.01,0.02])
plt.yticks([200,400,600,800,1000])

# VELOCITY EXTREME PROFILES
ax2 = plt.subplot2grid((2,2),(0,1))
ax2.plot(whi[0],press,color='blue',linewidth=1.25,label='EN - depth $\leq$ 65 K - 99th')
ax2.plot(whi[1],press,color='green',linewidth=1.25,label='EN - 65 K $<$ depth $\leq$ 85 K - 99th')
ax2.plot(whi[2],press,color='red',linewidth=1.25,label='EN - depth $>$ 85 K - 99th')
ax2.plot(whi[3],press,color='blue',linewidth=1.25,linestyle='--',label='LN - depth $\leq$ 65 K - 99th')
ax2.plot(whi[4],press,color='green',linewidth=1.25,linestyle='--',label='LN - 65 K $<$ depth $\leq$ 85 K - 99th')
ax2.plot(whi[5],press,color='red',linewidth=1.25,linestyle='--',label='LN - depth $>$ 85 K - 99th')

ax2.tick_params(axis='both',labelsize=fs+1)
plt.ylim([y2,y1])
plt.xlim([-5,0.25])
plt.text(0.05,0.92,'(b)',fontsize=fs+3,fontweight='bold',transform=ax2.transAxes)
plt.text(0.65,0.3,'least deep',color='blue',fontsize=fs,transform=ax2.transAxes)
plt.text(0.65,0.4,'intermediate',color='green',fontsize=fs,transform=ax2.transAxes)
plt.text(0.65,0.5,'deepest',color='red',fontsize=fs,transform=ax2.transAxes)
plt.xlabel(r'Extreme pressure velocity [Pa s$^{-1}$]',fontsize=fs+3)
plt.gca().invert_yaxis()
plt.yticks([200,400,600,800,1000])

# PROBABILITY DISTRIBUTION OF PRESSURE DRAG FOR LEAST DEEP SYSTEMS
ax3 = plt.subplot2grid((2,2),(1,0))
u = 4000
d = 7500
dragEN = np.load(basedir3 + 'dragintEN_d1.npy')
dragLN = np.load(basedir3 + 'dragintLN_d1.npy')
wgts = np.ones_like(dragEN)/float(len(dragEN))*100
h1,e1 = np.histogram(dragEN,bins=np.linspace(u,d,25),weights=wgts)
c1 = [(e1[i]+e1[i+1])/2. for i in np.arange(len(e1)-1)]
ww = (d-u)/25.5
plt.bar(c1,h1,color='red',edgecolor='k',align='center',width=ww,alpha=0.75)
wgts = np.ones_like(dragLN)/float(len(dragLN))*100
h2,e2 = np.histogram(dragLN,bins=np.linspace(u,d,25),weights=wgts)
c2 = [(e2[i]+e2[i+1])/2. for i in np.arange(len(e2)-1)]
plt.bar(c2,h2,color='blue',edgecolor='k',align='center',width=ww)
mask = np.argwhere(h1 < h2)
mask = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20])
mask.ravel()
juju = np.asarray([c1[m] for m in mask])
jojo = np.asarray([h1[m] for m in mask])
plt.bar(juju,jojo,color='red',edgecolor='k',align='center',width=ww,alpha=0.75)
plt.plot([np.nanmean(dragEN),np.nanmean(dragEN)],[0,18],color='red',linewidth=1.75,linestyle='-')
o = mean_confidence_interval(dragEN,0.99)
print np.nanmean(dragEN),np.nanmean(dragLN)
print (np.nanmean(dragEN)-np.nanmean(dragLN))/np.nanmean(dragEN)*100
plt.plot([np.nanmean(dragEN)-o,np.nanmean(dragEN)-o],[0,17],color='red',linewidth=1.25,linestyle='--')
plt.plot([np.nanmean(dragEN)+o,np.nanmean(dragEN)+o],[0,17],color='red',linewidth=1.25,linestyle='--')
plt.plot([np.nanmean(dragLN),np.nanmean(dragLN)],[0,18],color='blue',linewidth=1.75,linestyle='-')
o = mean_confidence_interval(dragLN,0.99)
plt.plot([np.nanmean(dragLN)-o,np.nanmean(dragLN)-o],[0,17],color='blue',linewidth=1.25,linestyle='--')
plt.plot([np.nanmean(dragLN)+o,np.nanmean(dragLN)+o],[0,17],color='blue',linewidth=1.25,linestyle='--')
ax3.tick_params(axis='both',labelsize=fs+1)
ax3.set_ylim([0,20])
plt.text(0.05,0.92,'(c)',fontsize=fs+3,fontweight='bold',transform=ax3.transAxes)
plt.ylabel('Probability',fontsize=fs+3)
plt.xlim([u,d])

# PROBABILITY DISTRIBUTION OF PRESSURE DRAG FOR DEEPEST SYSTEMS
ax4 = plt.subplot2grid((2,2),(1,1))
u = 6000
d = 9000
dragEN = np.load(basedir3 + 'dragintEN_d3_pmax_gt_10.npy')
dragLN = np.load(basedir3 + 'dragintLN_d3_pmax_gt_10.npy')
wgts = np.ones_like(dragEN)/float(len(dragEN))*100
h1,e1 = np.histogram(dragEN,bins=np.linspace(u,d,25),weights=wgts)
c1 = [(e1[i]+e1[i+1])/2. for i in np.arange(len(e1)-1)]
ww = (d-u)/25.5
plt.bar(c1,h1,color='red',edgecolor='k',align='center',width=ww,alpha=0.75)
wgts = np.ones_like(dragLN)/float(len(dragLN))*100
h2,e2 = np.histogram(dragLN,bins=np.linspace(u,d,25),weights=wgts)
c2 = [(e2[i]+e2[i+1])/2. for i in np.arange(len(e2)-1)]
plt.bar(c2,h2,color='blue',edgecolor='k',align='center',width=ww,alpha=0.75)
mask = np.argwhere(h1 < h2)
mask = np.array([11,12,13,14,15,16,17,18,19,20,21,22,23])
mask.ravel()
juju = np.asarray([c1[m] for m in mask])
jojo = np.asarray([h1[m] for m in mask])
plt.bar(juju,jojo,color='red',edgecolor='k',align='center',width=ww)
plt.plot([np.nanmean(dragEN),np.nanmean(dragEN)],[0,18],color='red',linewidth=1.75,linestyle='-')
o = mean_confidence_interval(dragEN,0.99)
print np.nanmean(dragEN),np.nanmean(dragLN)
print (np.nanmean(dragEN)-np.nanmean(dragLN))/np.nanmean(dragEN)*100
plt.plot([np.nanmean(dragEN)-o,np.nanmean(dragEN)-o],[0,17],color='red',linewidth=1.25,linestyle='--')
plt.plot([np.nanmean(dragEN)+o,np.nanmean(dragEN)+o],[0,17],color='red',linewidth=1.25,linestyle='--')
plt.plot([np.nanmean(dragLN),np.nanmean(dragLN)],[0,18],color='blue',linewidth=1.75,linestyle='-')
o = mean_confidence_interval(dragLN,0.99)
plt.plot([np.nanmean(dragLN)-o,np.nanmean(dragLN)-o],[0,17],color='blue',linewidth=1.25,linestyle='--')
plt.plot([np.nanmean(dragLN)+o,np.nanmean(dragLN)+o],[0,17],color='blue',linewidth=1.25,linestyle='--')
ax4.tick_params(axis='both',labelsize=fs+1)
ax4.set_ylim([0,20])
plt.text(0.05,0.92,'(d)',fontsize=fs+3,fontweight='bold',transform=ax4.transAxes)
fig.text(0.5, 0.04, r'Integrated pressure gradient force ($\dot{P}_{99}$) [m$^2$ s$^{-2}$]', ha='center',\
   fontsize=fs+2)
plt.xlim([u,d])

fig.subplots_adjust(wspace=0.2)
fig.savefig('./figures/fig3-upd-pgf.pdf',bbox_inches='tight')
plt.show()
