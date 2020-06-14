#!/usr/bin/env mcsplot
import sys,pickle,time
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm,LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from scipy.io import loadmat
import scipy.stats
import scipy.integrate as integrate

def mean_confidence_interval(data,confidence):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.median(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return h

fs = 13
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

fig = plt.figure(figsize=(9,6))
#basedir = '/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/temperature/core/'
#basedir2 = '/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/qv/core/'
basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/meteo_clim/temperature/core/'
basedir2 = '/work/bb1018/b380873/MCS_clim/ausgabe/meteo_clim/qv/core/'

ENDJFqv1 = np.load(basedir2 + 'qv_ENDJF1.npy')
LNDJFqv1 = np.load(basedir2 + 'qv_LNDJF1.npy')
ENDJFqv2 = np.load(basedir2 + 'qv_ENDJF2.npy')
LNDJFqv2 = np.load(basedir2 + 'qv_LNDJF2.npy')
ENDJFqv3 = np.load(basedir2 + 'qv_ENDJF3.npy')
LNDJFqv3 = np.load(basedir2 + 'qv_LNDJF3.npy')

ENDJFT1 = np.load(basedir + 'T_ENDJF1.npy')
LNDJFT1 = np.load(basedir + 'T_LNDJF1.npy')
ENDJFT2 = np.load(basedir + 'T_ENDJF2.npy')
LNDJFT2 = np.load(basedir + 'T_LNDJF2.npy')
ENDJFT3 = np.load(basedir + 'T_ENDJF3.npy')
LNDJFT3 = np.load(basedir + 'T_LNDJF3.npy')

# create a function to calculate saturation mixing ratio from temperature
def satMR(T):
    a1 = 54.842763; a2 = -6763.22; a3 = -4.21; a4 = 0.000367
    a5 = 0.0415; a6 = 218.8; a7 = 53.878; a8 = -1331.22
    a9 = -9.44523; a10 = 0.014025
    MWw = 0.01802
    MWa = 0.02897

    psat = a1 + a2/T + a3*np.log(T) + a4*T + np.arctan(a5*(T - a6))*(a7 + a8/T + \
           a9*np.log(T) + a10*T)
    psat = np.exp(psat)
    qvcsat = MWw/MWa*psat/(press*100. - psat)
    return qvcsat

RHavg = np.zeros((6,32)); Tavg = np.zeros((6,32))
qvsat = np.zeros((6,32))

# MEAN RH PROFILES
Tavg[0] = np.nanmean(ENDJFT1,axis=0)
Tavg[1] = np.nanmean(ENDJFT2,axis=0)
Tavg[2] = np.nanmean(ENDJFT3,axis=0)
Tavg[3] = np.nanmean(LNDJFT1,axis=0)
Tavg[4] = np.nanmean(LNDJFT2,axis=0)
Tavg[5] = np.nanmean(LNDJFT3,axis=0)
qvsat = satMR(Tavg)

RHavg[0] = np.nanmean(ENDJFqv1,axis=0)/qvsat[0]
RHavg[1] = np.nanmean(ENDJFqv2,axis=0)/qvsat[1]
RHavg[2] = np.nanmean(ENDJFqv3,axis=0)/qvsat[2]
RHavg[3] = np.nanmean(LNDJFqv1,axis=0)/qvsat[3]
RHavg[4] = np.nanmean(LNDJFqv2,axis=0)/qvsat[4]
RHavg[5] = np.nanmean(LNDJFqv3,axis=0)/qvsat[5]

# RH DIFFERENCE PROFILES
y1 = 1002; y2 = 122
ax2 = plt.subplot2grid((1,2),(0,0))
ax2.plot((RHavg[0]-RHavg[3])/RHavg[3]*100.,press,color='blue',linewidth=1.25)
ax2.plot((RHavg[1]-RHavg[4])/RHavg[4]*100.,press,color='green',linewidth=1.25)
ax2.plot((RHavg[2]-RHavg[5])/RHavg[5]*100.,press,color='red',linewidth=1.25)
ax2.plot([0,0],[y2,y1],color='black',linewidth=0.75,linestyle='--')
plt.text(0.05,0.92,'(a)',fontsize=fs+3,fontweight='bold',transform=ax2.transAxes)
plt.text(0.55,0.64,'Deepest',fontsize=fs,color='red',transform=fig.gca().transAxes)
plt.text(0.55,0.56,'Intermediate',fontsize=fs,color='green',transform=fig.gca().transAxes)
plt.text(0.55,0.48,'Least deep',fontsize=fs,color='blue',transform=fig.gca().transAxes)
plt.ylim([y2,y1]); plt.xlim([-25.,25.])
plt.xlabel('(EN-LN)/LN relative difference'
            '\n'
            'in mean RH [%]',fontsize=fs+2)
ax2.tick_params(axis='both',labelsize=fs+1)
plt.ylabel('Pressure [hPa]',fontsize=fs+3)
plt.gca().invert_yaxis()
plt.yticks([200,400,600,800,1000])

# ZBP CALCS
ax3 = plt.subplot2grid((2,2),(0,1))
data = loadmat('article9-reviews/buoyzbp_1ent_coincid_p99.mat')
buoy = data['buoy']
buoyT = data['buoyT']
buoyRH = data['buoyRH']
#95% confidence interval
data2 = loadmat('article9-reviews/buoyzbp_1ent_coincid_CI95.mat')
ci95 = data2['yCI95']
#ax3.fill_betweenx(press,buoy[0]-ci95[:,0],buoy[0]+ci95[:,0],color='blue',alpha=0.25)
#plt.plot(buoy[0],press,linewidth=1.25,color='blue')
#plt.plot(buoy[3],press,linewidth=1.25,color='blue',linestyle='--')
print(buoy[2]/buoy[5])
ax3.fill_betweenx(press,buoy[2]-ci95[:,2],buoy[2]+ci95[:,2],color='red',alpha=0.25)
ax3.fill_betweenx(press,buoy[5]-ci95[:,5],buoy[5]+ci95[:,5],color='blue',alpha=0.25)
#plt.plot(buoy[1],press,linewidth=1.25,color='green')
#plt.plot(buoy[4],press,linewidth=1.25,color='green',linestyle='--')
plt.plot(buoy[2],press,linewidth=1.25,color='red')
plt.plot(buoy[5],press,linewidth=1.25,color='blue')
plt.text(0.05,0.9,'(b)',fontsize=fs+3,fontweight='bold',transform=ax3.transAxes)
plt.text(0.5,0.1,'Deepest only',style='italic',fontsize=fs+1,transform=ax3.transAxes)
plt.xlabel(r'Buoyancy ($\dot{P}_{99}$) [m s$^{-2}$]',fontsize=fs+1)
plt.ylabel('Pressure [hPa]',fontsize=fs+1)
ax3.tick_params(axis='both',labelsize=fs)
plt.yticks([400,600,800,1000])
plt.ylim([300,1000])
plt.xlim([0,0.6])
plt.gca().invert_yaxis()

# ZBP CALC DIFFERENCES
ax4 = plt.subplot2grid((2,2),(1,1))
ax4.tick_params(axis='both',labelsize=fs)
#plt.plot(buoy[0]-buoy[3],press,linewidth=1.25,color=[0,0,1],linestyle='-')
#plt.plot(buoyRH[0]-buoyRH[3],press,linewidth=1.25,color='cyan',linestyle='--')
#plt.plot(buoyT[0]-buoyT[3],press,linewidth=1.25,color=[0,0,0.5],linestyle='-.')

plt.plot(buoyRH[2]-buoyRH[5],press,linewidth=1.25,color='pink',linestyle='--')
plt.plot(buoy[2]-buoy[5]-(buoyRH[2]-buoyRH[5]),press,linewidth=1.25,color='maroon',linestyle='-.')
plt.plot(buoy[2]-buoy[5],press,linewidth=1.25,color='red',linestyle='-')
#print buoyRH[2]-buoyRH[5]
#print buoyT[2]-buoyT[5]

#plt.plot([0,0],[300,1000],linewidth=0.75,color='black',linestyle='--')
plt.text(0.17,0.9,'(c)',fontsize=fs+3,fontweight='bold',transform=ax4.transAxes)
plt.text(0.14,750,'Combined',fontsize=fs,color='red')
plt.text(0.14,850,'Effect of RH',fontsize=fs,color='maroon')
plt.text(0.14,950,'Effect of T',fontsize=fs,color='pink')
plt.xlabel(r'El Ni$\~n$o-La Ni$\~n$a $\Delta$ Buoyancy ($\dot{P}_{99}$) [m s$^{-2}$]',fontsize=fs+1)
plt.ylabel('Pressure [hPa]',fontsize=fs+1)
#r1 = Rectangle((0.15,700),0.04,20,facecolor='blue')
#r3 = Rectangle((0.18,700),0.04,20,facecolor='red')
#r4 = Rectangle((0.15,800),0.04,20,facecolor=[0,0,0.5])
#r6 = Rectangle((0.18,800),0.04,20,facecolor=[0.5,0,0])
#r7 = Rectangle((0.15,900),0.04,20,facecolor='cyan')
#r9 = Rectangle((0.18,900),0.04,20,facecolor=[1,0.6,0.4])
#ax4.add_patch(r1); #ax4.add_patch(r2); 
#ax4.add_patch(r3)
#ax4.add_patch(r4); #ax4.add_patch(r5); 
#ax4.add_patch(r6)
#ax4.add_patch(r7); #ax4.add_patch(r8); 
#ax4.add_patch(r9)
#plt.ylim([300,1000])
plt.xlim([-0.02,0.23])
plt.gca().invert_yaxis()
#plt.yticks([400,600,800,1000])

#fig.subplots_adjust(wspace=0.2)
plt.tight_layout()
#fig.savefig('./figures/fig4-rh_zbp_ci95.pdf',bbox_inches='tight')
plt.show()
