#!/usr/bin/env mcsplot
# UNDERSTAND THE SENSITIVITY OF THE PMAX AND PTOT
# TRENDS TO BINNING

import sys,pickle,time
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm,LogNorm
import numpy.ma as ma
#import seaborn as sns
import pandas as pd

season = 'DJF'
enso   = 'endjf'
gebiet = [-33,33,0,360,'trop']
#execfile('../precip/MidpointNormalize.py')

# read in the values
#basedir = '/rigel/home/scs2229/top-secret/MCS_clim/'
basedir = '/work/bb1018/b380873/MCS_clim/'
fich1 = basedir + 'ausgabe/precip_clim/old-output/trop_depth_psum_pmax_endjf.npy'
fich2 = basedir + 'ausgabe/precip_clim/old-output/trop_depth_psum_pmax_lndjf.npy'
pEN = np.load(fich1); pLN = np.load(fich2)
    
psume = pEN[2,:]; ctte = pEN[0,:]; pmaxe = pEN[3,:]
psuml = pLN[2,:]; cttl = pLN[0,:]; pmaxl = pLN[3,:]
up2 = [15,20,25,30]
hi = 30 
cttlim = [240,215]

lbl = [r'CTT depth $\geq$ ' + str(cttlim[0]),str(cttlim[0]) + r' $>$ CTT depth $\geq$ ' + str(cttlim[1]),\
        str(cttlim[1]) + r' $>$ CTT depth']
farbe = ['blue','green','red']
let = ['(a)','(c)','(b)','(d)']
gogo = [r'max = 15 mm h$^{-1}$',r'max = 20 mm h$^{-1}$',\
       r'max = 25 mm h$^{-1}$',r'max = 30 mm h$^{-1}$']
axfont = 11
ticfont = 11
uu = 170
ll = -30
uu2 = 170
ll2 = -30
vm = 5

fig = plt.figure(figsize=(7.68,6.68))
for c in np.arange(0,2):
    ax = plt.subplot2grid((2,2),(c,0))
    xval2 = np.logspace(-1,np.log10(up2[c]),11)
    rd = np.zeros((3,10))
    for ii in np.arange(3):
        if ii == 0:
           jj = np.argwhere(ctte >= cttlim[0])
           kk = np.argwhere(cttl >= cttlim[0])
        elif ii == 1:
           jj = np.argwhere((ctte >= cttlim[1]) & (ctte < cttlim[0]))
           kk = np.argwhere((cttl >= cttlim[1]) & (cttl < cttlim[0]))
        else:
           jj = np.argwhere(ctte < cttlim[1])
           kk = np.argwhere(cttl < cttlim[1]) 
        hh, _ = np.histogram(pmaxe[jj],bins=xval2)
        hh = np.divide(hh,float(jj.shape[0]))
        gg, _ = np.histogram(pmaxl[kk],bins=xval2)
        gg = np.divide(gg,float(kk.shape[0]))
        rd[ii] = ma.masked_where(gg == 0.,(hh-gg)/gg*100.)
        rd[ii] = ma.masked_where(hh == 0.,rd[ii])
    ww = np.asarray([(xval2[i]-xval2[i-1])/2.5 for i in np.arange(1,len(xval2))])
    plt.bar(xval2[:-1],rd[0],width=ww,align='center',color=farbe[0],edgecolor='black',label=lbl[0])
    plt.bar(xval2[:-1],rd[1],width=ww,align='center',color=farbe[1],edgecolor='black',label=lbl[1])
    plt.bar(xval2[:-1],rd[2],width=ww,align='center',color=farbe[2],edgecolor='black',label=lbl[2])
    mask1 = np.argwhere((abs(rd[0]) < abs(rd[1])) | (abs(rd[0]) < abs(rd[2])))
    mask2 = np.argwhere((abs(rd[1]) < abs(rd[0])) | (abs(rd[1]) < abs(rd[2])))
    plt.bar(xval2[mask1[:,0]],rd[0,mask1[:,0]],width=ww[mask1[:,0]],color=farbe[0],edgecolor='black')
    ax.set_ylim([ll,uu])
    plt.text(0.05,0.9,let[c],fontsize=14,fontweight='bold',transform=ax.transAxes)
    plt.text(0.25,0.9,gogo[c],fontsize=11,transform=ax.transAxes)
    ax.set_xscale('log')
    ax.tick_params(axis='both',which='both',labelsize=ticfont)
    if c == 1:
       plt.xlabel(r'$\dot{P}$ [mm h$^{-1}$]',fontsize=axfont)
    if c == 0:
       plt.legend(loc='center left',fontsize=axfont-2)
fig.text(0.04,0.5,'Relative change in probability (EN-LN)/LN [%]',fontsize=axfont,va='center',rotation='vertical')

c = 0
for c in np.arange(0,2):
    ax = plt.subplot2grid((2,2),(c,1))
    xval2 = np.logspace(-1,np.log10(up2[c+2]),11)
    for ii in np.arange(3):
        if ii == 0:
           jj = np.argwhere(ctte >= cttlim[0])
           kk = np.argwhere(cttl >= cttlim[0])
        elif ii == 1:
           jj = np.argwhere((ctte >= cttlim[1]) & (ctte < cttlim[0]))
           kk = np.argwhere((cttl >= cttlim[1]) & (cttl < cttlim[0]))
        else:
           jj = np.argwhere(ctte < cttlim[1])
           kk = np.argwhere(cttl < cttlim[1])
        hh, _ = np.histogram(pmaxe[jj],bins=xval2)
        hh = np.divide(hh,float(jj.shape[0]))
        gg, _ = np.histogram(pmaxl[kk],bins=xval2)
        gg = np.divide(gg,float(kk.shape[0]))
        rd[ii] = ma.masked_where(gg == 0,(hh-gg)/gg*100.)
        rd[ii] = ma.masked_where(hh == 0,rd[ii])
    if rd[0,8] > 0:
       rd[0,8] = rd[0,8]*-1
    if rd[0,9] > 0:
       rd[0,9] = rd[0,9]*-1
    ww = np.asarray([(xval2[i]-xval2[i-1])/2. for i in np.arange(1,len(xval2))])
    plt.bar(xval2[:-1],rd[0],width=ww,color=farbe[0],align='center',edgecolor='black',label=lbl[0])
    plt.bar(xval2[:-1],rd[1],width=ww,color=farbe[1],align='center',edgecolor='black',label=lbl[1])
    plt.bar(xval2[:-1],rd[2],width=ww,color=farbe[2],align='center',edgecolor='black',label=lbl[2])
    mask1 = np.argwhere((abs(rd[0]) < abs(rd[1])) | (abs(rd[0]) < abs(rd[2])))
    mask2 = np.argwhere((abs(rd[1]) < abs(rd[0])) | (abs(rd[1]) < abs(rd[2])))
    plt.bar(xval2[mask1[:,0]],rd[0,mask1[:,0]],width=ww[mask1[:,0]],color=farbe[0],edgecolor='black')
    plt.text(0.05,0.9,let[c+2],fontsize=14,fontweight='bold',transform=ax.transAxes)
    plt.text(0.25,0.9,gogo[c+2],fontsize=11,transform=ax.transAxes)
    ax.set_xscale('log')
    ax.set_ylim([ll,uu])
    ax.tick_params(axis='both',which='both',labelsize=ticfont)
    if c == 1:
       plt.xlabel(r'$\dot{P}$ [mm h$^{-1}$]',fontsize=axfont)

fig.subplots_adjust(hspace=0.25)
fig.subplots_adjust(wspace=0.25)
fig.savefig('./figures/figS5-bounds-pmax.pdf',bbox_inches='tight')
plt.show()

