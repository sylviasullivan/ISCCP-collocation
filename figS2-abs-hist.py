#!/usr/bin/env mcsplot
# CREATE THE DIFFERENCE IN JOINT DISTRIBUTIONS OF MCS SIZE AND
# PTOT (PRECIPITATION ACCUMULATION) OR PMAX AND PTOT FOR A CERTAIN
# SUBDOMAIN AND SEASON.

import sys,pickle,time
import numpy as np
import numpy.ma as ma
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy.ma as ma
import pandas as pd

#execfile('../precip/MidpointNormalize.py')
basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/precip_clim/old-output/'
#fich1 = basedir + 'trop_depth-SST_psum_pmax_endjf.npy'
#fich2 = basedir + 'trop_depth-SST_psum_pmax_lndjf.npy'
fich1 = basedir + 'trop_depth_psum_pmax_endjf.npy'
fich2 = basedir + 'trop_depth_psum_pmax_lndjf.npy'

fs = 13
pEN = np.load(fich1); pLN = np.load(fich2)
    
# construct the joint distribution of the system depth and precip volume
mm = 6; nn = 8; small = 200; big = 250; bigpmax = 20
hEN3, xx3, yy3 = np.histogram2d(pEN[3,:],pEN[0,:],bins=(mm,nn),\
                  range=((0,bigpmax),(small,big)),normed=True)
hLN3, _, _ = np.histogram2d(pLN[3,:],pLN[0,:],bins=(mm,nn),\
                  range=((0,bigpmax),(small,big)),normed=True)
hEN3 = hEN3/np.nansum(np.nansum(hEN3,axis=1),axis=0)*100.
hLN3 = hLN3/np.nansum(np.nansum(hLN3,axis=1),axis=0)*100.

ctte = pEN[0,:]; pmaxe = pEN[3,:]
cttl = pLN[0,:]; pmaxl = pLN[3,:]
rd = np.zeros((6,10))
xval2 = np.logspace(-1,np.log10(bigpmax),11)
cttlim = np.asarray([240.,215.])
farbe = ['blue','forestgreen','red']
lbl = ['least deep','intermediate','deepest']
fs = 13; ticfont = 11

fig = plt.figure(figsize=(9,5.5))
plt.subplot2grid((1,3),(0,0))
for ii in np.arange(3):
    if ii == 0:
       jj = np.argwhere(ctte >= cttlim[ii])
       kk = np.argwhere(cttl >= cttlim[ii])
    elif ii == 1:
       jj = np.argwhere((ctte >= cttlim[ii]) & (ctte < cttlim[ii-1]))
       kk = np.argwhere((cttl >= cttlim[ii]) & (cttl < cttlim[ii-1]))
    else:
       jj = np.argwhere(ctte < cttlim[ii-1])
       kk = np.argwhere(cttl < cttlim[ii-1])
    hh, _ = np.histogram(pmaxe[jj],bins=xval2)
    hh = np.divide(hh,float(jj.shape[0]))
    gg, _ = np.histogram(pmaxl[kk],bins=xval2)
    gg = np.divide(gg,float(kk.shape[0]))

    rd[ii] = hh
    rd[ii+3] = gg

#rd[0,8] = -1*rd[0,8]
ww = np.asarray([(xval2[i]-xval2[i-1])/2 for i in np.arange(1,len(xval2))])
# where the least deep has the smallest relative difference
#mask1 = np.argwhere((abs(rd[0]) < abs(rd[1])) | (abs(rd[0]) < abs(rd[2])))
# where the intermediate has the smallest relative difference
#mask2 = np.argwhere((abs(rd[1]) < abs(rd[0])) | (abs(rd[1]) < abs(rd[2])))
ex = xval2[:-1]
plt.bar(ex,rd[0],width=ww,color=farbe[0],edgecolor='black',align='center',label=lbl[0])
plt.bar(ex,rd[1],width=ww,color=farbe[1],edgecolor='black',align='center',label=lbl[1])
plt.bar(ex,rd[2],width=ww,color=farbe[2],edgecolor='black',align='center',label=lbl[2])
plt.bar(ex[mask2[:,0]],rd[1,mask2[:,0]],width=ww[mask2[:,0]],color=farbe[1],edgecolor='black',align='center')
plt.text(0.05,0.9,'(a)',fontsize=fs+2,fontweight='bold',transform=plt.gca().transAxes)
plt.gca().set_xscale('log')
plt.ylim([-30,120])
plt.gca().tick_params(axis='both',which='both',labelsize=fs)
plt.legend(fontsize=ticfont,loc='center left')
plt.ylabel(r'Relative $\Delta p$ (EN-LN)/LN [%]',fontsize=fs)
plt.xlabel(r'$\dot{P}$ [mm h$^{-1}$]',fontsize=fs)

#plt.subplot2grid((2,2),(1,0))
#ww = np.asarray([(xval2[i]-xval2[i-1])/2 for i in np.arange(1,len(xval2))])
## where the least deep has the smallest relative difference
#mask1 = np.argwhere((abs(rd[3]) < abs(rd[4])) | (abs(rd[3]) < abs(rd[5])))
## where the intermediate has the smallest relative difference
#mask2 = np.argwhere((abs(rd[4]) < abs(rd[3])) | (abs(rd[4]) < abs(rd[5])))
#ex = xval2[:-1]
#plt.bar(ex,rd[3],width=ww,color=farbe[0],edgecolor='black',align='center')
#plt.bar(ex,rd[4],width=ww,color=farbe[1],edgecolor='black',align='center')
#plt.bar(ex,rd[5],width=ww,color=farbe[2],edgecolor='black',align='center')
#plt.bar(ex[mask2[:,0]],rd[4,mask2[:,0]],width=ww[mask2[:,0]],color=farbe[1],edgecolor='black',align='center')
#plt.text(0.05,0.9,'(b)',fontsize=fs+2,fontweight='bold',transform=plt.gca().transAxes)
#plt.xscale('log')
#plt.ylim([-2,2])
#plt.gca().tick_params(axis='both',which='both',labelsize=fs)
#plt.xlabel(r'$\dot{P}$ [mm h$^{-1}$]',fontsize=fs)
#plt.ylabel(r'Absolute $\Delta p$ (EN-LN) [%]',fontsize=fs)

basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/props_clim/'
half = 180
d1 = np.load(basedir + 'd1_map.npy')
d1 = np.concatenate((d1[half:-1,:],d1[:half,:]),axis=0)*2
d2 = np.load(basedir + 'd2_map.npy')
d2 = np.concatenate((d2[half:-1,:],d2[:half,:]),axis=0)
d3 = np.load(basedir + 'd3_map.npy')
d3 = np.concatenate((d3[half:-1,:],d3[:half,:]),axis=0)
field = np.stack((d1,d2,d3),axis=0)
field = ma.masked_where(np.isnan(field),field)

lowlat = -30; hilat = 30; lowlon = 0; hilon = 360
coul = [cm.Blues,cm.Greens,cm.Reds]
cguys = [[0,50],[0,150],[0,50]]
cgals = [np.linspace(cguys[0][0],cguys[0][1],5),\
         np.linspace(cguys[1][0],cguys[1][1],5),\
         np.linspace(cguys[2][0],cguys[2][1],5)]

titre = ['(b) Least deep systems','(c) Intermediate systems',\
         '(d) Deepest systems']
nx = d1.shape[0]; ny = d1.shape[1]
paral = np.asarray([-30,-15,15,30])
merid = np.asarray([-108,-36,36,108,180])
for i in np.arange(3):
    ax = plt.subplot2grid((3,2),(i,1))
    ax.set_title(titre[i],fontsize=fs)
    carte = Basemap(llcrnrlon=lowlon,llcrnrlat=lowlat,urcrnrlon=hilon,urcrnrlat=hilat,\
         resolution='i',projection='merc',lon_0=0.,lat_0=0.)
    lons, lats = carte.makegrid(nx, ny)
    x, y = carte(lons, lats)
    carte.drawmeridians(merid.astype(int),color='grey',linewidth=0.65,labels=[0,0,0,1],dashes=[1,1],fontsize=ticfont)
    carte.drawparallels(paral.astype(int),color='grey',linewidth=0.65,labels=[1,0,0,0],dashes=[1,1],fontsize=ticfont)
    carte.drawcoastlines(linewidth=0.75)
    cpc1 = carte.pcolor(x,y,np.flipud(np.transpose(field[i])),cmap=coul[i])   # shading='flat',
    cbar = carte.colorbar(cpc1,location='bottom',pad='40%',size='10%',ticks=cgals[i])
    if i == 2:
       cbar.set_label('# per year',fontsize=ticfont)
    cpc1.set_clim(cguys[i][0],cguys[i][1])
    cbar.ax.tick_params(labelsize=ticfont)
    ax.tick_params(axis='both',which='major',labelsize=5.5)

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0.2)
fig.savefig('./figures/fig1-precip-hist-map.pdf',bbox_inches='tight')
plt.show()

