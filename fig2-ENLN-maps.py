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
#import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import numpy.ma as ma
import pandas as pd

basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/props_clim/'
countsEN = np.load(basedir + 'ENLN-cores-CG/NEN_counts.npy_DJF')
countsLN = np.load(basedir + 'ENLN-cores-CG/NLN_counts.npy_DJF')

fs = 13
ticfont = 11
half = 180
d1EN = np.load(basedir + 'd1_mapENCG.npy')
#d1EN = np.concatenate((d1EN[half:-1,:],d1EN[:half,:]),axis=0)*2
d1LN = np.load(basedir + 'd1_mapLNCG.npy')
#d1LN = np.concatenate((d1LN[half:-1,:],d1LN[:half,:]),axis=0)*2
d2EN = np.load(basedir + 'd2_mapENCG.npy')
#d2EN = np.concatenate((d2EN[half:-1,:],d2EN[:half,:]),axis=0)
d2LN = np.load(basedir + 'd2_mapLNCG.npy')
#d2LN = np.concatenate((d2LN[half:-1,:],d2LN[:half,:]),axis=0)
d3EN = np.load(basedir + 'd3_mapENCG.npy')
#d3EN = np.concatenate((d3EN[half:-1,:],d3EN[:half,:]),axis=0)
d3LN = np.load(basedir + 'd3_mapLNCG.npy')
#d3LN = np.concatenate((d3LN[half:-1,:],d3LN[:half,:]),axis=0)
#field = np.stack((d1EN,d2EN,d3EN,d1LN,d2LN,d3LN),axis=0)
field = np.stack((countsEN,countsLN),axis=0)
field = ma.masked_where(np.isnan(field),field)

lowlat = -30; hilat = 30; lowlon = -180; hilon = 180
coul = [cm.Blues,cm.Greens,cm.Reds,cm.Blues,cm.Greens,cm.Reds]
#coul = [cm.coolwarm,cm.coolwarm,cm.coolwarm]
#cguys = [[0,5],[0,15],[0,5],[0,5],[0,15],[0,5]]
cguys = [[0,50],[0,150],[0,50],[0,50],[0,150],[0,5]]
#cguys = [[-10,10],[-50,50],[-10,10],[-10,10],[-50,50],[-10,10]]
cgals = [np.linspace(cguys[0][0],cguys[0][1],5),\
         np.linspace(cguys[1][0],cguys[1][1],5),\
         np.linspace(cguys[2][0],cguys[2][1],5),\
         np.linspace(cguys[3][0],cguys[3][1],5),\
         np.linspace(cguys[4][0],cguys[4][1],5),\
         np.linspace(cguys[5][0],cguys[5][1],5)]

#titre = ['(a) EN: Least deep systems','(b) EN: Intermediate systems',\
#         '(c) EN: Deepest systems','(d) LN: Least deep systems',\
#         '(e) LN: Intermediate systems','(f) LN: Deepest systems']
titre = ['(a) El Nino','(b) La Nina']
#nx = d1EN.shape[0]; ny = d1EN.shape[1]
nx = countsEN.shape[0]; ny = countsEN.shape[1]
paral = np.asarray([-30,-15,15,30])
merid = np.asarray([-108,-36,36,108,180])

c = 0
fig = plt.figure(figsize=(9.2,5.5))
for j in np.arange(1):
    for i in np.arange(2):
        ax = plt.subplot2grid((3,2),(i,j))
        ax.set_title(titre[c],fontsize=fs)
        carte = Basemap(llcrnrlon=lowlon,llcrnrlat=lowlat,urcrnrlon=hilon,urcrnrlat=hilat,\
            resolution='i',projection='merc',lon_0=0.,lat_0=0.)
        lons, lats = carte.makegrid(nx, ny)
        x, y = carte(lons, lats)
        carte.drawmeridians(merid.astype(int),color='grey',linewidth=0.65,labels=[0,0,0,1],dashes=[1,1],fontsize=ticfont)
        carte.drawparallels(paral.astype(int),color='grey',linewidth=0.65,labels=[1,0,0,0],dashes=[1,1],fontsize=ticfont)
        carte.drawcoastlines(linewidth=0.75)
        #field2 = (field[c]/41. - field[c+3]/46.)/(field[c+3]/46.)*100.  # relative difference
        cpc1 = carte.pcolor(x,y,np.transpose(field[c])/25.,cmap=coul[c])   # shading='flat',
        cbar = carte.colorbar(cpc1,location='bottom',pad='40%',size='10%',ticks=cgals[c])
        if i == 2:
           cbar.set_label('count per year',fontsize=ticfont)
        cpc1.set_clim(cguys[c][0],cguys[c][1])
        cbar.ax.tick_params(labelsize=ticfont)
        ax.tick_params(axis='both',which='major',labelsize=5.5)
        c += 1

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0.2)
#fig.savefig('./figures/fig1-precip-hist-map.pdf',bbox_inches='tight')
plt.show()

