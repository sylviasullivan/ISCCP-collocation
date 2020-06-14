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

execfile('MidpointNormalize.py')
basedir = '/work/bb1018/b380873/MCS_clim/ausgabe/props_clim/' #2deg-CG/'
fs = 13
ticfont = 11
half = 36 #half = 90
suffix = '.npy'   #'CG.npy'
d1EN = np.load(basedir + 'd1_mapEN' + suffix)
d1LN = np.load(basedir + 'd1_mapLN' + suffix)
d2EN = np.load(basedir + 'd2_mapEN' + suffix)
d2LN = np.load(basedir + 'd2_mapLN' + suffix)
d3EN = np.load(basedir + 'd3_mapEN' + suffix)
d3LN = np.load(basedir + 'd3_mapLN' + suffix)

suffix = '_JJA.npy'   #'.npy'
d1ENjja = np.load(basedir + 'd1_mapEN' + suffix)
d1LNjja = np.load(basedir + 'd1_mapLN' + suffix)
d2ENjja = np.load(basedir + 'd2_mapEN' + suffix)
d2LNjja = np.load(basedir + 'd2_mapLN' + suffix)
d3ENjja = np.load(basedir + 'd3_mapEN' + suffix)
d3LNjja = np.load(basedir + 'd3_mapLN' + suffix)

en = 41.
ln = 46.
enln = ln/en
djf = 41.
jja = 19.
jjaf = 41./19.
# Absolute difference only correcting for different number of months
#field = np.stack((d1EN*enln-d1LN,d2EN*enln-d2LN,d3EN*enln-d3LN,\
#             d1EN-d1ENjja*jjaf,d2EN-d2ENjja*jjaf,d3LN-d3LNjja*jjaf),axis=0)

# Relative differences
field = np.stack(((d1EN*enln-d1LN)/d1LN,(d2EN*enln-d2LN)/d2LN,\
              (d3EN*enln-d3LN)/d3LN,(d1EN-d1ENjja*jjaf)/(d1ENjja*jjaf),\
              (d2EN-d2ENjja*jjaf)/(d2ENjja*jjaf),(d3LN-d3LNjja*jjaf)/(d3LNjja*jjaf)),axis=0)

# Differences normalizing for month number
#field = np.stack((d1EN/en-d1LN/ln,d2EN/en-d2LN/ln,d3EN/en-d3LN/ln,\
#                  d1EN/djf-d1ENjja/jja,d2EN/djf-d2ENjja/jja,d3EN/djf-d3ENjja/jja),axis=0)
#field = field*10
field = ma.masked_where(np.isnan(field),field)

lowlat = -30; hilat = 30; lowlon = 0; hilon = 360
#coul = [cm.Blues,cm.Greens,cm.Reds,cm.Blues,cm.Greens,cm.Reds]
cguys = [[-5,5],[-15,15],[-5,5],[-5,5],[-15,15],[-5,5]]
cguys = [[-1,1],[-1,2],[-1,1],[-1,1],[-2,2],[-1,1]]
#cguys = [[-10,10],[-50,50],[-10,10],[-10,10],[-50,50],[-10,10]]
cgals = [np.linspace(cguys[0][0],cguys[0][1],5),\
         np.linspace(cguys[1][0],cguys[1][1],5),\
         np.linspace(cguys[2][0],cguys[2][1],5),\
         np.linspace(cguys[3][0],cguys[3][1],5),\
         np.linspace(cguys[4][0],cguys[4][1],5),\
         np.linspace(cguys[5][0],cguys[5][1],5)]

titre = ['(a) EN-LN: Least deep systems','(b) EN-LN: Intermediate systems',\
         '(c) EN-LN: Deepest systems','(d) DJF-JJA: Least deep systems',\
         '(e) DJF-JJA: Intermediate systems','(f) DJF-JJA: Deepest systems']
nx = d1EN.shape[0]; ny = d1EN.shape[1]
paral = np.asarray([-30,-15,15,30])
merid = np.asarray([-108,-36,36,108,180])

c = 0
fig = plt.figure(figsize=(9.2,5.5))
for j in np.arange(2):
    for i in np.arange(3):
        ax = plt.subplot2grid((3,2),(i,j))
        ax.set_title(titre[c],fontsize=fs)
        carte = Basemap(llcrnrlon=lowlon,llcrnrlat=lowlat,urcrnrlon=hilon,urcrnrlat=hilat,\
            resolution='i',projection='merc',lon_0=0.,lat_0=0.)
        lons, lats = carte.makegrid(nx, ny)
        x, y = carte(lons, lats)
        carte.drawmeridians(merid.astype(int),color='grey',linewidth=0.65,labels=[0,0,0,1],dashes=[1,1],fontsize=ticfont)
        carte.drawparallels(paral.astype(int),color='grey',linewidth=0.65,labels=[1,0,0,0],dashes=[1,1],fontsize=ticfont)
        carte.drawcoastlines(linewidth=0.75)
        cpc1 = carte.pcolor(x,y,(np.transpose(field[c])),cmap=cm.coolwarm,\
                  norm=MidpointNormalize(midpoint=0,vmin=cguys[c][0],vmax=cguys[c][1]))
        cbar = carte.colorbar(cpc1,location='bottom',pad='40%',size='10%',ticks=cgals[c])
        if i == 2:
           cbar.set_label('relative diff. count per month',fontsize=ticfont)
        cpc1.set_clim(cguys[c][0],cguys[c][1])
        cbar.ax.tick_params(labelsize=ticfont)
        ax.tick_params(axis='both',which='major',labelsize=5.5)
        c += 1

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0.2)
fig.savefig('./figures/fig-ENLN-maps.pdf',bbox_inches='tight')
plt.show()

