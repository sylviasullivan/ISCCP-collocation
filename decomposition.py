#!/usr/bin/env mcsplot
# DECOMPOSE CHANGES IN PRECIPITATION (totalchange) INTO THOSE DUE TO MCS
# FREQUENCY / PRECIPITATION OCCURRENCE (freqchange) AND THOSE DUE TO PRECIP
# QUANTITTY (precipchange). resid REPRESENTS THE RESIDUAL OF totalchange -
# freqchange + precipchange. THIS ANALYSIS IS AS IN JACKSON ET AL. 2017
# Increases in tropical rainfall driven by change in frequency of organized
# deep convection
# PRECIPITATION CHANGES ARE TAKEN EITHER BETWEEN EARLY (1983-1995) AND LATE 
# (1996-2008) (earlylate = True) OR BETWEEN WARM-COLD PHASE (ENLN = True).
import time,sys,pickle,warnings
import numpy as np
import matplotlib
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from netCDF4 import Dataset,num2date

################# SET THESE VARIABLES BEFORE RUNNING ########################
totalchange = True
freqchange = False
precipchange = False
resid = False
earlylate = False
ENLN = True
lowlat = -33; hilat  = 33; lowlon = -180; hilon  = 180
var = 4
#############################################################################

# load the ERA total and convective precipitation values into tp and cp
basedir = '../../data/ERA/decomposition/'
dd0 = np.load('../../data/ISCCP/CT-allDataFilteredEdges.npy')
dd1 = Dataset(basedir + 'ERAint_monthly_mean_tp1.nc','r+')
dd2 = Dataset(basedir + 'ERAint_monthly_mean_tp2.nc','r+')
dd3 = Dataset(basedir + 'ERAint_monthly_mean_tp3.nc','r+')
tp = dd1.variables['tp'][:] +  dd2.variables['tp'][:] + dd3.variables['tp'][:]
tp = np.concatenate((tp[90:,:],tp[:90,:]))

dd1 = Dataset(basedir + 'ERAint_monthly_mean_cp1.nc','r+')
dd2 = Dataset(basedir + 'ERAint_monthly_mean_cp2.nc','r+')
dd3 = Dataset(basedir + 'ERAint_monthly_mean_cp3.nc','r+')
cp = dd1.variables['cp'][:] + dd2.variables['cp'][:] + dd3.variables['cp'][:]
cp = np.concatenate((cp[90:,:],cp[:90,:]))
if earlylate == True:
   tim = dd1.variables['time'][:]
elif ENLN == True:
   tim = num2date(dd1.variables['time'][:],'hours since 1900-01-01 00:00:0.0')
else:
   print('Neither earlylate nor ENLN are designated for taking precip differences.')

print('Size of precipitation data: ' + str(cp.shape))
del dd1
del dd2
del dd3

if earlylate == True:
   # separate the precip values 1983-1995 and 1996-2008
   print('Precipitation changes taken between 1983-1995 and 1996-2008.')
   indx = np.argwhere(tim > 841512)
   indx = np.asscalar(indx[0])
   cp1 = np.average(cp[3:indx,29:-28,:],axis=0)
   cp2 = np.average(cp[indx:,29:-28,:],axis=0)
   tp1 = np.average(tp[3:indx,29:-28,:],axis=0)
   tp2 = np.average(tp[indx:,29:-28,:],axis=0)
   print('Precipitation data - first half: ' + str(tp1.shape))
   print('Precipitation data - second half: ' + str(tp2.shape))
elif ENLN == True:
   # separate the precip values El Nino versus La Nina
   # EN = (1983, 1987, 1988, 1992) versus (1995, 1998, 2003, 2007)
   indx11 = [tim[ii].year == 1989 for ii in range(len(tim))]
   indx21 = [tim[ii].year == 1999 or tim[ii].year == 2000 or tim[ii].year == 2008 \
                  for ii in range(len(tim))]
   indx12 = [tim[ii].month == 1 or tim[ii].month == 2 for ii in range(len(tim))]
   indx1  = [indx11[ii] and indx12[ii] for ii in range(len(indx11))]
   indx2  = [indx21[ii] and indx12[ii] for ii in range(len(indx21))]
   cp1 = np.average(cp[indx1,29:-28,:],axis=0)
   cp2 = np.average(cp[indx2,29:-28,:],axis=0)
   tp1 = np.average(tp[indx1,29:-28,:],axis=0)
   tp2 = np.average(tp[indx2,29:-28,:],axis=0)
   print('Precipitation data - first half: ' + str(tp1.shape))
   print('Precipitation data - second half: ' + str(tp2.shape))

DD = np.zeros((dd0.shape[0],var))
DD[:,0] = dd0[:,12]         # lat
DD[:,1] = dd0[:,13]         # lon
if (lowlon < 0):
   indx = np.argwhere(DD[:,1] > 180)
   DD[indx,1] = DD[indx,1] - 360
DD[:,2] = dd0[:,8]          # month
DD[:,3] = dd0[:,7]          # year
print('Size of MCS data: '+ str(DD.shape))
del dd0

if earlylate == True:
   # filter those MCS 1983-1995 and 1996-2008
   indx1 = []; indx2 = []
   for jj in range(0,DD.shape[0]):
       if (DD[jj,3] < 1996):
          indx1.append(jj)
       else:
          indx2.append(jj) 
   DD1 = DD[indx1,:]; DD2 = DD[indx2,:]
   DD1 = np.reshape(DD1,[len(indx1),var])
   DD2 = np.reshape(DD2,[len(indx2),var])
   print('MCS data - first half: ' + str(DD1.shape))
   print('MCS data - second half: ' + str(DD2.shape))
   del DD
elif ENLN == True:
   # La Nina = (1989) versus (1999, 2000,2008)
   # El Nino = (1983, 1987, 1988, 1992) versus (1995, 1998, 2003, 2007)
   indx11 = [DD[jj,3] == 1989 for jj in range(DD.shape[0])]
   indx21 = [DD[jj,3] == 1999 or DD[jj,3] == 2000 or DD[jj,3] == 2008 for jj in range(DD.shape[0])]
   indx12 = [DD[jj,2] == 1 or DD[jj,2] == 2 for jj in range(DD.shape[0])]
   indx1  = [indx11[ii] and indx12[ii] for ii in range(len(indx11))]
   indx2  = [indx21[ii] and indx12[ii] for ii in range(len(indx21))]
   DD1 = DD[indx1,:]; DD2 = DD[indx2,:]
   print('MCS data - first half: ' + str(DD1.shape))
   print('MCS data - second half: ' + str(DD2.shape))
   del DD

# density maps with dimensions of the tropical region 
dic1 = np.zeros((int(hilon-lowlon)+1,int(hilat-lowlat)+1))
dic2 = np.zeros((int(hilon-lowlon)+1,int(hilat-lowlat)+1))

for jj in DD1:
    if (jj[0] > lowlat and jj[0] < hilat and jj[1] >= lowlon and jj[1] <= hilon):  
        dic1[int(jj[1])+abs(lowlon), int(jj[0])+abs(lowlat)] += 1

for jj in DD2:
    if (jj[0] > lowlat and jj[0] < hilat and jj[1] >= lowlon and jj[1] <= hilon):
        dic2[int(jj[1])+abs(lowlon), int(jj[0])+abs(lowlat)] += 1

# Coarse-grain the data to a 2-degree grid
dic1cg = np.zeros((dic1.shape[0]/2 + 1,dic1.shape[1]/2 + 1))
dic2cg = np.zeros((dic2.shape[0]/2 + 1,dic2.shape[1]/2 + 1))
lat = np.arange(0,dic1.shape[0],2); lon = np.arange(0,dic1.shape[1],2)
for ii in range(0,lat.shape[0]):
    for jj in range(0,lon.shape[0]):
        subset1 = dic1[ii*2:(ii+1)*2,jj*2:(jj+1)*2]
        subset2 = dic2[ii*2:(ii+1)*2,jj*2:(jj+1)*2]
        dic1cg[ii,jj] = sum(sum(subset1))
        dic2cg[ii,jj] = sum(sum(subset2))
dic1cg = dic1cg/147.        # calculate a monthly avg from a cumulative frequency
dic2cg = dic2cg/147.

cpmean = np.transpose(np.average(cp,axis=0))
cpmean = cpmean[:,29:-28]
dic1cg = dic1cg[:180,:]; dic2cg = dic2cg[:180,:]
freqmean = (dic1cg + dic2cg)/2.               # summed frequency over 25 years --> mon avg

deltaFMCS = np.transpose(cpmean*1000.*0.4*(dic2cg - dic1cg))       # 1000/24 to convert m to mm (day-1)
deltaP = (tp2 - tp1)*1000.
deltaCP = np.transpose(freqmean)*(cp2 - cp1)*1000.*0.4
residual = deltaP - deltaFMCS - deltaCP

######################### Precipitation decomposition ############################
fig = plt.figure()
nx = dic1cg.shape[0]; ny = dic2cg.shape[1]
paral=np.linspace(lowlat,hilat,4); merid=np.linspace(lowlon,hilon,6)
axes = plt.subplot2grid((1,1),(0,0))

if totalchange == True:
   titre = 'Total change in precipitation'
   field = np.concatenate((deltaP[90:,:],deltaP[:90,:]))
   cbl   = ''
   cl    = [-1.6,1.6]
elif freqchange == True:
   titre = 'Mean MCS precip x change in MCS frequency'
   field = np.concatenate((deltaFMCS[90:,:],deltaFMCS[:90,:]))
   cbl   = '[mm day-1]'
   cl    = [-1.6,1.6]
elif precipchange == True:
   titre = 'Mean MCS frequency x change in MCS precip'
   field = np.concatenate((deltaCP[90:,:],deltaCP[:90,:]))
   cbl   = '[mm day-1]'
   cl    = [-1.6,1.6]
elif resid == True:
   titre = 'Residual'
   field = np.concatenate((residual[90:,:],residual[:90,:]))
   cbl   = '[mm day-1]'
   cl    = [-1.6,1.6]
else:
   print('None of the visualization booleans were specified.')
   sys.exit()

axes.set_title(titre)
carte = Basemap(llcrnrlon=lowlon,urcrnrlon=hilon,llcrnrlat=lowlat,urcrnrlat=hilat,\
                   resolution='i',projection='merc',lon_0=180.,ax=axes)
lons, lats = carte.makegrid(nx, ny)
x, y = carte(lons, lats)
carte.drawparallels(paral.astype(int),labels=[1,1,0,0],dashes=[1,1])
carte.drawmeridians(merid.astype(int),labels=[0,0,0,1],dashes=[1,1])
carte.drawcoastlines(linewidth=1)
cpc = carte.pcolor(x,y,field,shading='flat',cmap=cm.RdBu_r)
cbar = carte.colorbar(cpc,location='bottom',pad='20%',size='15%')
cbar.set_label(cbl)
cpc.set_clim(cl[0],cl[1])

plt.show()
