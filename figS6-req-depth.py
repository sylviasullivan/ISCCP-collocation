import numpy as np
import matplotlib.pyplot as plt
import sys

basedir = '/work/bb1018/b380873/MCS_clim/data/ISCCP/'
data = np.load(basedir + 'CT-allDataLocalTimes.npy')
writedata = False
if writedata == True:
   LNDJFtups = [(1983,10),(1983,11),(1983,12),(1984,11),(1984,12),(1985,1),(1985,2),\
             (1985,3),(1985,4),(1988,10),(1988,11),(1988,12),(1989,1),(1989,2),\
             (1989,3),(1989,4),(1995,10),(1995,11),(1995,12),(1996,1),(1996,2),\
             (1998,10),(1998,11),(1998,12),(1999,1),(1999,2),(1999,3),(1999,4),\
             (1999,10),(1999,11),(1999,12),(2000,1),(2000,2),(2000,3),(2000,4),\
             (2000,10),(2000,11),(2000,12),(2001,1),(2005,12),(2006,1),(2006,2),\
             (2006,2),(2007,10),(2007,11),(2007,4)]
   ENDJFtups = [(1986,10),(1986,11),(1986,12),(1987,1),(1987,2),(1987,3),(1987,4),\
             (1984,10),(1987,11),(1987,12),(1988,1),(1991,10),(1991,11),(1991,12),\
             (1992,1),(1992,2),(1992,3),(1992,4),(1994,10),(1994,11),(1994,12),\
             (1995,1),(1995,2),(1997,10),(1997,11),(1997,12),(1998,1),(1998,2),\
             (1998,3),(1998,4),(2002,10),(2002,11),(2002,12),(2003,1),(2004,10),\
             (2004,11),(2004,12),(2005,1),(2006,10),(2006,11),(2006,12)]

   dataEN = []; dataLN = []
   for row in data:
       if((row[7],row[8]) in LNDJFtups):
          dataLN.append(row)
       if((row[7],row[8]) in ENDJFtups):
          dataEN.append(row)
   np.save(basedir + 'CT-allDataLocalTimes_EN.npy',np.asarray(dataEN))
   np.save(basedir + 'CT-allDataLocalTimes_LN.npy',np.asarray(dataLN))
else:
   dataEN = np.load(basedir + 'CT-allDataLocalTimes_EN.npy')
   dataLN = np.load(basedir + 'CT-allDataLocalTimes_LN.npy')

indx = np.random.random_integers(0,dataEN.shape[0],1000)
reqEN = dataEN[indx,11]   # the 11th element is the equivalent radius, 2
depEN = dataEN[indx,25]   # the 25th element is the cloud top temp, 3
indx = np.random.random_integers(0,dataLN.shape[0],1000)
reqLN = dataLN[indx,11]   # the 11th element is the equivalent radius
depLN = dataLN[indx,25]   # the 25th element is the cloud top temp

fs = 13
ax1 = plt.subplot2grid((2,2),(0,0))
ax1.scatter(depEN,reqEN,s=15,edgecolor='k',color='r')
ax1.set_ylabel('Equivalent radius [km]',fontsize=fs)
ax1.set_xlabel('Cloud top temperature [K]',fontsize=fs)
ax1.text(0.1,0.9,'(a)',fontweight='bold',fontsize=fs+2,transform=ax1.transAxes)

ax2 = plt.subplot2grid((2,2),(0,1))
ax2.scatter(depLN,reqLN,s=15,edgecolor='k',color='b')
ax2.set_xlabel('Cloud top temperature [K]',fontsize=fs)
ax2.text(0.1,0.9,'(b)',fontweight='bold',fontsize=fs+2,transform=ax2.transAxes)

# binning part
reqEN = dataEN[:,2]; reqLN = dataLN[:,2]   #11,2
cttEN = dataEN[:,3]; cttLN = dataLN[:,3]   #25,3
ctt = np.arange(180,244,2)
rMLN = np.zeros(len(ctt-1),)
rMEN = np.zeros(len(ctt-1),)
cttM = np.zeros(len(ctt-1),)
for i in np.arange(len(ctt)-1):
    j = np.argwhere((cttEN > ctt[i]) & (cttEN <= ctt[i+1]))
    rMEN[i] = np.nanmean(reqEN[j[:,0]])
    j = np.argwhere((cttLN > ctt[i]) & (cttLN <= ctt[i+1]))
    rMLN[i] = np.nanmean(reqLN[j[:,0]])
    cttM[i] = (ctt[i] + ctt[i+1])/2.
print(rMEN)
print(rMLN)
print(cttM)

ax = plt.subplot2grid((2,1),(1,0))
ax.plot(cttM[:-2],rMEN[:-2],color='r',linewidth=1.25)
ax.plot(cttM[:-1],rMLN[:-1],color='b',linewidth=1.25)
ax.set_xlim([180,245])
ax.set_ylabel('Equivalent radius [km]',fontsize=fs)
ax.set_xlabel('Cloud top temperature [K]',fontsize=fs)
ax.text(0.1,0.9,'(c)',fontweight='bold',fontsize=fs+2,transform=ax.transAxes)

plt.tight_layout()
#plt.savefig('./figures/figS6-req-depth.pdf')
plt.show()

