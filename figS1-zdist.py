#!/usr/bin/env mcsplot
# PLOT THE DISTRIBUTIONS OF ALTITUDES ASSOCIATED WITH THE TEMPERATURES
import sys,pickle,time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

z65EN = np.load('./zdist/z240EN.npy')
z65LN = np.load('./zdist/z240LN.npy')
z85EN = np.load('./zdist/z215EN.npy')
z85LN = np.load('./zdist/z215LN.npy')

print 'Max, min of z65EN: ' + str(np.nanmax(z65EN)), str(np.nanmin(z65EN))   
print 'Max, min of z65LN: ' + str(np.nanmax(z65LN)), str(np.nanmin(z65LN))
print 'Max, min of z85EN: ' + str(np.nanmax(z85EN)), str(np.nanmin(z85EN))
print 'Max, min of z85LN: ' + str(np.nanmax(z85LN)), str(np.nanmin(z85LN))
print 'Median of z65EN, z85EN: ' + str(np.nanmedian(z65EN)), str(np.nanmedian(z85EN))
print 'Median of z65LN, z85LN: ' + str(np.nanmedian(z65LN)), str(np.nanmedian(z85LN))
print 'Mean of z65EN, z85EN: ' + str(np.nanmean(z65EN)), str(np.nanmean(z85EN))
print 'Mean of z65LN, z85LN: ' + str(np.nanmean(z65LN)), str(np.nanmean(z85LN))
print 'St. dev. of z65EN, z85EN: ' + str(np.nanstd(z65EN)), str(np.nanstd(z85EN))
print 'St. dev. of z65LN, z85LN: ' + str(np.nanstd(z65LN)), str(np.nanstd(z85LN))
 
# construct the histograms of El Nino and La Nina associated altitudes
xval = np.linspace(np.nanmin(z65EN),np.nanmax(z65EN),20)
xval2 = np.linspace(np.nanmin(z85EN),np.nanmax(z85EN),20)
lbl = ['EN z(240 K)','LN z(240 K)','EN z(215 K)','LN z(215 K)']
farbe = ['seagreen','lightgreen','red','pink']
axfont = 10
ticfont = 9
#uu = 120
#ll = -30

fig = plt.figure(figsize=(7,4.3))
wgts1 = np.ones_like(z65EN)/float(len(z65EN))*100.
h1,_ = np.histogram(z65EN,bins=xval,weights=wgts1)
wgts2 = np.ones_like(z65LN)/float(len(z65LN))*100.
h2,_ = np.histogram(z65LN,bins=xval,weights=wgts2)
mask1 = np.argwhere(h1 < h2)
ww = np.asarray([(xval[i]-xval[i-1]) for i in np.arange(1,len(xval))])
plt.bar(xval[:-1],h1,width=ww,color=farbe[0],edgecolor='black',label=lbl[0])
plt.bar(xval[:-1],h2,width=ww,color=farbe[1],edgecolor='black',label=lbl[1])
plt.bar(xval[mask1[:,0]],h1[mask1[:,0]],width=ww[mask1[:,0]],color=farbe[0],edgecolor='black')

wgts1 = np.ones_like(z85EN)/float(len(z85EN))*100.
h1,_ = np.histogram(z85EN,bins=xval2,weights=wgts1)
wgts2 = np.ones_like(z85LN)/float(len(z85LN))*100.
h2,_ = np.histogram(z85LN,bins=xval2,weights=wgts2)
mask1 = np.argwhere(h1 < h2)
ww = np.asarray([(xval2[i]-xval2[i-1]) for i in np.arange(1,len(xval2))])
plt.bar(xval2[:-1],h1,width=ww,color=farbe[2],edgecolor='black',label=lbl[2])
plt.bar(xval2[:-1],h2,width=ww,color=farbe[3],edgecolor='black',label=lbl[3])
plt.bar(xval2[mask1[:,0]],h1[mask1[:,0]],width=ww[mask1[:,0]],color=farbe[2],edgecolor='black')

plt.ylabel('Probability [%]')
plt.xlabel('Depth [km]')
plt.legend(loc='center right')

#basedir = '/rigel/home/scs2229/top-secret/MCS_clim/scripts/figs/article9-reviews/'
fig.savefig('figures/figS1-zdist.pdf',bbox_inches='tight')
plt.show()
sys.exit()

