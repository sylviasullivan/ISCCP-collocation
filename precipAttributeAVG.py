#!/usr/bin/env mcsplot
# EVALUATE THE MEAN ATTRIBUTION FIELD FROM THE MONTHLY ATTRIBUTIONS
# ADJUST phase TO CHANGE EL NINO v LA NINA AND SEASON

import time,sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset,num2date
from datetime import datetime,timedelta
from functools import reduce

################# SET THESE VARIABLES BEFORE RUNNING ########################
lowlat = -30; hilat  = 30; lowlon = 0; hilon  = 360    # subdomain
phase  = 'LNJJA'                         # which phase and season?
basedir2 = '/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/precip_clim/MSWEP/'
#############################################################################

if phase == 'LNDJF':
   yrs = [1983,1983,1983,1984,1984,1985,1985,1985,1985,1988,1988,\
          1988,1989,1989,1989,1989,1995,1995,1995,1996,1996,1998,\
          1998,1998,1999,1999,1999,1999,1999,1999,1999,2000,2000,\
          2000,2000,2000,2000,2000,2001,2005,2006,2006,2007,2007,\
          2007,2008,2008,2007,2007]
   mon = ['10','11','12','11','12','01','02','03','04','10','11',\
          '12','01','02','03','04','10','11','12','01','02','10',\
          '11','12','01','02','03','04','10','11','12','01','02','03',\
          '04','10','11','12','01','12','01','02','10','11','12','01',\
          '02','03','04']
   filename = 'ATT_LNDJF'
   basedir  = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-la-nina-yrs/DJF/'
elif phase == 'ENDJF':
     yrs = [1986,1986,1986,1987,1987,1987,1987,1987,1987,1987,1988,1991,\
            1991,1991,1992,1992,1992,1992,1994,1994,1994,1995,1995,1997,\
            1997,1997,1998,1998,1998,1998,2002,2002,2002,2003,2004,2004,\
            2004,2005,2006,2006,2006]
     mon = ['10','11','12','01','02','03','04','10','11','12','01','10',\
            '11','12','01','02','03','04','10','11','12','01','02','10',\
            '11','12','01','02','03','04','10','11','12','01','10','11',\
            '12','01','10','11','12']
     filename = 'ATT_ENDJF'
     basedir  = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-el-nino-yrs/DJF/'
elif phase == 'LNJJA':
     yrs = [1985,1985,1985,1985,1988,1988,1988,1988,1989,1995,1998,1998,\
            1999,1999,1999,1999,1999,2000,2000,2000,2000,2000,2007]
     mon = ['05','06','07','08','06','07','08','09','05','09','08','09',\
            '05','06','07','08','09','05','06','07','08','09','09']
     filename = 'ATT_LNJJA'
     basedir  = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-la-nina-yrs/JJA/'
elif phase == 'ENJJA':
     yrs = [1987,1987,1987,1987,1987,1991,1991,1991,1991,1997,1997,1997,\
            1997,2002,2002,2002,2004,2004]
     mon = ['05','06','07','08','09','06','07','08','09','06','07','08',\
            '09','07','08','09','08','09']
     filename = 'ATT_ENJJA'
     basedir  = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-el-nino-yrs/JJA/'

# extract all the precipitation values collocated with MCS
# for attributions
pr = np.zeros((14,72,len(yrs)))
for ii in range(len(yrs)):
    mswep = np.load(basedir + 'ATT_' + str(yrs[ii]) + mon[ii] + '_3h_CG2.npy')
    pr[:,:,ii] = mswep

pr = np.nanmedian(pr,axis=2)
print basedir2 + filename
np.save(basedir2 + filename + '_frac.npy',pr)

sys.exit()

# for absolute precipitation values
pr = np.zeros((140,720,len(yrs)))
for ii in range(len(yrs)):
    mswep = Dataset(basedir + 'MSWEP_prec' + str(yrs[ii]) + mon[ii] + '_3h_halfdeg.nc','r+')
    pr[:,:,ii] = mswep
    temp = mswep['precipitation']
    temp = np.nanmean(temp,axis=0)
    pr[:,:,ii] = temp[110:250,:]


np.save('mean_precip_ENDJF.npy',pr)
