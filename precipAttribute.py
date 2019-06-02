#!/usr/bin/env ncplot
# CALCULATE THE MATRIX OF MCS-COLLOCATED PRECIPITATION VALUES
# AND DIVIDE IT BY THE TOTAL PRECIPITATION TO GENERATE THE MCS
# ATTRIBUTION SAVED IN A "ATT_*" TYPE FILE

# adjust yrs, mon, and basedir for (ENDJF, ENJJA, LNDJF versus LNJJA)
import time,sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset,num2date
from datetime import datetime,timedelta,date
from functools import reduce

################# SET THESE VARIABLES BEFORE RUNNING ########################
# subdomain for which you want to filter
lowlat = -30; hilat  = 30; lowlon = 0; hilon  = 360
phase = 'LNDJF'  #'ENDJF','LNDJF','ENJJA'
#############################################################################

if phase == 'LNDJF':
   # La Nina DJF periods
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
   basedir = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-la-nina-yrs/DJF/'

elif phase == 'ENDJF':
     # El Nino DJF periods
     yrs = [1986,1986,1986,1987,1987,1987,1987,1987,1987,1987,1988,1991,\
            1991,1991,1992,1992,1992,1992,1994,1994,1994,1995,1995,1997,\
            1997,1997,1998,1998,1998,1998,2002,2002,2002,2003,2004,2004,\
            2004,2005,2006,2006,2006]
     mon = ['10','11','12','01','02','03','04','10','11','12','01','10',\
            '11','12','01','02','03','04','10','11','12','01','02','10',\
            '11','12','01','02','03','04','10','11','12','01','10','11',\
            '12','01','10','11','12']
     basedir = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-el-nino-yrs/DJF/'

elif phase == 'LNJJA':
     # La Nina JJA periods
     yrs = [1985,1985,1985,1985,1988,1988,1988,1988,1989,1995,1998,1998,\
            1999,1999,1999,1999,1999,2000,2000,2000,2000,2000,2007]
     mon = ['05','06','07','08','06','07','08','09','05','09','08','09',\
            '05','06','07','08','09','05','06','07','08','09','09']
     basedir = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-la-nina-yrs/JJA/'

elif phase == 'ENJJA':
     # El Nino JJA periods
     yrs = [1987,1987,1987,1987,1987,1991,1991,1991,1991,1997,1997,1997,\
           1997,2002,2002,2002,2004,2004]
     mon = ['05','06','07','08','09','06','07','08','09','06','07','08',\
            '09','07','08','09','08','09']
     basedir = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-el-nino-yrs/JJA/'

# load CT data
datei  = "/rigel/home/scs2229/top-secret/MCS_clim/data/ISCCP/CT-allDataLocalTimes.npy"
CTdata = np.load(datei)
DD = np.zeros((CTdata.shape[0],9))
DD[:,0]  = CTdata[:,7]       # years
DD[:,1]  = CTdata[:,8]       # months
DD[:,2]  = CTdata[:,9]       # days
DD[:,3]  = CTdata[:,10]      # hours
DD[:,4]  = CTdata[:,40]      # min lat
DD[:,5]  = CTdata[:,41]      # max lat
DD[:,6]  = CTdata[:,42]      # min lon
DD[:,7]  = CTdata[:,43]      # max lon
DD[:,8]  = CTdata[:,17]      # convective fraction
# center lat = 12, center lon = 13

# filter for the subdomain here
indx1 = np.argwhere(DD[:,5] <= hilat)
indx2 = np.argwhere(DD[:,4] >= lowlat)
indx3 = np.argwhere(DD[:,7] <= hilon)
indx4 = np.argwhere(DD[:,6] >= lowlon)
indx5 = np.argwhere(DD[:,7] != 360)
indx  = reduce(np.intersect1d,(indx1,indx2,indx3,indx4,indx5))
DD2   = DD[indx,:]
del DD

# calculate and save an array of precipitation intensities collocated with MCS
for ii in range(len(yrs)):
    print(str(yrs[ii]) + mon[ii])
    cumprec = np.zeros((140,720))
    mswep = Dataset(basedir + 'MSWEP_prec' + str(yrs[ii]) + mon[ii] + '_3h_halfdeg.nc','r+')
    prec  = mswep.variables['precipitation']
    prec  = prec[:,110:251,:]                          # filter for the tropics from 35S to 35N
    prec  = np.multiply(prec,3.)                       # multiply by 3 hours to get accumulation
    zeit  = num2date(mswep.variables['time'][:],'days since 1899-12-31 00:00:00')

    # filter for the year and month now
    if mon[ii][0] == '0':
       monat = float(mon[ii][1])
    else:
       monat = float(mon[ii])
    indx1 = np.argwhere((DD2[:,0] == yrs[ii]) & (DD2[:,1] == monat))
    DD3 = DD2[indx1[:,0]]

    # ensure that the hours exist on a 3-hrly grid
    CTtimes = np.zeros((DD3.shape[0]),dtype='datetime64[h]')
    for jj,row in enumerate(DD3):
        if int(row[3])%3 != 0:
           newhour = int(3 * round(float(row[3])/3))
           if newhour > 21:
              tomorrow = date(int(row[0]),int(row[1]),int(row[2])) + timedelta(days=1)
              CTtimes[jj] = datetime(tomorrow.year,tomorrow.month,tomorrow.day,0)
           else:
              CTtimes[jj] = datetime(int(row[0]),int(row[1]),int(row[2]),newhour)
        else:
           CTtimes[jj] = datetime(int(row[0]),int(row[1]),int(row[2]),int(row[3]))
    DD3[:,0] = np.asarray([CTtimes[jj].astype(object).year for jj,_ in enumerate(CTtimes)])
    DD3[:,1] = np.asarray([CTtimes[jj].astype(object).month for jj,_ in enumerate(CTtimes)])
    del CTtimes

    # collocate MCS and MSWEP values at this year / month
    for jj in range(zeit.shape[0]):
        indx3 = np.argwhere(DD3[:,2] == zeit[jj].day)
        indx4 = np.argwhere(DD3[:,3] == zeit[jj].hour) # | \
                    # (DD2[:,3] == (zeit[jj] - timedelta(hours=3)).hour))
                    # include precip instances that proceed the MCS by 3 hrs
        indx  = reduce(np.intersect1d,(indx3,indx4))
        subset = DD3[indx,:]

        lat = np.asarray([int(round(subset[kk,4]*2+70)) for kk in range(len(indx))])
        lat = np.vstack((lat,np.asarray([int(round(subset[kk,5]*2+70)) \
                 for kk in range(len(indx))])))
        lat = np.vstack((lat,[int(round(subset[kk,8])) for kk in range(len(indx))]))

        lon = np.asarray([int(round(subset[kk,6]*2)) for kk in range(len(indx))])
        lon = np.vstack((lon,np.asarray([int(round(subset[kk,7]*2)) \
                 for kk in range(len(indx))])))
        lon[lon == 720] = 0

        pr = prec[jj,:,:]
        if(lat.shape[1] != 0 and lon.shape[1] != 0):
           for qq in range(len(indx)):
               rr = lat[0,qq]; ss = lon[0,qq]; ff = lat[2,qq]
               while rr <= lat[1,qq]:
                     ss = lon[0,qq]
                     while ss <= lon[1,qq]:
                           cumprec[rr,ss] += pr[rr,ss]*ff/100. # weight by conv fraction
                           ss += 1

                     rr += 1
               while rr <= lat[1,qq]:
                     ss = lon[0,qq]
                     while ss <= lon[1,qq]:
                           cumprec[rr,ss] += pr[rr,ss]*ff/100. # weight by conv fraction
                           ss += 1

                     rr += 1

    # total cumulative precipitation over the month, again coarse-graining
    totprec = np.nansum(prec,axis=0)   #ma.masked_less(np.nansum(prec,axis=0),0.1)
    totprec[totprec < 0] = 0

    # coarse grain the cumulative and total precipitation
    cumprecCG = np.zeros((cumprec.shape[0]/10,cumprec.shape[1]/10))
    totprecCG = np.zeros((totprec.shape[0]/10,totprec.shape[1]/10))

    lll = np.arange(0,cumprec.shape[0],10); mmm = np.arange(0,cumprec.shape[1],10)
    for kk in range(0,lll.shape[0]):
        for pp in range(0,mmm.shape[0]):
            subset2 = cumprec[kk*10:(kk+1)*10,pp*10:(pp+1)*10]
            cumprecCG[kk,pp] = sum(sum(subset2))
            subset2 = totprec[kk*10:(kk+1)*10,pp*10:(pp+1)*10]
            totprecCG[kk,pp] = sum(sum(subset2))

    # calculate the ratio of cumulative MCS-collocated to total precip
    precatt = cumprecCG/totprecCG*100.
    print(np.nanmax(precatt))
    toomuch = precatt[precatt > 100]
    print(len(toomuch))
    print(basedir + 'ATT_' + str(yrs[ii]) + mon[ii] + '_3h_CG.npy')   #  _nofrac
    np.save(basedir + 'ATT_' + str(yrs[ii]) + mon[ii] + '_3h_CG2.npy',precatt)   # _nofrac

