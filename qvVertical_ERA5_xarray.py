#!/usr/bin/env hodograph
# COLLOCATE VERTICAL PROFILES IN THE SPECIFIC HUMIDITY WITH MCS OCCURENCE. 
# jahr AND mois ARE THE YEAR AND MONTH FOR WHICH THE NETCDF FILE EIXSTS.
# delay IS THE TIME DELAY FROM MCS OCCURRENCE
# PRODUCES A LIST OF MCS-COLLOCATED VERTICAL PROFILES OF QV.

def qvVertical_ERA5_xarray(jahr,mois,delay):
    import time,sys
    import numpy as np
    import xarray as xr
    from datetime import datetime,timedelta
    from to_datetime import to_datetime

    # tropical subdomain
    lowlat  = -30; hilat  = 30; lowlon = 0; hilon  = 360
    
    # load ERA-5 temperature data
    # ttt [=] (744 times, 26 pressure levels, 241 latitudes, 1440 longitudes)
    dataDIR = '/rigel/home/scs2229/top-secret/MCS_clim/ERA5_qv.nc'
    DS = xr.open_dataset(dataDIR)
    qvqv = DS.q
    zeit = DS.time
    print 'File start and end times: ' + str(zeit[0].values) + ' , ' + str(zeit[-1].values)

    # convert the ncfile times to datetime objects for comparison
    zeit2 = np.zeros((qvqv.shape[0]),dtype='datetime64[h]')
    for ii in np.arange(len(zeit2)):
        zeit2[ii] = to_datetime(zeit[ii].values)
    del zeit

    # load CT data
    datei  = '/rigel/home/scs2229/top-secret/MCS_clim/data/ISCCP/CT-allDataLocalTimes.npy'
    CTdata = np.load(datei)
    DD = np.zeros((CTdata.shape[0],16))
    DD[:,0] = CTdata[:,7]      # years
    DD[:,1] = CTdata[:,8]      # months
    DD[:,2] = CTdata[:,9]      # days
    DD[:,3] = CTdata[:,10]     # hours
    DD[:,4] = CTdata[:,21]     # 21 = core lat, 12 = sys lat
    DD[:,5] = CTdata[:,22]     # 22 = core lon, 13 = sys lon
    DD[:,6] = CTdata[:,37]     # land/water flag
    DD[:,7] = CTdata[:,2]      # maximum radius
    DD[:,8] = CTdata[:,3]      # minimum temp
    DD[:,9] = CTdata[:,5]      # lifetime 
    DD[:,10] = CTdata[:,11]     # CS radius
    DD[:,11] = CTdata[:,12]     # center lat
    DD[:,12] = CTdata[:,13]     # center lon
    DD[:,13] = CTdata[:,17]     # conv fraction
    DD[:,14] = CTdata[:,18]     # number of cores
    DD[:,15] = CTdata[:,24]     # CS temp
    del CTdata

    # filter out instance where the core position was not recorded (~19% of values)
    tossout = []
    cc = 1
    for ii in range(DD.shape[0]):
        if DD[ii,5] < 0:
           tossout.append(ii)
        elif int(DD[ii,5]) == 360:
           tossout.append(ii)
    DD = np.delete(DD,tossout,axis=0)

    # locate based upon cores
    lala = DD[:,4]; lolo = DD[:,5]
    # locate based on MCS center
    #lala = DD[:,11]; lolo = DD[:,12]

    indx = np.argwhere((lala < hilat) & (lala >= lowlat) & (lolo < hilon) & \
            (lolo >= lowlon) & (DD[:,0] == jahr) & (DD[:,1] == mois))
    DD = DD[indx[:,0],:]
 
    # transform the year / month / day / hr format to datetime
    # 32635 instances in 1983 for example
    CTtimes = np.zeros((DD.shape[0]),dtype='datetime64[h]')
    CTtimes = np.asarray([datetime(int(DD[ii,0]),int(DD[ii,1]),int(DD[ii,2]),int(DD[ii,3]),0) for ii in range(DD.shape[0])])
    print 'Starting with ' + str(CTtimes.shape[0]) + ' MCS occurrences in the month of interest.'

    # append each collocated 26-level profile to qvARR
    qvARR = []
    for ii in np.arange(DD.shape[0]):
        obj = CTtimes[ii]

        # identify the vertical profiles preceding by 6 hrs MCS occurrence
        if((obj - timedelta(hours=delay)).month != mois):
           continue
        else:
           colloc = np.argwhere((zeit2 == obj - timedelta(hours=delay)))
           # collocating with system here (11,12) versus core (4,5)
           lat25 = round((hilat - DD[ii,4])*4)
           # ERA longitudes are stored from -180 to 180
           if DD[ii,5] > 180:
              DDlon = DD[ii,5] - 360
           else:
              DDlon = DD[ii,5]
           # find the nearest longitude to 0.25
           lon25 = round((180 + DDlon)*4)
           if lon25 == 1440:
              lon25 = 0
           if(len(colloc) != 0):
              qvprof = qvqv[colloc[0][0],:,int(lat25),int(lon25)]
              qvARR.append(np.concatenate((qvprof,DD[ii,:]),axis=0))
           # ERA-5 profiles should be available every hour so that colloc always exists
           # but verify this
           else:
              print 'No collocation made for ' + str(obj - timedelta(hours=delay))
    print 'Final size of the specific humidity profile list is ' + str(np.asarray(qvARR).shape)
    return np.asarray(qvARR)

