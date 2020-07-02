#!/usr/bin/env hodograph
# COLLOCATE VERTICAL PROFILES IN SPECIFIC HUMIDITY WITH MCS OCCURENCE. 
# anfang AND ende SHOULD BE DATETIME OBJECTS THAT ARE THE BEGINNING AND END 
# OF THE YEAR OF INTEREST. PRODUCES A LIST OF COLLOCATED SPECIFIC HUMIDITY 
# PROFILES FOR THE GIVEN YEAR.
# UNITS OF kg kg-1

def qvVertical_local(jahr):
    import time,sys
    import numpy as np
    from netCDF4 import Dataset,num2date
    from datetime import datetime,timedelta,date

    # tropical subdomain
    lowlat  = -30; hilat  = 30; lowlon = 0; hilon  = 360
    
    # load ERA Interim vertical velocity data
    # div [=] (1460 times, 32 vertical levels, 61 latitudes, 360 longitudes)
    eradata = Dataset('/rigel/home/scs2229/top-secret/MCS_clim/ERAint_qv.nc','r+')
    q = np.asarray(eradata.variables['q'])
    zeit = num2date(eradata.variables['time'][:],'hours since 1900-01-01 00:00:0.0')

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

    # filter out instance where the core position was not recorded (~19% of values)
    tossout = []
    cc = 1
    for ii in range(DD.shape[0]):
        if DD[ii,5] < 0:
           tossout.append(ii)
        elif int(DD[ii,5]) == 360:
           tossout.append(ii)

    DD = np.delete(DD,tossout,axis=0)

    # filter for the subdomain and year here
    hilat = 30; lowlat = -30; hilon = 360; lowlon = 0
    # locate based upon cores
    #lala = DD[:,4]; lolo = DD[:,5]
    # locate based on MCS center
    lala = DD[:,11]; lolo = DD[:,12]

    indx = np.argwhere((lala < hilat) & (lala >= lowlat) & (lolo < hilon) & \
            (lolo >= lowlon) & (DD[:,0] == jahr))
    DD = DD[indx[:,0],:]
 
    # transform the year / month / day / hr format to datetime
    # 32635 instances in 1983 for example
    CTtimes = np.zeros((DD.shape[0]),dtype='datetime64[h]')
    CTtimes = np.asarray([datetime(int(DD[ii,0]),int(DD[ii,1]),int(DD[ii,2]),int(DD[ii,3]),0) for ii in range(DD.shape[0])])
    print 'Starting with ' + str(CTtimes.shape[0]) + ' MCS occurrences in the year of interest.'

    # append each 32 length profile to divARR
    qARR = []
    counter  = 0.0
    for ii in np.arange(DD.shape[0]):
        #if int(counter/CTtimes.shape[0]*100.)%10 == 0:
        #   print(float(counter/CTtimes.shape[0]))
        
        if int(CTtimes[ii].hour)%3 != 0:
           newhour = int(3 * round(float(CTtimes[ii].hour)/3))
           if newhour > 21:
              tomorrow = date(int(CTtimes[ii].year),int(CTtimes[ii].month),\
                         int(CTtimes[ii].day)) + timedelta(days=1)
              obj = datetime(tomorrow.year,tomorrow.month,tomorrow.day,0)
           else:
              obj = datetime(int(CTtimes[ii].year),int(CTtimes[ii].month),\
                          int(CTtimes[ii].day),newhour)
        else:
           obj = CTtimes[ii]

        
        # identify the vertical profiles collocated in time with MCS occurrence
        colloc = np.argwhere((zeit == obj))
        if(len(colloc) != 0):
           qprof  = q[np.asscalar(colloc),:,int(DD[ii,11])-lowlat,int(DD[ii,12])]
           qARR.append(np.concatenate((qprof,DD[ii,:]),axis=0))
           counter += 1.0
        # ERA-Interim profiles are only available at 00:00,06:00,12:00,18:00
        # For 03:00,09:00,15:00,21:00 use the preceding profile
        elif(len(np.argwhere((zeit == obj - timedelta(hours=3)))) != 0):
           colloc = np.argwhere((zeit == obj - timedelta(hours=3)))
           qprof  = q[np.asscalar(colloc),:,int(DD[ii,11])-lowlat,int(DD[ii,12])]
           qARR.append(np.concatenate((qprof,DD[ii,:]),axis=0))
           counter += 1.0
    
    print 'Final size of the specific humidity profile list is ' + str(np.asarray(qARR).shape)
    return np.asarray(qARR)

