#!/usr/bin/env ncplot
# COLLOCATE THE CAPE for the INPUT YEAR

def capeCollocate_ERA5(jahr,mois):
    import time,sys
    import numpy as np
    from netCDF4 import Dataset,num2date
    from datetime import datetime,timedelta,date

    # tropical subdomain
    lowlat  = -30; hilat  = 30; lowlon = 0; hilon = 360
    
    # load ERA-5 CAPE data
    # www [=] (1460 times, 32 vertical levels, 61 latitudes, 360 longitudes)
    eradata = Dataset('/rigel/home/scs2229/top-secret/MCS_clim/ERA5_cape.nc','r+')
    alal = np.asarray(eradata.variables['latitude'])
    olol = np.asarray(eradata.variables['longitude'])
    cape = np.asarray(eradata.variables['cape'])
    zeit = num2date(eradata.variables['time'][:],'hours since 1900-01-01 00:00:0.0')

    # load CT data
    datei  = '/rigel/home/scs2229/top-secret/MCS_clim/data/ISCCP/CT-allDataLocalTimes.npy'
    CTdata = np.load(datei)
    DD = np.zeros((CTdata.shape[0],20))
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
    DD[:,16] = CTdata[:,40]     # min latitude
    DD[:,17] = CTdata[:,41]     # max latitude
    DD[:,18] = CTdata[:,42]     # min longitude
    DD[:,19] = CTdata[:,43]     # max longitude

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

    # filter for systems in the tropics and in the year and month of interest
    indx = np.argwhere((lala < hilat) & (lala >= lowlat) & (lolo < hilon) & \
            (lolo >= lowlon) & (DD[:,0] == jahr) & (DD[:,1] == mois))
    DD = DD[indx[:,0],:]

    # transform the year / month / day / hr format to datetime
    # 32635 instances in 1983 for example
    CTtimes = np.zeros((DD.shape[0]),dtype='datetime64[h]')
    CTtimes = np.asarray([datetime(int(DD[ii,0]),int(DD[ii,1]),int(DD[ii,2]),int(DD[ii,3]),0) for ii in range(DD.shape[0])])
    print 'Starting with ' + str(CTtimes.shape[0]) + ' MCS occurrences in the year of interest.'

    # append the four temperatures to tARR
    tARR = []
    counter  = 0.0
    for ii in np.arange(DD.shape[0]):
        #print(np.float(ii)/np.float(DD.shape[0]))
        # the hour no longer has to be on a 3-hourly grid 
        obj = CTtimes[ii]

        # identify the vertical profiles preceding by 6 hrs MCS occurrence
        if((obj - timedelta(hours=6)).month != mois):
           continue
        else:
           colloc = np.argwhere((zeit == obj - timedelta(hours=6)))
           #if(len(colloc) != 0):
           # collocating with system (11,12) versus core (4,5) here
           # find the nearest latitude to 0.25
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
           ###### SOLUTION IF LONGITUDES ARE FROM 0 AT INDX 0 to 360 AT INDX 1440 #######
           ## find the nearest longitude to 0.25; u is rounding up, d is rounding down
           #u = np.ceil(DD[ii,5]*4)/4; d = np.floor(DD[ii,5]*4)/4
           #if np.abs(u - DD[ii,5]) < np.abs(d - DD[ii,5]):
           #   lon25 = np.ceil(DD[ii,5]*4)
           #else:
           #   lon25 = np.floor(DD[ii,5]*4)
           ###### SOLUTION IF LONGITUDES ARE FROM 0 AT INDX 0 to 360 AT INDX 1440 #######
           #print np.asscalar(colloc),DD[ii,4],alal[int(lat25)],DD[ii,5],DDlon,olol[int(lon25)]
           cc = cape[np.asscalar(colloc),int(lat25),int(lon25)]
           tARR.append(np.insert(DD[ii,:],0,cc,axis=0))
           counter += 1.0
           # ERA-5 profiles should be available every hour so that colloc always exists
           #elif(len(np.argwhere((zeit == obj - timedelta(hours=6)))) != 0):
           #   print colloc
           #   colloc = np.argwhere((zeit == obj - timedelta(hours=7)))
           #   # collocating with system here (11,12) versus core (4,5)
           #   cc = cape[np.asscalar(colloc),int(DD[ii,4])-lowlat,int(DD[ii,5])]
           #   tARR.append(np.insert(DD[ii,:],0,cc,axis=0))
           #   counter += 1.0

    print 'Final size of the CAPE profile list is ' + str(np.asarray(tARR).shape)
    return np.asarray(tARR)

