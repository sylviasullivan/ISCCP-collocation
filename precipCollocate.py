#!/usr/bin/env mcsplot
# SAVE THE VECTOR OF (1) MCS DEPTH [K] (2) LOWEST TEMPERATURE LEVEL FROM TEMP PROFILES [K] 
# (3) PRECIP INTENSITY IN CONVECTIVE CORE [mm h-1] 
# IN A VECTOR FOR THE PHASE DESIGNATED
# BY ENSO AND FOR THE SUBDOMAIN INDICATED BY lowlat-hilat, lowlon-hilon

# ENSO = 'ENDJF', 'ENJJA', 'LNDJF', 'LNJJA'
# domain = 'am', 'wa', 'cp'
# [-5, 5, 165, 190]; [-25, 5, 290, 325]; [0, 20, -10, 50]

def precipCollocate6_local(phase):
    import time,sys,math
    import numpy as np
    import numpy.ma as ma
    from netCDF4 import Dataset,num2date
    from datetime import datetime,timedelta,date
    from functools import reduce
   
    if phase == 'lndjf':
       monperyr = [[10,11,12],[11,12],[1,2,3,4],[],[],[10,11,12],[1,2,3,4],[],[],\
                  [],[],[],[10,11,12],[1,2],[],[10,11,12],[1,2,3,4,10,11,12],\
                  [1,2,3,4,10,11,12],[1],[],[],[],[12],[1,2],[4,10,11],[]]
    elif phase == 'endjf':
         monperyr = [[],[],[],[10,11,12],[1,2,3,4,10,11,12],[1],[],[],[10,11,12],\
                     [1,2,3,4],[],[10,11,12],[1,2],[],[10,11,12],[1,2,3,4],\
                     [],[],[],[10,11,12],[1],[10,11,12],[1],[10,11,12],[],[]]
    basedir = '/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/surftemp/'
    if phase == 'endjf':
       basedir2 = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-el-nino-yrs/DJF/'
    elif phase == 'lndjf':
       basedir2 = '/rigel/home/scs2229/top-secret/MCS_clim/data/MSWEP-la-nina-yrs/DJF/'
    else:
       print 'Phase mis-specified.'
       sys.exit()

    # calculate and save an array of precipitation intensities collocated with these systems
    t2m = []; skt = []; sst = []; dpt = []; precip = []; ctt = []; psum = []
    for jahr in np.arange(1983,2009):
        # temperature and convective system data
        tprof = np.load(basedir + 'SurfTemp_' + str(jahr) + '.npy')

        # put the hours on 3-hourly grid
        for j,row in enumerate(tprof):
            if int(row[7])%3 != 0:
               newhour = int(3 * round(float(row[7])/3))
               if newhour > 21:
                  tomorrow = date(int(row[4]),int(row[5]),int(row[6])) + timedelta(days=1)
                  obj = datetime(tomorrow.year,tomorrow.month,tomorrow.day,0)
               else:
                  obj = datetime(int(row[4]),int(row[5]),int(row[6]),newhour)
               tprof[:,4] = obj.year
               tprof[:,5] = obj.month
               tprof[:,6] = obj.day
               tprof[:,7] = obj.hour

        for mon in monperyr[jahr-1983]:
            if mon < 10:
               monstr = '0' + str(mon)
            else:
               monstr = str(mon)
            print basedir2 + 'MSWEP_prec' + str(jahr) + monstr + '_3h_halfdeg.nc'
            mswep = Dataset(basedir2 + 'MSWEP_prec' + str(jahr) + monstr + '_3h_halfdeg.nc','r+')
            prec  = mswep.variables['precipitation']
            #prec  = np.divide(prec,3.)                 # divide by 3 hours for hourly intensity
            zeit  = num2date(mswep.variables['time'][:],'days since 1899-12-31 00:00:00')       
 
            # find the convective systems with the same month, day, hour in both
            # the temp profile and sst values
            for jj in range(zeit.shape[0]):
                indx1 = np.argwhere(tprof[:,5] == zeit[jj].month)
                indx2 = np.argwhere(tprof[:,6] == zeit[jj].day)
                indx3 = np.argwhere(tprof[:,7] == zeit[jj].hour)
                indx = reduce(np.intersect1d,(indx1,indx2,indx3))  
                subset = tprof[indx,:]
         
                # extract the values at the closest 0.5 deg to the core lat, lon
                for kk in range(len(indx)):
                    #lat = int(round(subset[kk,8]*2.0)) + 180
                    #lon = int(round(subset[kk,9]*2.0))# + 180
                    #if lon != 720:
                    #   pr = prec[jj,lat,lon]
                    lat = [int(round(subset[kk,20]*2)+180),int(round(subset[kk,21]*2)+180)]
                    lon = [int(round(subset[kk,22]*2)),int(round(subset[kk,23]*2))]
                
                    # define a grid with [lat] x [lon] dimensions
                    grid = np.zeros((len(lat),len(lon)))
 
                    # and store precipitation values for the max/min-lat/lon domain
                    pr = prec[jj,:,:]
                    for xx in np.arange(len(lat)):
                        rr = lat[xx]
                        if rr == 360:
                           rr = 0
                        ss = lon[0]
                        for yy in np.arange(len(lon)):
                            ss = lon[yy]
                            if ss == 720:
                               ss = 0
                            grid[xx,yy] = pr[rr,ss]

                    # save the 2-m, skin, sea surface, and dew point temperatures
                    t2m.append(subset[kk,0])
                    skt.append(subset[kk,1])
                    sst.append(subset[kk,2])
                    dpt.append(subset[kk,3])
                    # append the cloud top temperature to ctt
		    ctt.append(subset[kk,19])
                    # append the precipitation max   
                    precip.append(np.nanmax(grid[:]))
                    # append the precipitation accumulation
                    psum.append(np.nansum(grid[:])*3025*10/36.*len(lat)*len(lon))

    # store the system cloud top temp [K], SST [K], low-level temperature [K], max precip rate [mm h-1]
    # precip accumulation [m3 s-1]
    mcsprec = np.zeros((7,len(ctt)))
    for ii in range(len(ctt)):
        mcsprec[0,ii] = t2m[ii] 
        mcsprec[1,ii] = skt[ii]
        mcsprec[2,ii] = sst[ii]
        mcsprec[3,ii] = dpt[ii]
        mcsprec[4,ii] = ctt[ii]
        mcsprec[5,ii] = precip[ii]
        mcsprec[6,ii] = psum[ii]
    print(mcsprec.shape)

    print('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/precip_clim/MSWEP/' + \
          'trop_temps_pmax_psum_' + str(phase) + '.npy')
    np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/precip_clim/MSWEP/' + \
          'trop_temps_pmax_psum_' + str(phase) + '.npy',mcsprec)

