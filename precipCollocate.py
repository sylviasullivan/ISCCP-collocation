#!/usr/bin/env mcsplot
# SAVE THE VECTOR OF (1) MCS SIZE [km^2] (2) DEPTH [K] 
# (3) ACCUMULATED PRECIP [m^3 s-1] AND (4) PRECIP INTENSITY [mm h-1] 
# IN A VECTOR FOR THE PHASE DESIGNATED
# BY ENSO AND FOR THE SUBDOMAIN INDICATED BY lowlat-hilat, lowlon-hilon

# ENSO = 'ENDJF', 'ENJJA', 'LNDJF', 'LNJJA'
# domain = 'am', 'wa', 'cp'
# [-5, 5, 165, 190]; [-25, 5, 290, 325]; [0, 20, -10, 50]

def precipCollocate3(ENSO):
    import time,sys,math
    import numpy as np
    import numpy.ma as ma
    from netCDF4 import Dataset,num2date
    from datetime import datetime,timedelta
    from functools import reduce
   
    lowlat = -33; hilat = 33; lowlon = 0; hilon = 360; domain = 'trop' 
    print(ENSO)
    print(domain)

    # La Nina DJF periods
    if(ENSO == 'lndjf'):
       yrs = [1983,1983,1983,1984,1984,1985,1985,1985,1985,1988,1988,\
              1988,1989,1989,1989,1989,1995,1995,1995,1996,1996,1998,\
              1998,1998,1999,1999,1999,1999,1999,1999,1999,2000,2000,\
              2000,2000,2000,2000,2000,2001,2005,2006,2006,2007,2007,\
              2007,2008,2008,2007,2007]
       mon = ['10','11','12','11','12','01','02','03','04','10','11',\
              '12','01','02','03','04','10','11','12','01','02','10',\
              '11','12','01','02','03','04','10','11','12','01','02',\
              '03','04','10','11','12','01','12','01','02','10','11',\
              '12','01','02','03','04']
       basedir = '../../data/MSWEP-la-nina-yrs/DJF/'
       period  = 'lndjf'

    # El Nino DJF periods
    elif(ENSO == 'endjf'):
         yrs = [1986,1986,1986,1987,1987,1987,1987,1987,1987,1987,1988,1991,\
                1991,1991,1992,1992,1992,1992,1994,1994,1994,1995,1995,1997,\
                1997,1997,1998,1998,1998,1998,2002,2002,2002,2003,2004,2004,\
                2004,2005,2006,2006,2006]
         mon = ['10','11','12','01','02','03','04','10','11','12','01','10',\
                '11','12','01','02','03','04','10','11','12','01','02','10',\
                '11','12','01','02','03','04','10','11','12','01','10','11',\
                '12','01','10','11','12']
         basedir = '../../data/MSWEP-el-nino-yrs/DJF/'
         period  = 'endjf'

    # La Nina JJA periods
    elif(ENSO == 'lnjja'):
         yrs = [1985,1985,1985,1985,1988,1988,1988,1988,1989,1995,1998,1998,\
                1999,1999,1999,1999,1999,2000,2000,2000,2000,2000,2007]
         mon = ['05','06','07','08','06','07','08','09','05','09','08','09',\
                '05','06','07','08','09','05','06','07','08','09','09']
         basedir = '../../data/MSWEP-la-nina-yrs/JJA/'
         period  = 'lnjja'


    # El Nino JJA periods
    else:
         yrs = [1987,1987,1987,1987,1987,1991,1991,1991,1991,1997,1997,1997,\
                1997,2002,2002,2002,2004,2004]
         mon = ['05','06','07','08','09','06','07','08','09','06','07','08',\
                '09','07','08','09','08','09']
         basedir = '../../data/MSWEP-el-nino-yrs/JJA/'
         period  = 'enjja'
    
    # SST values
    erasst = Dataset('/rigel/home/scs2229/top-secret/MCS_clim/data/ERA/ERAint_SST_clim.nc','r+')
    sst = np.asarray(erasst.variables['sst'])   # size of (294, 181, 360)
    tiempo = np.asarray(num2date(erasst.variables['time'][:],'hours since 1900-01-01 00:00:0.0'))
    
    # load CT data
    datei  = "../../data/ISCCP/CT-allDataFilteredEdges.npy"
    CTdata = np.load(datei)
    DD = np.zeros((CTdata.shape[0],12))
    DD[:,0]  = CTdata[:,7]       # years
    DD[:,1]  = CTdata[:,8]       # months
    DD[:,2]  = CTdata[:,9]       # days
    DD[:,3]  = CTdata[:,10]      # hours
    DD[:,4]  = CTdata[:,12]      # lat
    DD[:,5]  = CTdata[:,13]      # lon
    DD[:,6]  = CTdata[:,24]      # depth
    DD[:,7]  = CTdata[:,18]      # number of convective cores
    DD[:,8]  = CTdata[:,40]      # min lat
    DD[:,9]  = CTdata[:,41]      # max lat
    DD[:,10] = CTdata[:,42]      # min lon
    DD[:,11] = CTdata[:,43]      # max lon

    # filter for the subdomain here
    indx1 = np.argwhere(DD[:,4] <= hilat)
    indx2 = np.argwhere(DD[:,4] >= lowlat)
    indx3 = np.argwhere(DD[:,5] <= hilon)
    indx4 = np.argwhere(DD[:,5] >= lowlon)
    indx5 = np.argwhere(DD[:,5] != 360)
    indx  = reduce(np.intersect1d,(indx1,indx2,indx3,indx4,indx5))
    DD2   = DD[indx,:]
    del DD

    # calculate and save an array of precipitation intensities collocated with MCS
    depth = []; psum = []; pmax = []; ncc = []   # szs = []; life = []; ecc = []
    for ii in range(len(yrs)):
        print(str(yrs[ii]) + mon[ii])
        ppp   = []
        mswep = Dataset(basedir + 'MSWEP_prec' + str(yrs[ii]) + mon[ii] + '_3h_halfdeg.nc','r+')
        prec  = mswep.variables['precipitation']
        prec  = np.divide(prec,3.)                 # divide by 3 hours for hourly intensity
        zeit  = num2date(mswep.variables['time'][:],'days since 1899-12-31 00:00:00')       
 
        # find the MCS that correspond to this MSWEP time period
        for jj in range(zeit.shape[0]):
            indx1  = np.argwhere(DD2[:,0] == zeit[jj].year)
            indx2  = np.argwhere(DD2[:,1] == zeit[jj].month)
            indx3  = np.argwhere(DD2[:,2] == zeit[jj].day)
            indx4  = np.argwhere(DD2[:,3] == zeit[jj].hour)
            indx   = reduce(np.intersect1d,(indx1,indx2,indx3,indx4))  
            subset = DD2[indx,:]

            # extract the maximum and minimum lat and lon of the system
            for kk in range(len(indx)):
                lat = [int(round(subset[kk,8]*2)+180), int(round(subset[kk,9]*2)+180)]
                lon = [int(round(subset[kk,10]*2)), int(round(subset[kk,11]*2))]
                lon[lon == 720] = 0
              
                # define a grid with [lat] x [lon] dimensions
                grid = np.zeros((len(lat),len(lon)))

                # and store precipitation values for the max/min-lat/lon domain
                pr = prec[jj,:,:] 
                for xx in range(len(lat)):
                    rr = lat[xx]
                    if rr == 360:
                       rr = 0
                    ss = lon[0]
                    for yy in range(len(lon)):
                        ss = lon[yy]
                        if(ss == 720):
                           ss = 0
                        grid[xx,yy] = pr[rr,ss]
                
                # calculate the SST from the center lat/lon and time
                sstval = sst[np.argwhere(tiempo == datetime(int(subset[kk,0]),\
                         int(subset[kk,1]),1)),int(subset[kk,4])+90,\
                         int(subset[kk,5])]

                # invalid sea surface temperatures are stored as -32767
                if sstval and np.asscalar(sstval) > 200.:
                   # calculate cumulative and max precip from the grid, store these and MCS values
		   # assume that each 0.5 deg x 0.5 deg cell has an area of 55^2 km2
	           # 10/36. converts from mm h-1 km2 to m3 s-1
		   ppp = np.nansum(grid[:])*3025*10/36.*len(lat)*len(lon)
		   qqq = np.nanmax(grid[:])
		   if int(ppp) != 0:
	              depth.append(sstval - subset[kk,6])
		      ncc.append(subset[kk,7])
                      psum.append(ppp)
		      pmax.append(qqq)   

    # store the system depth - SST [K], precip vol [m3 s-1], and max precip rate [mm h-1]
    # store the system conv fraction, number of conv cores, precip vol [m3 s-1], and max precip rate [mm h-1]
    mcsprec = np.zeros((4,len(depth)))
    for ii in range(len(depth)):
        mcsprec[0,ii] = depth[ii] 
        mcsprec[1,ii] = ncc[ii]
        mcsprec[2,ii] = psum[ii]
        mcsprec[3,ii] = pmax[ii]
    print(mcsprec.shape)

    print('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/precip_clim/MSWEP/' + \
          domain + '_depth-SST_psum_pmax_' + period + '.npy')
    np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/precip_clim/MSWEP/' + \
          domain + '_depth-SST_psum_pmax_' + period + '.npy',mcsprec)

