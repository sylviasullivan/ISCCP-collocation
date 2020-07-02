#!/usr/bin/env python
import sys
import os
import numpy as np
from datetime import datetime,timedelta

basedir = '/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/'
sys.path.insert(1,basedir)
execfile(basedir + 'qvVertical_ERA5_xarray.py')
yr = 2006
mon = 3
if mon < 10:
   monstr = '0' + str(mon)
else:
   monstr = str(mon)
basedir = '/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/qv/'
qv = qvVertical_ERA5_xarray(yr,mon,0)
fichier = basedir + 'ERA5/qv_' + str(yr) + monstr + '.npy'
np.save(fichier,qv)
qv = qvVertical_ERA5_xarray(yr,mon,3)
fichier = basedir + 'ERA5_pre/qv_' + str(yr) + monstr + '.npy'
np.save(fichier,qv)
qv = qvVertical_ERA5_xarray(yr,mon,9)
fichier = basedir + 'ERA5_pre9/qv_' + str(yr) + monstr + '.npy'
np.save(fichier,qv)
os.remove('/rigel/home/scs2229/top-secret/MCS_clim/ERA5_qv.nc')
