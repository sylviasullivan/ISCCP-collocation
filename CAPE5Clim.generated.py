#!/usr/bin/env python
import os
import numpy as np
from datetime import datetime,timedelta

execfile('/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/capeCollocate_ERA5.py')
yr = 1995
mon = 5
w = capeCollocate_ERA5(yr,mon)
if mon < 10:
   monstr = '0' + str(mon)
else:
   monstr = str(mon)
np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/cape/ERA5_pre6/CAPE_' + str(yr) + monstr,w)
