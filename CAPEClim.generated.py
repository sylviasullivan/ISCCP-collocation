#!/usr/bin/env python
import os
import numpy as np
from datetime import datetime,timedelta

execfile('/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/capeCollocate_local_pre.py')
yr = 2008
w = capeCollocate_local_pre(yr)
np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/cape/post9/CAPE_' + str(yr),w)
