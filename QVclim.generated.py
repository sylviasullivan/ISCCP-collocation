#!/usr/bin/env python
import os
import numpy as np
from datetime import datetime,timedelta

execfile('/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/qvVertical_local_pre.py')
yr = 1983
w = qvVertical_local_pre(yr)
np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/qv/core_pre9/qv_' + str(yr),w)
