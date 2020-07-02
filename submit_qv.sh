#!/bin/sh
##
#SBATCH --account=glab
#SBATCH --job-name=humidity
#SBATCH -c 1
#SBATCH --time=10:00:00
#SBATCH --mem=20gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scs2229@columbia.edu

source activate ncplot

for year in 1995 1994 1993 1992 1991 1990 1989 1988 1987 1986 1985 1984 1983; do
    echo Starting for year $year
    cat > ERAint_QVRetrieve.generated.py <<EOL
#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import os

yr = ${year}
datestring = str(yr) + "-01-01/to/" + str(yr) + "-12-31"
print(datestring)
server = ECMWFDataServer()

server.retrieve({
   'class'      : "ei",
   'dataset'    : "interim",
   'date'       : datestring,
   'expver'     : "1",
   'grid'       : "1/1",
   'stream'     : "oper",
   'levtype'    : "ml",
   'levellist'  : "26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57",
   'param'      : "133.128",
   'step'       : "0",
   'time'       : "00:00:00/06:00:00/12:00:00/18:00:00",
   'type'       : "an",
   'target'     : "ERAint_qv.nc",
   'format'     : "netcdf",
   'area'       : "30/-180/-30/180"
})
EOL

    sbatch --wait --output=/dev/null /rigel/home/scs2229/top-secret/MCS_clim/submit_QVrequest.sh
    # once this job is done
    cat > QVclim.generated.py <<EOL
#!/usr/bin/env python
import os
import numpy as np
from datetime import datetime,timedelta

execfile('/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/qvVertical_local_pre.py')
yr = ${year}
w = qvVertical_local_pre(yr)
np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/qv/core_pre9/qv_' + str(yr),w)
EOL
    python QVclim.generated.py
done
