#!/bin/sh
##
#SBATCH --account=glab
#SBATCH --job-name=capePRE
#SBATCH -c 1
#SBATCH --time=8:00:00
#SBATCH --mem=4gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scs2229@columbia.edu

source activate ncplot

for year in 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008; do    
    echo Starting for year $year
    cat > ERAint_CAPERetrieve.generated.py <<EOL
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
   'levtype'    : "sfc",
   'param'      : "59.128",
   'step'       : "3/6/9/12",
   'time'       : "00:00:00/12:00:00",
   'type'       : "fc",
   'target'     : "ERAint_cape.nc",
   'format'     : "netcdf",
   'area'       : "30/-180/-30/180"
})
EOL

    sbatch --wait --output=/dev/null /rigel/home/scs2229/top-secret/MCS_clim/submit_CAPErequest.sh
    # once this job is done

    cat > CAPEClim.generated.py <<EOL
#!/usr/bin/env python
import os
import numpy as np
from datetime import datetime,timedelta

execfile('/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/capeCollocate_local_pre.py')
yr = ${year}
w = capeCollocate_local_pre(yr)
np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/cape/post9/CAPE_' + str(yr),w)
EOL
    python CAPEClim.generated.py 
done
