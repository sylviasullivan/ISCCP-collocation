#!/bin/sh
##
#SBATCH --account=glab
#SBATCH --job-name=cape5
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scs2229@columbia.edu

source activate ncplot

for year in 1995; do
#1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008; do
    for month in 5; do
        echo Starting for year $year
        cat > ERA5_CAPERetrieve.generated.py <<EOL
import cdsapi
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'variable':'convective_available_potential_energy',
        'year':${year},
        'month':'0${month}',
        'day':['01','02','03','04','05','06','07',
               '08','09','10','11','12','13','14',
               '15','16','17','18','19','20','21',
               '22','23','24','25','26','27','28',
               '29','30','31'],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'],
        'area':'30/-180/-30/180',
        'format':'netcdf'
    },
 './ERA5_cape.nc')
EOL

    sbatch --wait --output=/dev/null /rigel/home/scs2229/top-secret/MCS_clim/submit_CAPE5request.sh
    # once this job is done

    cat > CAPE5Clim.generated.py <<EOL
#!/usr/bin/env python
import os
import numpy as np
from datetime import datetime,timedelta

execfile('/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/capeCollocate_ERA5.py')
yr = ${year}
mon = ${month}
w = capeCollocate_ERA5(yr,mon)
if mon < 10:
   monstr = '0' + str(mon)
else:
   monstr = str(mon)
np.save('/rigel/home/scs2229/top-secret/MCS_clim/ausgabe/meteo_clim/cape/ERA5_pre6/CAPE_' + str(yr) + monstr,w)
EOL
    python CAPE5Clim.generated.py 
done
done
