#!/bin/sh
#
#SBATCH --account=glab
#SBATCH --job-name=qv5all
#SBATCH -c 1
#SBATCH --time=10:00:00
#SBATCH --mem=20gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scs2229@columbia.edu

source activate ncplot

for year in 2006; do
    for month in 3; do
        echo Starting for year $year , month $month
        cat > ERA5_QVRetrieve.generated.py <<EOL
import cdsapi
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'variable':'specific_humidity',
        'pressure_level':['125','150','175','200','225',
                          '250','300','350','400','450',
                          '500','550','600','650','700',
                          '750','775','800','825','850',
                          '875','900','925','950','975',
                          '1000'],
        'year':${year},
        'month':'${month}',
        'day':['01','02','03','04','05','06','07','08','09',
               '10','11','12','13','14','15','16','17','18',
               '19','20','21','22','23','24','25','26','27',
               '28','29','30','31'],
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
 './ERA5_qv.nc')
EOL

    sbatch --wait /rigel/home/scs2229/top-secret/MCS_clim/submit_QV5request.sh
    # once this job is done

    cat > QV5Clim.generated.py <<EOL
#!/usr/bin/env python
import sys
import os
import numpy as np
from datetime import datetime,timedelta

basedir = '/rigel/home/scs2229/top-secret/MCS_clim/scripts/meteo/'
sys.path.insert(1,basedir)
execfile(basedir + 'qvVertical_ERA5_xarray.py')
yr = ${year}
mon = ${month}
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
EOL
    python QV5Clim.generated.py 
done
wait
done
