# Collocation scripts for ISCCP / MSWEP / ERA-I analyses

The procedure to collocate ERA reanalysis values with the mesoscale convective systems in the ISCCP CT database is as follows.

1. *submit_\<var\>\<vers\>.sh* - This script iterates through years from 1983 to 2008 to produce *ERA\<vers\>_\<var\>.nc*. It uses the Climate Data Store API to make the request for variable *var* and for reanalysis version *vers*. Two examples are uploaded: submit_cape5.sh as a single-level variable example and submit_qv5.sh as a pressure-level example. Scripts to download from ERA-Interim have a blank \<vers\>, and those to download from ERA-5 have a \<vers\> = 5. After this request is submitted in the slurm system, a python file is generated to do collocation for this year before overwriting the *ERA\<var\>_\<var\>.nc* for the next year.

2. *submit_\<var\>\<vers\>request.sh* - Within *submit_\<var\>\<vers\>.sh*, this script submits the Climate Data Store API request each year. 
  
3. *\<var\>\<vers\>Clim.generated.py* - A python file written within *submit_\<var\>\<vers\>.sh* to perform collocation between the full-year ERA file and the intermittent systems in the ISCCP CT database.
  
4. *\<var\>Vertical_ERA\<vers\>.py* - A python file that identifies values within the ERA nc file corresponding spatiotemporally to the convective system in the ISCCP CT database. xarray is required for the ERA-5 variables on pressure levels as they are ~10 GB per year.

# Visualization scripts 

Collocation of MSWEP precipitation data with ISCCP convective systems

*fig1v12-precip-hist.py* -- Histograms of relative difference in precipitation intensities stratified by system depth and El Nino phase

*fig2-SF-GMS-winds.py* -- Surface wind profiles, surface flux and net gross moist stability distributions

*fig3-upd-pgf.py* -- Updraft profiles and distributions of the pressure gradient force during El Nino versus La Nina

*fig4-rh-zbp.py* -- RH and zero-buoyancy plume profiles

*figS1-zdist.py* -- Convective depth distributions 

*figS2-abs-hist.py* -- Histogram of precipitation intensities stratified by system depth and El Nino phase

*figS3-thresholds.py*, *figS4-binning.py*, *figS5-bounds.py* -- Sensitivity of histogram of relative difference in precipitation intensities to depth thresholds, binning, and min/max bounds

*figS6-req-depth.py* -- How does equivalent radius covary with system depth in El Nino versus La Nina?
