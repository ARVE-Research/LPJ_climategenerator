README FOR GENERATING LPJ CLIMATE INPUT FILES FOR A PALEO TIME-SLICE SIMULATION

1. Collect GCM output for the following variables (monthly mean unless specified)
  a. 2m air temperature (degC)
  b. diurnal temperature range (degC)
  c. monthly total precipitation (mm)
  d. days with precip. > 0.1 mm (optional) (days)
  e. total cloud cover (percent)
  f. 10m windspeed (m s-1)
  g. lightning stroke density (optional as long as you have a precursor variable) (km-2 d-1)

1a. GCM output is delivered as a timeseries, calculate the climatological mean (multi-year monthly mean) of the above variables (12 months) (cdo ymonmean)

1b. Where necessary, calculate the derived variables: (all steps calcwetdays.sh)
    - total monthly precipitation (cdo muldpm)
    - wet days (using calcwetdays.sh and calcwetdays.f90)
    - lightning on the basis of convective mass flux (Magi parameterization - magi_lightning.py; or other formula)

2. Generate de-biased GCM anomalies by subtracting GCM paleoclimate climatology from GCM baseline (2nd half of 20th century) climatology (12 months) (make_LPJ_climate.sh)

3. Interpolate GCM anomalies to 0.5 degree, this is the paleoclimate anomaly (paleo-anomalies, 12 months) (interpolate.sh) # This happens in make_LPJ_climate.sh

4. Add the paleo-anomalies to the hi-res present-day baseline climatology (climate_wwna_wglc_shelves.nc), this becomes the "basefile" (12 months) (NB check units between GCM and baseline) (addanom.f90 in addanom.sh)

5. Create a timeseries of interannually variable climate by adding the basefile to detrended interannual variability from 20CR (makeclimate.f90 use ./makeclimate $jobfile $output)
