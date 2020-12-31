# Generating climate input files for LPJ-LMfire

1. Collect GCM output for the following variables (monthly mean unless specified):
    1. 2m air temperature (degC)
    2. diurnal temperature range (degC)
    3. monthly total precipitation (mm)
    4. days with precip. > 0.1 mm (optional) (days)
    5. total cloud cover (percent)
    6. 10m windspeed (m s-1)
    7. lightning stroke density (optional as long as you have a precursor variable) (km-2 d-1)

2. If the GCM output is delivered as a timeseries, calculate the climatological mean (multi-year monthly mean) of the above variables (12 months) `cdo ymonmean`

3. Where necessary, calculate the derived variables: (all steps `calcwetdays.sh`)
    - total monthly precipitation (`cdo muldpm`)
    - wet days (using `calcwetdays.sh` and `calcwetdays.f90`)
    - lightning on the basis of convective mass flux (Magi parameterization - `magi_lightning.py`; or other formula)

4. Generate de-biased GCM anomalies by subtracting GCM paleoclimate climatology from GCM baseline (2nd half of 20th century) climatology (12 months) `make_LPJ_climate.sh`

5. Interpolate GCM anomalies to 0.5 degree, this is the paleoclimate anomaly (paleo-anomalies, 12 months) (`interpolate.sh`) # This happens in `make_LPJ_climate.sh`

6. Add the paleo-anomalies to the hi-res present-day baseline climatology `climate_wwna_wglc_shelves.nc`, this becomes the "basefile" (12 months) (NB check units between GCM and baseline) (`addanom.f90` in `addanom.sh`)

7. Create a timeseries of interannually variable climate by adding the basefile to detrended interannual variability from 20CR (`makeclimate.f90` use `./makeclimate $jobfile $output`)
