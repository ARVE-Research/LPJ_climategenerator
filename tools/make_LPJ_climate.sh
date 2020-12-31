#!/bin/bash

# variables we need

# tsurf - surface air temperature (degC)
# dtdiurn : diurnal temperature range (degC)
# prec- precipitation (mm d-1)
# wsurf - surface wind speed (m s-1)
# pcldt - total cloud cover (percent)

# variables that are derived

# preacc - total monthly precipitation (mm mon-1)
# lght - lightning density (from convective mass flux Magi parameterization) (strokes km-2 day-1)
# wetd - number of days with > 0.1mm precipitation (days)

# paleo file names
varnamesp=(tas dtr preacc wetd sfcWind clt lightning)

# baseline file names
varnamesb=(tas dtr preacc wetd sfcWind clt lightning)

prefixp=/group/esd_kaplan/datasets/ModelE/climategenerator_data/CMIP6_lgm/MPI-ESM1-2-LR/rawdata/climean/
prefixb=/group/esd_kaplan/datasets/ModelE/climategenerator_data/climatology_1961-1990/cmip6/
prefixo=/group/esd_kaplan/datasets/ModelE/climategenerator_data/CMIP6_lgm/MPI-ESM1-2-LR/anomalies_1961-1990/

suffix6190=_MPI-ESM1-2-LR_climatology_1961-1990.nc
suffixpaleo=_MPI-ESM1-2-LR_climatology.nc
suffixo=_MPI-ESM1-2-LR_anomaly.nc

# generate anomalies

for ((i=0;i<=6;i++))
do
  echo ${varnamesp[$i]} ${varnamesb[$i]}

  # subtract the 1961-1990 climatology

  cdo sub $prefixp/${varnamesp[$i]}$suffixpaleo $prefixb/${varnamesb[$i]}$suffix6190 $prefixo/nativeres/${varnamesb[$i]}$suffixo

done

# interpolate anomalies to 30 minute (bilinear interpolation ok)

sh /home/akoch/scripts/LPJ_climategenerator/tools/interpolate.sh $prefixo $suffixo

# run script interpolate.sh
