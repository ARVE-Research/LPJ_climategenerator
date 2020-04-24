#!/bin/bash

# GISS E2-R variable names

# variables we don't need

# pdcld : Deep convective cloud cover
# mccltp : Convective cloud top pressure
# mccldbs = convective cloud base pressure
# TMAXC - SURFACE AIR TEMPERATURE (DIURNAL MAX)
# TMINC - SURFACE AIR TEMPERATURE (DIURNAL MIN)
# TMNMX - SURFC AIR TEMPERATURE (LOWEST DIURNAL HIGH)

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
varnamesp=(tsurf dtdiurn preacc wetd wsurf pcldt Magi_lightning)

# baseline file names
varnamesb=(tas dtr preacc wetdays sfcWind clt lightning)

prefixp=/home/ucfaako/Documents/Allegra_GISS/data/giss_climate/4xCO2/rawdata/climean
prefixb=/home/ucfaako/Documents/Allegra_GISS/data/giss_climate/climatology_1961-1990
prefixo=/home/ucfaako/Documents/Allegra_GISS/data/giss_climate/4xCO2/anomalies_1961-1990

suffix6190=_GISS-E2-R_1961-1990_climmean.nc
suffixpaleo=_CEN_Clim_2900-2999.aijE2p1_anl_4xCO2.nc 
suffixo=_GISS-E2-R_anomaly.nc

# interpolate anomalies to 30 minute (bilinear interpolation ok)

sh interpolate.sh $prefixo $suffixo

# add to baseline
for ((i=0;i<=6;i++))
do
  echo ${varnamesp[$i]} ${varnamesb[$i]}
 
  cdo ymon $prefixp/${varnamesp[$i]}$suffixpaleo $prefixb/${varnamesb[$i]}$suffix6190 $prefixo/nativeres/${varnamesb[$i]}$suffixo
  
done

