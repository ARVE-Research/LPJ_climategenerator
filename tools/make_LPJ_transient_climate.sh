#!/bin/bash
# To be run in /rawdata/climean;
# for transient files
exp=$1
clim=$2

# 1961-1990 anomalies
for var in tas dtr preacc wetdays clt sfcWind lightning
do

  infile=${var}${exp}
  climfile=${var}${clim}

  baseline=../../../climatology_1961-1990/$climfile

  outfile=../../anomalies_1961-1990/nativeres/$infile

  echo $var

  cdo -s -f nc4 ymonsub $infile $baseline $outfile

done
