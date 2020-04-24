#!/bin/bash

# extract 1961-1990 from baseline simulations and do ymonmean

for exp in RCP26 RCP45 RCP60 RCP85
do

  for var in tas dtr preacc wetdays clt sfcWind lightning
  do
  
      infile=$var"_GISS-E2-R_"$exp"_all".nc
  
      transient=original/$infile
      baseline=climatology_1961-1990/$infile
      
      outfile=anomalies/$infile
      
      echo $exp $var

      cdo -s -f nc4 ymonsub $transient $baseline $outfile
  
  done
done
