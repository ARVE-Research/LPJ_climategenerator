#!/bin/bash

# extract 1961-1990 from baseline simulations and do ymonmean

for exp in RCP26 RCP45 RCP60 RCP85
do
		for var in tas dtr preacc wetdays clt sfcWind lightning
		do

				infile=original/$var"_GISS-E2-R_"$exp"_all".nc

				outfile=climatology_1961-1990/${infile##*/}

				echo $var $infile $outfile

				cdo -s -f nc4 ymonmean -seldate,1961-01-01,1990-12-31 $infile $outfile
		
		done
done
