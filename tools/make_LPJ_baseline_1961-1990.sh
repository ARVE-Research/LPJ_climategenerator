#!/bin/bash
# extract 1961-1990 from baseline simulations and do ymonmean
exp=$1
clim=$2
cmip=$3

# paleo file names
varnamesp=(tas dtr preacc wetd sfcWind clt lightning)

# baseline file names
varnamesb=(tas dtr preacc wetdays sfcWind clt lightning)

for i in "${!varnamesp[@]}"
	do

			infile=${varnamesp[$i]}$exp
			outfile=${varnamesb[$i]}$clim

			outfile=/group/esd_kaplan/datasets/ModelE/climategenerator_data/climatology_1961-1990/$cmip/$outfile

			echo $infile $outfile

			cdo -s -f nc4 ymonmean -seldate,1961-01-01,1990-12-31 $infile $outfile

	done
