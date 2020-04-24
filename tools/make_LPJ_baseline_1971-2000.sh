#!/bin/bash
# to be run from Timeseries folder

# extract 1971-2000 from baseline simulations and do ymonmean
exp=$1
clim=$2

# paleo file names
varnamesp=(tsurf_aij dtdiurn_aij wsurf_aij pcldt_aij MCamFX_aijl Pressure_aij)
# varnamesp=(tsurf dtdiurn preacc wetd wsurf pcldt lighting)

# baseline file names
varnamesb=(tas dtr sfcWind clt MCamFX Pressure)
# varnamesb=(tas dtr preacc wetdays sfcWind clt lightning)

for ((i=0;i<=6;i++))
	do

			infile=${varnamesp[$i]}$exp
			outfile=${varnamesb[$i]}$clim

			outfile=../../climategenerator_data/climatology_1971-2000/$outfile

			echo $infile $outfile

			cdo -s -f nc4 ymonmean -seldate,1971-01-01,2000-12-31 $infile $outfile

	done
