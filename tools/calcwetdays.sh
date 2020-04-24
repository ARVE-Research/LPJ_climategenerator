#!/bin/bash

# ---

datadir=/group/esd_kaplan/datasets/ModelE/climategenerator_data/climatology_1961-1990
codedir=/data/akoch/lpjlm/LPJ_climategenerator/LPJ_climategenerator_code/tools
coeffile=$codedir/pre2wet_coeffs_rotated.nc

#for exp in RCP26 RCP45 RCP60 RCP85
#do
#
#  echo $exp
#
#  prefile=$datadir/"preacc_GISS-E2-R_"$exp"_all".nc
#  outfile=$datadir/"wetdays_GISS-E2-R_"$exp"_all".nc
#
#  ncgen -k nc4 -o $outfile wetdays.cdl
#
#  ./calcwetdays $prefile $coeffile $outfile
#
#done

prefile=$datadir/prec_GISS-E2-R_climatology_1961-1990.nc
preaccfile=$datadir/preacc_GISS-E2-R_climatology_1961-1990.nc
outfile=$datadir/wetd_GISS-E2-R_climatology_1961-1990.nc

cdo muldpm $prefile $preaccfile

ncgen -k nc4 -o $outfile $codedir/wetdays.cdl

$codedir/calcwetdays $preaccfile $coeffile $outfile
