#!/bin/bash

# ---

datadir=/group/esd_kaplan/datasets/ModelE/climategenerator_data/4xCO2/rawdata/climean
codedir=/home/akoch/scripts/LPJ_climategenerator/tools
coeffile=$codedir/pre2wet_coeffs_rotated.nc
suffix=_CEN_Clim_2900-2999.aijE2p1_anl_4xCO2.nc
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

prefile=$datadir/prec$suffix
preaccfile=$datadir/preacc$suffix
outfile=$datadir/wetd$suffix

cdo muldpm $prefile $preaccfile

ncgen -k nc4 -o $outfile $codedir/wetdays.cdl

$codedir/calcwetdays $preaccfile $coeffile $outfile
