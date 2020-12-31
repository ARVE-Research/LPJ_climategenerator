#!/bin/bash

# ---
for MODEL in MIROC-ES2L;
do
datadir=/group/esd_kaplan/datasets/ModelE/climategenerator_data/CMIP6_lgm/${MODEL}/rawdata/climean/
codedir=/home/akoch/LPJ/climategenerator/tools
coeffile=$codedir/pre2wet_coeffs_rotated.nc
suffix=_${MODEL}_climatology.nc
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

prefile=$datadir/pr$suffix
preaccfile=$datadir/preacc$suffix
outfile=$datadir/wetd$suffix

#cdo muldpm -mulc,86400 $prefile $preaccfile # OR ##
cdo muldpm $prefile $preaccfile
ncrename -v pr,prec $preaccfile
ncgen -k nc4 -o $outfile $codedir/wetdays.cdl

$codedir/calcwetdays $preaccfile $coeffile $outfile
done
