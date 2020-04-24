#!/bin/bash

for exp in RCP26 RCP45 RCP60 RCP85
do

		cdo -f nc4 sub tasmax_GISS-E2-R_$exp"_all".nc tasmin_GISS-E2-R_$exp"_all".nc dtr_GISS-E2-R_$exp"_all".nc 

done
