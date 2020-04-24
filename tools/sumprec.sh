#!/bin/bash

for exp in RCP26 RCP45 RCP60 RCP85
do
  
		infile="pr_GISS-E2-R_"$exp"_all".nc

		original=original/$infile
		outfile=original/"preacc_GISS-E2-R_"$exp"_all".nc
								
		echo $exp $var

		cdo mulc,86400 -muldpm $original $outfile

		ncatted -a units,pr,o,c,"mm mon-1" $outfile 

done
