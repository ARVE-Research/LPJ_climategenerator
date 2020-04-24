#!/bin/bash
# interpolate anomalies to 0.5 degree

prefix=$1

for var in *.nc
do

    original=$prefix/nativeres/$var
    outfile=$prefix/30m/$var
          
    echo $var

    cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

done
