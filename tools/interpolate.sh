#!/bin/bash

prefix=$1
suffix=$2

# interpolate anomalies to 0.5 degree

for var in tas dtr preacc wetdays clt sfcWind lightning
do

    infile=$var$suffix

    original=$prefix/nativeres/$infile
    outfile=$prefix/30m/$infile
          
    echo $var

    cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

done
