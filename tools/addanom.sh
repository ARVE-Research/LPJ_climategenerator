#!/bin/bash

title=${1}

jobfile=${2}
outfile=${3}

# generate an empty output file with the same format as the baseline file

sed "s/TITLE/$title/g" lpj_climate_input_30m.cdl > tmp.cdl

ncgen -o $outfile tmp.cdl

# add the paleoclimate anomaly climatology to the baseline climatology

./addanom $jobfile $outfile
