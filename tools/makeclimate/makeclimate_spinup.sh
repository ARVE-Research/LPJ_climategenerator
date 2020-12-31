#!/bin/bash

template=climate_template.cdl

tlen=1020

jobfile=${1}
output=${2}

echo $basefile $output

sed 's/TCHUNK/'$tlen'/g' $template | ncgen -k nc4 -o $output

./makeclimate $jobfile $output
