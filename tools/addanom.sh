#!/bin/bash

title=GISS_E2R_climatology_021ka

jobfile=/group/esd_kaplan//datasets/ModelE/climategenerator_data/006ka/giss021ka.namelist
outfile=/group/esd_kaplan//datasets/ModelE/climategenerator_data/006ka/baseline/$title.nc

# generate an empty output file with the same format as the baseline file

sed 's/TITLE/$title/g' lpj_climate_input_30m.cdl > tmp.cdl

ncgen -o $outfile tmp.cdl

# add the paleoclimate anomaly climatology to the baseline climatology

./addanom $jobfile $outfile
