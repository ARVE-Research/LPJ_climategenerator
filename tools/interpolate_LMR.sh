#!/bin/bash
# interpolate anomalies to 0.5 degree

INDIR=/group/esd_kaplan/datasets/ModelE/climategenerator_data/lmr_giss_climon/anomalies_1961-1990/nativeres/
OUTDIR=/group/esd_kaplan/datasets/ModelE/climategenerator_data/lmr_giss_climon/anomalies_1961-1990/30m/

fname=clt_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

fname=dtdiurn_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

fname=preacc_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

fname=sfcWind_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

fname=lightning_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

fname=tas_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile

fname=wetdays_Amonsynth_GISS-E2-R_anomalies.nc
original=${INDIR}${fname}
outfile=${OUTDIR}${fname}
cdo -s -f nc4 -P 8 remapbic,global_0.5 $original $outfile
