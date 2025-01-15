#!/bin/bash
#
# To run this script in "batch mode":
#
#   ./run_benchmarks.sh >run_benchmarks.log 2>&1 &
#

set -e  # exit immediately if a simple command exits with a non-zero status

## Try to use R to obtain the path to Rscript (requires R in the PATH).
RSCRIPT=`R -s --vanilla -e 'cat(file.path(R.home("bin"), "Rscript"))'`

## Manually set RSCRIPT here if R is not in the PATH.
#RSCRIPT=path/to/Rscript

NORMALIZE_AND_PCA_R=`$RSCRIPT -e 'suppressPackageStartupMessages(library(HDF5Array)); cat(system.file(package="HDF5Array", "scripts", "normalize_and_PCA.R", mustWork=TRUE))'`

date >> normalize_and_PCA_timings.txt
$RSCRIPT $NORMALIZE_AND_PCA_R sparse  12500 250 100
$RSCRIPT $NORMALIZE_AND_PCA_R dense   12500 100 100
$RSCRIPT $NORMALIZE_AND_PCA_R sparse  25000 250 100
$RSCRIPT $NORMALIZE_AND_PCA_R dense   25000 100 100
$RSCRIPT $NORMALIZE_AND_PCA_R sparse  50000 250 100
$RSCRIPT $NORMALIZE_AND_PCA_R dense   50000 100 100
$RSCRIPT $NORMALIZE_AND_PCA_R sparse 100000 250 100
$RSCRIPT $NORMALIZE_AND_PCA_R dense  100000 100 100
$RSCRIPT $NORMALIZE_AND_PCA_R sparse 200000 250 100
$RSCRIPT $NORMALIZE_AND_PCA_R dense  200000 100 100
date >> normalize_and_PCA_timings.txt

