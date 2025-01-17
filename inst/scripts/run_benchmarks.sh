#!/bin/bash
#
# To run this script in "batch mode":
#
#   cd path/to/HDF5Array/inst/scripts/timings_db/<machine-name>
#   time ../../run_benchmarks.sh >run_benchmarks.log 2>&1 &
#
# It takes between 5 and 10 hours to complete, depending on the machine!

set -e  # exit immediately if a simple command exits with a non-zero status

## Try to use R to obtain the path to Rscript (requires R in the PATH).
RSCRIPT=`R -s --vanilla -e 'cat(file.path(R.home("bin"), "Rscript"))'`

## Manually set RSCRIPT here if R is not in the PATH.
#RSCRIPT=path/to/Rscript

NORMALIZE_AND_PCA_R=`$RSCRIPT -e 'suppressPackageStartupMessages(library(HDF5Array)); cat(system.file(package="HDF5Array", "scripts", "normalize_and_PCA.R", mustWork=TRUE))'`

echo "Starting run_benchmarks.sh on `date`."
echo ""

# --------------------------- ncells format norm_block_size pca_block_size
$RSCRIPT $NORMALIZE_AND_PCA_R  12500 sparse              40             40
$RSCRIPT $NORMALIZE_AND_PCA_R  12500 sparse             100            100
$RSCRIPT $NORMALIZE_AND_PCA_R  12500 sparse             250            250
$RSCRIPT $NORMALIZE_AND_PCA_R  12500 dense               40             40
$RSCRIPT $NORMALIZE_AND_PCA_R  12500 dense              100            100
$RSCRIPT $NORMALIZE_AND_PCA_R  12500 dense              250            250
# --------------------------- ncells format norm_block_size pca_block_size
$RSCRIPT $NORMALIZE_AND_PCA_R  25000 sparse              40             40
$RSCRIPT $NORMALIZE_AND_PCA_R  25000 sparse             100            100
$RSCRIPT $NORMALIZE_AND_PCA_R  25000 sparse             250            250
$RSCRIPT $NORMALIZE_AND_PCA_R  25000 dense               40             40
$RSCRIPT $NORMALIZE_AND_PCA_R  25000 dense              100            100
$RSCRIPT $NORMALIZE_AND_PCA_R  25000 dense              250            250
# --------------------------- ncells format norm_block_size pca_block_size
$RSCRIPT $NORMALIZE_AND_PCA_R  50000 sparse              40             40
$RSCRIPT $NORMALIZE_AND_PCA_R  50000 sparse             100            100
$RSCRIPT $NORMALIZE_AND_PCA_R  50000 sparse             250            250
$RSCRIPT $NORMALIZE_AND_PCA_R  50000 dense               40             40
$RSCRIPT $NORMALIZE_AND_PCA_R  50000 dense              100            100
$RSCRIPT $NORMALIZE_AND_PCA_R  50000 dense              250            250
# --------------------------- ncells format norm_block_size pca_block_size
$RSCRIPT $NORMALIZE_AND_PCA_R 100000 sparse              40             40
$RSCRIPT $NORMALIZE_AND_PCA_R 100000 sparse             100            100
$RSCRIPT $NORMALIZE_AND_PCA_R 100000 sparse             250            250
$RSCRIPT $NORMALIZE_AND_PCA_R 100000 dense               40             40
$RSCRIPT $NORMALIZE_AND_PCA_R 100000 dense              100            100
$RSCRIPT $NORMALIZE_AND_PCA_R 100000 dense              250            250
# --------------------------- ncells format norm_block_size pca_block_size
$RSCRIPT $NORMALIZE_AND_PCA_R 200000 sparse              40             40
$RSCRIPT $NORMALIZE_AND_PCA_R 200000 sparse             100            100
$RSCRIPT $NORMALIZE_AND_PCA_R 200000 sparse             250            250
$RSCRIPT $NORMALIZE_AND_PCA_R 200000 dense               40             40
$RSCRIPT $NORMALIZE_AND_PCA_R 200000 dense              100            100
$RSCRIPT $NORMALIZE_AND_PCA_R 200000 dense              250            250

echo "Completed run_benchmarks.sh on `date`."
echo ""

dest_file="timings-`date +\%Y\%m\%d`.dcf"
mv timings.dcf $dest_file
echo "See timings in '$dest_file'."
echo ""

