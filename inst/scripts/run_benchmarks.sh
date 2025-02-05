#!/bin/bash
#
# To run this script in "batch mode":
#
#   cd path/to/HDF5Array/inst/scripts/timings_dbs/<machine-name>
#   time ../../run_benchmarks.sh >run_benchmarks.log 2>&1 &
#
# Takes about 32 hours to complete on xps15 with HDF5Array 1.35.12!
# Should take between 20 and 50 hours to complete, depending on the machine.
#
set -e  # exit immediately if a simple command exits with a non-zero status

## Try to use R to obtain the path to Rscript (requires R in the PATH).
RSCRIPT=`R -s --vanilla -e 'cat(file.path(R.home("bin"), "Rscript"))'`

## Manually set RSCRIPT here if R is not in the PATH.
#RSCRIPT=path/to/Rscript

NORMALIZE_AND_PCA_R=`$RSCRIPT -e 'suppressPackageStartupMessages(library(HDF5Array)); cat(system.file(package="HDF5Array", "scripts", "normalize_and_PCA.R", mustWork=TRUE))'`

normalize_and_PCA()
{
	ncells="$1"
	num_var_genes="$2"
	format="$3"
	norm_block_size="$4"
	realize_block_size="$4"
	pca_block_size="$4"
	$RSCRIPT $NORMALIZE_AND_PCA_R "$ncells" "$num_var_genes" "$format" "$norm_block_size" "$realize_block_size" "$pca_block_size"
}

echo "Starting run_benchmarks.sh on `date`."
echo ""

for ncells in 12500 25000 50000 100000 200000; do
	for num_var_genes in 1000 2000; do
		for format in s D Ds; do
			for block_size in 16 40 100 250; do
				normalize_and_PCA "$ncells" "$num_var_genes" "$format" "$block_size"
			done
		done
	done
done

echo "Completed run_benchmarks.sh on `date`."
echo ""

dest_file="timings-`date +\%Y\%m\%d`.dcf"
mv timings.dcf $dest_file
echo "See timings in '$dest_file'."
echo ""

