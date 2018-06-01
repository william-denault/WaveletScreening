#!/bin/bash

# Wrapper for parallel calling of make_slice.sh.
# Requires bcftools and GNU parallel.
# Input files: slice definitions and imputed VCFs.

# (Can also be moved into bw_defineslices.R)

for chr in {1..22} "X"
do
	# starts 25 jobs of make_slice.sh for this chromosome
	# passes a row from slicedef.txt as args 1-3
	# and adds M or F as arg 4 to the job
	parallel -j 25 --progress -a slicedef_moms_${chr}.txt --colsep ' ' ./make_slice.sh {} M
	parallel -j 25 --progress -a slicedef_fets_${chr}.txt --colsep ' ' ./make_slice.sh {} F
done
