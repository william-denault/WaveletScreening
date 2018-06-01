#!/bin/bash

# Wrapper for parallel calling of make_slice.sh.
# Requires bcftools and GNU parallel.
# Input files: slice definitions and imputed VCFs.

# (Can also be moved into bw_defineslices.R)

cat slicedef_moms_{1..22}.txt > slicedef_moms_autosomes.txt
cat slicedef_fets_{1..22}.txt > slicedef_fets_autosomes.txt

# starts 25 jobs of make_slice.sh for each definition file
# passes a row from slicedef.txt as args 1-3
# and adds M or F as arg 4 to the job
parallel -j 25 --progress -a slicedef_moms_autosomes.txt --colsep ' ' ./make_slice.sh {} M
parallel -j 25 --progress -a slicedef_fets_autosomes.txt --colsep ' ' ./make_slice.sh {} F
parallel -j 25 --progress -a slicedef_moms_X.txt --colsep ' ' ./make_slice.sh {} M
parallel -j 25 --progress -a slicedef_fets_X.txt --colsep ' ' ./make_slice.sh {} F
