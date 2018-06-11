#!/bin/bash

# Extracts one small VCF slice from a large imputed VCF.
# Requires bcftools.
# USAGE:
# ./make_slice.sh chromosome startPos endPos M/F
# (last argument selects source genome file)


set -e

BCFTOOLS=/media/local-disk2/jjuod/erc-genotypes/bin/bcftools
INDIR=/media/local-disk2/jjuod/merging/
OUTDIR=../slicebw/

echo "Running slice ${1}:${2}-${3}..."

if [ "$4" == "M" ]
then
	echo "using maternal genotypes..."
	$BCFTOOLS query -H -f '%CHROM\t%POS\t [ %DS \t]\n' \
		-r ${1}:${2}-${3} \
		${INDIR}/moms_${1}.vcf.gz | \
		gzip -c > ${OUTDIR}/moms_${1}_${2}-${3}.vcf.gz
else
	echo "using fetal genotypes..."
	$BCFTOOLS query -H -f '%CHROM\t%POS\t [ %DS \t]\n' \
		-r ${1}:${2}-${3} \
		${INDIR}/fets_${1}.vcf.gz | \
		gzip -c > ${OUTDIR}/fets_${1}_${2}-${3}.vcf.gz
fi
echo "Completed slice ${1}:${2}-${3}"
