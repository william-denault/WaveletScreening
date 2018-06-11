#!/usr/bin/Rscript

# Extracts observed SNP positions and defines slices.
# Requires bcftools.
# Input: full, imputed VCFs.
# Arguments - chromosomes to be processed.

# Usage: ./defineslices.R 1 2 3 22 X

options(stringsAsFactors = F)
library(WaveletScreaming)

setwd("/media/local-disk2/jjuod/merging/")
outdir = "/media/local-disk2/jjuod/other/wt/bw/"
bcftools = "/media/local-disk2/jjuod/erc-genotypes/bin/bcftools"

args = commandArgs(TRUE)

for(chr in args){
	print(chr)
	cmd = paste0(bcftools, " query -f '%POS\n' -e 'AC/AN<0.01 | AC/AN>0.99' moms_", chr, ".vcf.gz")
	bp = readLines(pipe(cmd))
	if(length(bp)==0) stop("ERROR: no positions read!")
	print(sprintf("for maternal chr %s, %i snps read", chr, length(bp)))
	bp = as.numeric(bp)
	print(head(bp))
	s = slice_definition(bp, 1e6, 1e4, chr)
	write.table(s, paste0(outdir, "slicedef_moms_", chr, ".txt"), quote=F, col.names=F, row.names=F)
	
	cmd = paste0(bcftools, " query -f '%POS\n' -e 'AC/AN<0.01 | AC/AN>0.99' fets_", chr, ".vcf.gz")
	bp = readLines(pipe(cmd))
	if(length(bp)==0) stop("ERROR: no positions read!")
	print(sprintf("for fetal chr %s, %i snps read", chr, length(bp)))
	bp = as.numeric(bp)
	print(head(bp))
	s = slice_definition(bp, 1e6, 1e4, chr)
	write.table(s, paste0(outdir, "slicedef_fets_", chr, ".txt"), quote=F, col.names=F, row.names=F)
}


