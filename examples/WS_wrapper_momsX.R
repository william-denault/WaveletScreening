#!/usr/bin/Rscript

# Main analysis wrapper.

# ARGUMENTS
# 1. current batch, indexed from 0
# 2. slices to process in each batch, e.g. 50
# Input files are specified below.

# USAGE
# Can be launched manually on a single slice as:
#    ./WS_wrapper.R 125 1
# as a single-thread process over all slices as:
#    ./WS_wrapper.R 0 `wc -l slicedef.txt`
# or in parallel as:
#    seq 0 99 | parallel -j 10 ./WS_wrapper.R {} 50
# (so total (batch+1)*batchSize = 100*50 = 5000 slices would be processed)

options(stringsAsFactors = F)

library(EbayesThresh)
library(readr)
library(wavelets)
library(wavethresh)
library(data.table)
library(corpcor)
library(rARPACK)
library(WaveletScreaming)
library(RhpcBLASctl)

blas_set_num_threads(1)

args = commandArgs(T)
batch = as.numeric(args[1])
batchSize = as.numeric(args[2])

#########################
### SET INPUT HERE

workdir <- "/media/local-disk2/jjuod/other/wt/"
resdir <- "results/momsX/"

# read phenotype and covariate data
pheno <- read.table(paste0(workdir, "bw/bw_moms.txt"), h=T)
Yall <- pheno$VEKT
Confounderall <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/bw_moms_aa46.txt"), h=T)
Yall46 <- pheno$VEKT
Confounderall46 <- as.matrix(pheno[,3:ncol(pheno)])

# identify slices to read
slices = read.table(paste0(workdir, "bw/slicedef_moms_X.txt"), h=F, sep=" ")
slicestem = paste0(workdir, "slicesbw/moms_")

### END OF INPUT SETTINGS
#########################


print(paste("Specified model: Y ~", paste(colnames(Confounderall), collapse=" + ")))
print(paste("Specified model: Y ~", paste(colnames(Confounderall46), collapse=" + ")))
print(sprintf("%i slices read from definition file", nrow(slices)))

start = batch*batchSize+1
end = min(nrow(slices), (batch+1)*batchSize)
if(start>end){
	stop("ERROR: starting position outside the range of provided slice definitions")
} else {
	print(sprintf("processing slices between %i and %i", start, end))
}


res_coefc = res_coefd  = res_coefc_aa46 = res_coefd_aa46 = data.frame()

# loop over slices
for(i in start:end){
  ptm <- proc.time()
  print(sprintf("working on locus %i of %i for chromosome %s", i, nrow(slices), slices[i,1]))

  #Importing
  my_slice  <- fread(paste0("gzip -dc ", slicestem,
  					  slices[i,1], "_", slices[i,2], "-", slices[i,3], ".vcf.gz"))
  bp <- unlist(my_slice[,2], use.names=F)
  my_slice <- as.matrix(my_slice[,3:(ncol(my_slice)-1)])

  print("locus imported successfully")
  print(proc.time()-ptm)


  ## SETTINGS 1
  # Run the WaveletScreaming function
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screaming(Y=Yall,loci=my_slice,bp=bp,confounder=Confounderall,
  						  lev_res=9, sigma_b=0.2, coeftype="c", para=FALSE)

  # Save the results into a result matrix after each run.
  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefc <- rbind(res_coefc, t(res))
  print("analysis with coeftype C, no AA46, completed successfully")
  print(proc.time()-ptm)

  ## SETTINGS 2
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screaming(Y=Yall,loci=my_slice,bp=bp,confounder=Confounderall,
  						  lev_res=9, sigma_b=0.2, coeftype="d", para=FALSE)

  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefd <- rbind(res_coefd, t(res))
  print("analysis with coeftype D, no AA46, completed successfully")
  print(proc.time()-ptm)

  ## SETTINGS 3
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screaming(Y=Yall46,loci=my_slice,bp=bp,confounder=Confounderall46,
  						  lev_res=9, sigma_b=0.2, coeftype="c", para=FALSE)

  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefc_aa46 <- rbind(res_coefc_aa46, t(res))
  print("analysis with coeftype C, with AA46, completed successfully")
  print(proc.time()-ptm)

  ## SETTINGS 4
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screaming(Y=Yall46,loci=my_slice,bp=bp,confounder=Confounderall46,
  						  lev_res=9, sigma_b=0.2, coeftype="d", para=FALSE)

  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefd_aa46 <- rbind(res_coefd_aa46, t(res))
  print("analysis with coeftype D, with AA46, completed successfully")
  print(proc.time()-ptm)
}

write.table(res_coefc,
			paste0(workdir, resdir, "res_coefc_", start, "-", end, ".txt"),
			row.names=F, col.names=T, quote=F)
write.table(res_coefd,
			paste0(workdir, resdir, "res_coefd_", start, "-", end, ".txt"),
			row.names=F, col.names=T, quote=F)
write.table(res_coefc_aa46,
			paste0(workdir, resdir, "res_coefc_aa46_", start, "-", end, ".txt"),
			row.names=F, col.names=T, quote=F)
write.table(res_coefd_aa46,
			paste0(workdir, resdir, "res_coefd_aa46_", start, "-", end, ".txt"),
			row.names=F, col.names=T, quote=F)
