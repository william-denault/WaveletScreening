#!/usr/bin/Rscript

# Various additional post-run calculations.

options(stringsAsFactors = F)
library(WaveletScreaming)

workdir <- "/mnt/HUNT/other/wt/"

pheno <- read.table(paste0(workdir, "bw/bw_moms.txt"), h=T)
Yall <- pheno$VEKT
Confounderall <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/moms_aa46.txt"), h=T)
Yall46 <- pheno$AA46
Confounderall46 <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/bw_fets.txt"), h=T)
Yfets <- pheno$VEKT
Confounderfets <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/bw_fetsx.txt"), h=T)
YfetsX <- pheno$VEKT
ConfounderfetsX <- as.matrix(pheno[,3:ncol(pheno)])


# LAMBDA VALUES
# standard method:
Lmoms = get_lambda1(Yall, Confounderall, 0.2)
Lmoms46 = get_lambda1(Yall46, Confounderall46, 0.2)
Lfets = get_lambda1(Yfets, Confounderfets, 0.2)
LfetsX = get_lambda1(YfetsX, ConfounderfetsX, 0.2)

output = c("Lmoms"=Lmoms, "Lmoms46"=Lmoms46, "Lfets"=Lfets, "LfetsX"=LfetsX)
print(output, digits = 17)

save("output", file=paste0(workdir, "results/lambdas_bw.RData"))

# fast method:
Lmoms2 = get_lambda2(Yall, Confounderall, 0.2)
Lmoms462 = get_lambda2(Yall46, Confounderall46, 0.2)
Lfets2 = get_lambda2(Yfets, Confounderfets, 0.2)
LfetsX2 = get_lambda2(YfetsX, ConfounderfetsX, 0.2)

output2 = c("Lmoms"=Lmoms2, "Lmoms46"=Lmoms462, "Lfets"=Lfets2, "LfetsX"=LfetsX2)
print(output2, digits = 17)


# FUNCTION FOR WRITING OUTPUT
formatter = function(df){
  left = df[,1:3]
  lambda = format(df[,4], digits=10)
  right = format(df[,-c(1:4)], digits=4)
  return(cbind(left, lambda, right))
}


# GATHER MATERNAL RESULTS
workdir <- "/mnt/HUNT/other/wt/"
resdir <- "results/momsAUTO/"
slices = read.table(paste0(workdir, "bw/slicedef_moms_autosomes.txt"))

goodnames = names(res) # had to recreate names for now because they were lost in rbind

batchSize = 50
res_coefc_full = res_coefd_full = data.frame()
for(b in 0:100){
  start = b*batchSize+1
  end = min(nrow(slices), (b+1)*batchSize)
  print(sprintf("reading batch %i out of 100", b))

  res_coefc = read.table(paste0(workdir, resdir, "res_coefc_", start, "-", end, ".txt"), h=T)
  colnames(res_coefc) = goodnames
  res_coefc_full = rbind(res_coefc_full, res_coefc)

  res_coefd = read.table(paste0(workdir, resdir, "res_coefd_", start, "-", end, ".txt"), h=T)
  colnames(res_coefd) = goodnames
  res_coefd_full = rbind(res_coefd_full, res_coefd)
}

resdir <- "results/momsX/"
slices = read.table(paste0(workdir, "bw/slicedef_moms_X.txt"))

batchSize = 15
for(b in 0:12){
  start = b*batchSize+1
  end = min(nrow(slices), (b+1)*batchSize)
  print(sprintf("reading batch %i out of 12", b))

  res_coefc = read.table(paste0(workdir, resdir, "res_coefc_", start, "-", end, ".txt"), h=T)
  colnames(res_coefc) = goodnames
  res_coefc_full = rbind(res_coefc_full, res_coefc)

  res_coefd = read.table(paste0(workdir, resdir, "res_coefd_", start, "-", end, ".txt"), h=T)
  colnames(res_coefd) = goodnames
  res_coefd_full = rbind(res_coefd_full, res_coefd)
}

# 5198 rows each
nrow(res_coefc_full)
nrow(res_coefd_full)

write.table(formatter(res_coefc_full),
            gzfile(paste0(workdir, "bw/full_res_coefc.txt.gz"), "w"),
            col.names=T, row.names=F, quote=F)
write.table(formatter(res_coefd_full),
            gzfile(paste0(workdir, "bw/full_res_coefd.txt.gz"), "w"),
            col.names=T, row.names=F, quote=F)
closeAllConnections()


# GATHER FETAL RESULTS
workdir <- "/mnt/HUNT/other/wt/"
resdir <- "results/fetsAUTO/"
slices = read.table(paste0(workdir, "bw/slicedef_fets_autosomes.txt"))

batchSize = 50
res_coefc_full = res_coefd_full = data.frame()
for(b in 0:100){
  start = b*batchSize+1
  end = min(nrow(slices), (b+1)*batchSize)
  print(sprintf("reading batch %i out of 100", b))

  res_coefc = read.table(paste0(workdir, resdir, "res_coefc_", start, "-", end, ".txt"), h=T)
  res_coefc_full = rbind(res_coefc_full, res_coefc)

  res_coefd = read.table(paste0(workdir, resdir, "res_coefd_", start, "-", end, ".txt"), h=T)
  res_coefd_full = rbind(res_coefd_full, res_coefd)
}

resdir <- "results/fetsX/"
slices = read.table(paste0(workdir, "bw/slicedef_fets_X.txt"))

batchSize = 10
for(b in 0:19){
  start = b*batchSize+1
  end = min(nrow(slices), (b+1)*batchSize)
  print(sprintf("reading batch %i out of 19", b))

  res_coefc = read.table(paste0(workdir, resdir, "res_coefc_", start, "-", end, ".txt"), h=T)
  res_coefc_full = rbind(res_coefc_full, res_coefc)

  res_coefd = read.table(paste0(workdir, resdir, "res_coefd_", start, "-", end, ".txt"), h=T)
  res_coefd_full = rbind(res_coefd_full, res_coefd)
}

# 5202 rows each
nrow(res_coefc_full)
nrow(res_coefd_full)

write.table(formatter(res_coefc_full),
            gzfile(paste0(workdir, "bw/full_fets_coefc.txt.gz"), "w"),
            col.names=T, row.names=F, quote=F)
write.table(formatter(res_coefd_full),
            gzfile(paste0(workdir, "bw/full_fets_coefd.txt.gz"), "w"),
            col.names=T, row.names=F, quote=F)
closeAllConnections()


# GATHER MATERNAL AA46 RESULTS
workdir <- "/mnt/HUNT/other/wt/"
resdir <- "results/momsAA46/"
slices = read.table(paste0(workdir, "bw/slicedef_moms_autosomes.txt"))

batchSize = 50
res_coefc_full = res_coefd_full = data.frame()
for(b in 0:100){
  start = b*batchSize+1
  end = min(nrow(slices), (b+1)*batchSize)
  print(sprintf("reading batch %i out of 100", b))

  res_coefc = read.table(paste0(workdir, resdir, "res_coefc_aa46_", start, "-", end, ".txt"), h=T)
  res_coefc_full = rbind(res_coefc_full, res_coefc)

  res_coefd = read.table(paste0(workdir, resdir, "res_coefd_aa46_", start, "-", end, ".txt"), h=T)
  res_coefd_full = rbind(res_coefd_full, res_coefd)
}

resdir <- "results/momsAA46X/"
slices = read.table(paste0(workdir, "bw/slicedef_moms_X.txt"))

batchSize = 10
for(b in 0:19){
  start = b*batchSize+1
  end = min(nrow(slices), (b+1)*batchSize)
  print(sprintf("reading batch %i out of 19", b))

  res_coefc = read.table(paste0(workdir, resdir, "res_coefc_aa46_", start, "-", end, ".txt"), h=T)
  res_coefc_full = rbind(res_coefc_full, res_coefc)

  res_coefd = read.table(paste0(workdir, resdir, "res_coefd_aa46_", start, "-", end, ".txt"), h=T)
  res_coefd_full = rbind(res_coefd_full, res_coefd)
}

# 5198 rows each
nrow(res_coefc_full)
nrow(res_coefd_full)

write.table(formatter(res_coefc_full),
            gzfile(paste0(workdir, "bw/full_res_coefca.txt.gz"), "w"),
            col.names=T, row.names=F, quote=F)
write.table(formatter(res_coefd_full),
            gzfile(paste0(workdir, "bw/full_res_coefda.txt.gz"), "w"),
            col.names=T, row.names=F, quote=F)
closeAllConnections()
