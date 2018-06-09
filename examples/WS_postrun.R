#!/usr/bin/Rscript

# Various additional post-run calculations.

options(stringsAsFactors = F)
library(WaveletScreaming)

workdir <- "/mnt/HUNT/other/wt/"

pheno <- read.table(paste0(workdir, "bw/bw_moms.txt"), h=T)
Yall <- pheno$VEKT
Confounderall <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/bw_moms_aa46.txt"), h=T)
Yall46 <- pheno$VEKT
Confounderall46 <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/bw_fets.txt"), h=T)
Yfets <- pheno$VEKT
Confounderfets <- as.matrix(pheno[,3:ncol(pheno)])

pheno <- read.table(paste0(workdir, "bw/bw_fetsx.txt"), h=T)
YfetsX <- pheno$VEKT
ConfounderfetsX <- as.matrix(pheno[,3:ncol(pheno)])


# LAMBDA VALUES
Lmoms = get_lambda1(Yall, Confounderall, 0.2)
Lmoms46 = get_lambda1(Yall46, Confounderall46, 0.2)
Lfets = get_lambda1(Yfets, Confounderfets, 0.2)
LfetsX = get_lambda1(YfetsX, ConfounderfetsX, 0.2)

save(list(Lmoms, Lmoms46, Lfets, LfetsX), file=paste0(workdir, "results/lambdas_bw.RData"))
