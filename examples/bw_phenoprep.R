#!/usr/bin/Rscript

############### COMMON INIT #####################

options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)

harvdir = "~/data/geno/harvest-aux/"
outdir = "~/Documents/benchmarking/wt/realpheno_bw/"

flags = read.table(paste0(harvdir, "harvest-flag-list.txt"), h=T)
link = read.table("~/data/mobaqs/p1724/harvest_linkage.csv", sep=";", h=T)
linkm = coremoms = corekids = corepaired = data.frame()


# Usage: attachPheno(mfr) or attachPheno(q1)
attachPheno = function(phenodf){
	linkm <<- inner_join(link, phenodf, by="PREG_ID_1724")
	print(sprintf("%i individuals with phenotype info found", nrow(linkm)))
}

# Usage: getCore(linkm)
getCore = function(pheno){
	coremoms <<- filter(flags, genotypesOK, phenotypesOK, coreOK) %>%
		semi_join(pheno, ., by=c("SentrixID_1"="IID")) %>%
		filter(Role=="Mother")
	# 7065, incl. multiple pregs
	print(sprintf("%i core mothers found", nrow(coremoms)))
	
	corekids <<- filter(flags, genotypesOK, phenotypesOK, coreOK) %>%
		semi_join(pheno, ., by=c("SentrixID_1"="IID")) %>%
		filter(Role=="Child")
	print(sprintf("%i core kids found", nrow(corekids)))
}

# removes repeated pregnancies,
# keeping the one included in the other core set.
# Usage: removeRepeated(coremoms, corekids) to clean mothers
# 		removeRepeated(corekids, coremoms) to clean kids
#			- obviously not needed, because kids don't have multiple MFR rows
removeRepeated = function(dftoclean, basedon){
	corepaired = group_by(dftoclean, SentrixID_1) %>%
		mutate(hasPair = PREG_ID_1724 %in% basedon$PREG_ID_1724) %>%
		top_n(1, hasPair) %>%
		filter(rank(PREG_ID_1724)==1)
	print(sprintf("after removing repeated pregnancies, %i remain", nrow(corepaired)))
	return(corepaired)
}

# Usage: attachCovariates AA87 moms 6
attachCovariates = function(corepaired, dataset, numPCs){
	## add batch covariate
	corepaired = select(flags, one_of(c("IID", "BATCH"))) %>%
		left_join(corepaired, ., by=c("SentrixID_1" = "IID")) %>%
		mutate(BATCH = as.numeric(BATCH=="M24"))
	
	# read in correct ibd-exclusion list and pca-covar file
	if(dataset=="moms"){
		ibd = read.table(paste0(harvdir, "ibd_pihat_exclude_mothers"))
		pcs = read.table(paste0(harvdir, "plink_covar_mothers"))
	} else if (dataset=="fets") {
		ibd = read.table(paste0(harvdir, "ibd_pihat_exclude_offspring"))
		pcs = read.table(paste0(harvdir, "plink_covar_offspring"))
	} else {
		return()
	}
	
	colnames(ibd) = c("FID1", "IID1", "IID2", "PIHAT")
	colnames(pcs)[1:2] = c("FID", "SentrixID_1")
	
	corepaired = inner_join(corepaired, pcs[,2:(numPCs+2)], by="SentrixID_1") %>%
		anti_join(ibd, by=c("SentrixID_1"="IID1"))
	print(sprintf("after attaching covariates, %i remain", nrow(corepaired)))
	return(corepaired)
}


# For SNPTEST
# Usage: makeOutputs moms_height.txt AA87 moms
makeOutputs = function(corepaired, phenofile, dataset){
	if(dataset=="moms"){
		samplelist = paste0(harvdir, "allmoms.txt")
	} else if (dataset=="fets"){
		samplelist = paste0(harvdir, "allfets.txt")
		samplelistx = paste0(harvdir, "allfetsx.txt")
	} else {
		return()
	}
	
	samplelist = read.table(samplelist, h=F)
	out = left_join(samplelist, corepaired, by=c("V1"="SentrixID_1"))[,-2]
	write.table(out, col.names=T, row.names=F, quote=F,
				file=paste0(outdir, phenofile, ".txt"))
	
	if(dataset=="fets"){
		samplelistx = read.table(samplelistx, h=F)
		outx = left_join(samplelistx, corepaired, by=c("V1"="SentrixID_1"))[,-2]
		write.table(outx, col.names=T, row.names=F, quote=F,
					file=paste0(outdir, phenofile, "x.txt"))
	}
}



################# PHENOTYPE-SPECIFIC PARTS ###################


### BW
m = read.table("~/data/mobaqs/p1724/harvest_mfr.csv", sep=";", h=T)
q1 = read.table("~/data/mobaqs/p1724/q1_pdb1724_v9.csv", sep=";", h=T)
dim(m)
dim(q1)
q1 = q1[,c("PREG_ID_1724", "AA46")]
q1 = filter(q1, AA46!=0) %>%
	mutate(AA46 = 2-AA46)

final = filter(m, FLERFODSEL==0,
			   DODKAT<6 | DODKAT>10,
			   is.na(DIABETES_MELLITUS),
			   C00_MALF_ALL==0,
			   !is.na(SVLEN_DG))

nrow(final) 

attachPheno(final)
getCore(linkm)

## MATERNAL
corepaired = removeRepeated(coremoms, corekids)
corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "VEKT", "SVLEN_DG", "KJONN")]

corepaired = attachCovariates(corepaired, "moms", 10)
makeOutputs(corepaired, "bw_moms", "moms") 

## FETAL
corepaired = removeRepeated(corekids, coremoms)
corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "VEKT", "SVLEN_DG", "KJONN")]

corepaired = attachCovariates(corepaired, "fets", 10)
makeOutputs(corepaired, "bw_fets", "fets") 


## WITH AA46 covariate
finalq = inner_join(final, q1, by="PREG_ID_1724")

attachPheno(finalq)
getCore(linkm)
## MATERNAL
corepaired = removeRepeated(coremoms, corekids)
corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "VEKT", "SVLEN_DG", "KJONN", "AA46")]

corepaired = attachCovariates(corepaired, "moms", 10)
makeOutputs(corepaired, "bw_moms_aa46", "moms") 

## FETAL
corepaired = removeRepeated(corekids, coremoms)
corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "VEKT", "SVLEN_DG", "KJONN", "AA46")]

corepaired = attachCovariates(corepaired, "fets", 10)
makeOutputs(corepaired, "bw_fets_aa46", "fets") 

