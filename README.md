
## Wavelet Screening: a tutorial





# Installation #
---

Using the 'devtools' package:

```{r , eval=FALSE}
install.packages( "devtools" )
library( devtools )
install_github( "william-denault/WaveletScreening" )
```

Soon on CRAN



# Usage #

See help from Wavelet_screening

### Generating simulated data:
Here, a dummy example.
```{r , echo=FALSE}
#########################################
#Generate a randomly sampled locus size of 1Mb
#########################################

#5000 Randomly chosen bp
my_bp <- sort(sample(1:1000000, size=5000,replace = FALSE))
#############################
#Three different bump signals
#############################
my_functions <-data.frame(f0 = c(rep(0,400000),rep(0,200000),rep(0,400000)),
                         f1 = c(rep(0,400000),rep(1,200000),rep(0,400000)) ,
                         f2=c(rep(0,400000),rep(2,200000),rep(0,400000)))



###########################
#Minor allele frequency of 30%
###########################
MAF=0.3
sampl_schem <- c((1-MAF)^2,2*MAF*(1-MAF),MAF^2)
#######################################
#Sampling assuming Hardy-Weinberg equilibrium
#######################################
#Assigning class

#sample size =4000
n_size=4000
type_fn <-sample(0:2,replace = TRUE,size=n_size,  prob=  sampl_schem  )


genotype <-  matrix(my_functions[my_bp,2 ], ncol=1 ) %*%t(matrix(type_fn,ncol=1))
#dim(genotype)= nSNP, nind

###############################################################
#Generate a phenotype with variance explained by genotype  of 0.5%
###############################################################
varexp=0.005
var_noise <- (1-varexp)*var(sample(0:2,replace = TRUE,size=10000,
                                  prob=sampl_schem ))/varexp
Y <-  rnorm(n=n_size,sd=sqrt(var_noise)) +type_fn

```

### Performing the analysis


```{r , echo=FALSE}
res <- Wavelet_screening( Y,loci=genotype,bp=my_bp,
                         lev_res=6,sigma_b = 0.2)

```



### Assessing the significance

The test statistics of the Wavelet Screening method relies on a combination of two statistics. The first one L_h is the basis for the test, and the second is used as a penalization.

(it takes 5-7 minutes)
```{r , echo=FALSE}
Sim <- Simu_null(Y,lev_res = 6,sigma_b = 0.2,size=100000)

```
####Calibration of the hyperparameter
```{r , echo=FALSE}
lambda <- Search_lambda(Sim,plot=TRUE)

```



#### P-value 

We use the simulated observations of the null distribution and the penalty parameter to fit a normal distribution.  We then use this fitted distribution to estimate the associated quantile of the observation and  the p-value.
```{r , echo=FALSE}
Th <- Sim[,c("L_h")]+lambda*Sim[,c("min_ph_pv")]
mu <- median(Th,na.rm = TRUE)
sdv <- mad(Th,na.rm = TRUE)
####################################
#Test value of the locus to be tested
####################################
th <-  res[c("L_h")]+lambda*res["min_ph_pv"]
#######
#P-value
#######
pval <- 1-pnorm(th,mean=muv,sd=sdv)
pval
df <- data.frame(Th = Th,type = factor( c(rep("Null",length(Th)))) )
ggplot(df,aes(Th,fill=type))+
  xlim(c(min(c(Th,th)),max(Th,th)))+
  geom_density()+
  guides(fill=FALSE)+
  geom_point(aes(x=th, y=0), colour="red")+theme(legend.position="none")+
  geom_text(label="Value of the test statistics", x=th, y=0.001)+
  geom_text(label="Null distribution", x=mean(Th), y=0.001)+
  theme_bw()
```

### Visualization

You can visualize your results using the function *plot_WS*.
It represents the Betas for the different scales of the wavelet decomposition. The size and the darkness of the points that represent the Betas are scaled by the value of the Betas. If a Beta is not thresholded, then the region that represents the Beta is highlighted in order to give an overview of the size and the origin of the genetic signal.


```{r , echo=FALSE}
bp <- c(min(my_bp),max(my_bp))
plot_WS(res=res,bp=bp,lev_res=6)
```




# Genome-wide screening

You can find in the folder *examples*, our script for the full-genome screening.

## Slicing of the genome into regions

As the Wavelet Screening is made to analyze regions of the genome, the first step is to divide your data in *screenable regions*. The function *slice_definition* is meant for that. The user has to specify the size of the region to analyze (Loci_size) *AND* the maximum gap between two SNPs (thresh option) (in our setup, Loci_size=1e6, thresh=1e4).

Below is a script that slices the vcf file and stores the slices in a specified directory. These slices will be analyzed in the next step.

```{r , echo=FALSE}
library(WaveletScreening)

setwd("yourworking path")
outdir = "your output path"
bcftools = "whereisyour/bcftools"

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

```




## Phenotype preparation for analysis

IMPORTANT: The Wavelet Screening function assumes that your covariate, your confounder matrix, and your genotype matrix are sorted in the same order.



## Genome-wide run
Below is a code for dealing with genome-wide data.
In the *example* folder, you will find detailed code for the Wavelet Screening package that you can use via a terminal in combination with *bcftools*.
We choose to store the slices in a given directory instead of working directly in the memory. This might not be optimal. However, it is handier when you want to examine various models, different confounding structures, and different phenotypes.
The provided scripts are already parallelized.


```{r , echo=FALSE}
#!/usr/bin/Rscript

# Main analysis wrapper.

# ARGUMENTS
# 1. Current batch, indexed from 0
# 2. Slices to process in each batch, e.g., 50
# Input files are specified below.

# USAGE
# Can be launched manually on a single slice as follows:
#    ./WS_wrapper.R 125 1
# As a single-thread process over all slices as:
#    ./WS_wrapper.R 0 `wc -l slicedef.txt`
# Or in parallel as:
#    seq 0 99 | parallel -j 10 ./WS_wrapper.R {} 50
# (So, the total (batch+1)*batchSize = 100*50 = 5000 slices would be processed)

options(stringsAsFactors = F)

library(EbayesThresh)
library(readr)
library(wavelets)
library(wavethresh)
library(data.table)
library(corpcor)
library(rARPACK)
library(WaveletScreening)
library(RhpcBLASctl)

blas_set_num_threads(1)

args = commandArgs(T)
batch = as.numeric(args[1])
batchSize = as.numeric(args[2])

#########################
### Set input here

workdir <- "/media/local-disk2/jjuod/other/wt/"

# Read phenotypes and covariate data:
pheno <- read.table(paste0(workdir, "bw/bw_moms.txt"), h=T)
Yall <- pheno$VEKT #Our phenoytpe

#Confounding matrix:
Confounderall <- as.matrix(pheno[,3:ncol(pheno)])#Where we put our  


# Identify slices to read:
slices = read.table(paste0(workdir, "bw/slicedef_moms_autosomes.txt"), h=F, sep=" ")
slicestem = paste0(workdir, "slicesbw/moms_")

### End of input
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


res_coefc = res_coefd  = data.frame()

# Loop over slices
for(i in start:end){
  ptm <- proc.time()
  print(sprintf("working on locus %i of %i for chromosome %s", i, nrow(slices), slices[i,1]))

  #Importing slices
  my_slice  <- fread(paste0("gzip -dc ", slicestem,
         slices[i,1], "_", slices[i,2], "-", slices[i,3], ".vcf.gz"))
  bp <- unlist(my_slice[,2], use.names=F)
  my_slice <- as.matrix(my_slice[,3:(ncol(my_slice)-1)])

  print("Locus imported successfully")
  print(proc.time()-ptm)


  ## Coeftype D
  # Run the Wavelet Screening function
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screening(Y=Yall,loci=my_slice,bp=bp,confounder=Confounderall,
          lev_res=9, sigma_b=0.2, coeftype="d", para=FALSE)

  # Save the results in a result matrix after each run.
  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefc <- rbind(res_coefc, t(res))
  print("Analysis with coeftype C, no AA46, completed successfully.")
  print(proc.time()-ptm)

  ##  Coeftype C
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screening(Y=Yall,loci=my_slice,bp=bp,confounder=Confounderall,
          lev_res=9, sigma_b=0.2, coeftype="c", para=FALSE)
  # Save the results in a result matrix after each run.
  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefd <- rbind(res_coefd, t(res))
  print("Analysis with coeftype D, no AA46, completed successfully.")
  print(proc.time()-ptm)


}

write.table(res_coefc,
   paste0(workdir, "results/res_coefc_", start, "-", end, ".txt"),
   row.names=F, col.names=T, quote=F)
write.table(res_coefd,
   paste0(workdir, "results/res_coefd_", start, "-", end, ".txt"),
   row.names=F, col.names=T, quote=F)

```
