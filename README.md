
## Wavelet Screaming: a tutorial





# Installation #
---

Using the 'devtools' package:

```{r , eval=FALSE}
install.packages( "devtools" )
library( devtools )
install_github( "william-denault/WaveletScreaming" )
```

Soon on CRAN



# Usage #

See help from Wavelet_screaming

### Generating fake data:

```{r , echo=FALSE}
#########################################
#Generate a randomly sample loci size=1Mb
#########################################

#5000 Randomly choosen bp
my_bp <- sort(sample(1:1000000, size=5000,replace = FALSE))
#############################
#Three different bump signals
#############################
my_functions <-data.frame(f0 = c(rep(0,400000),rep(0,200000),rep(0,400000)),
                         f1 = c(rep(0,400000),rep(1,200000),rep(0,400000)) ,
                         f2=c(rep(0,400000),rep(2,200000),rep(0,400000)))



###########################
#Minor allele frequency 30%
###########################
MAF=0.3
sampl_schem <- c((1-MAF)^2,2*MAF*(1-MAF),MAF^2)
#######################################
#sampling at Hardy Weinberg equilibrium
#######################################
#Assigning class

#sample size =4000
n_size=4000
type_fn <-sample(0:2,replace = TRUE,size=n_size,  prob=  sampl_schem  )


genotype <-  matrix(my_functions[my_bp,2 ], ncol=1 ) %*%t(matrix(type_fn,ncol=1))
#dim(genotype)= nSNP, nind

###############################################################
#Generate a phenotype with variance explained by genotype  0.5%
###############################################################
varexp=0.005
var_noise <- (1-varexp)*var(sample(0:2,replace = TRUE,size=10000,
                                  prob=sampl_schem ))/varexp
Y <-  rnorm(n=n_size,sd=sqrt(var_noise)) +type_fn

```

### Performing the analysis


You can specify the wavelet coefficients that you want to use. From our experiences (see our paper), the d decomposition is more robust to SNP coding misspecification. So we advise to use *d* in general. However *c* coefficients analysis are more powerfull if the coding is consistent with the direction of the effect (see our paper section simulation).


```{r , echo=FALSE}
res <- Wavelet_screaming( Y,loci=genotype,bp=my_bp,
                         lev_res=6,sigma_b = 0.2)
# or:
genotype_df <- as.data.frame(genotype)
res <- Wavelet_screaming( Y,loci=genotype_df,bp=my_bp,coeftype="c",
                         lev_res=6,sigma_b = 0.2)
```



### Assessing the significance

The test statistics of the Wavelet Screaming method relly of the distribution of the computed Bayes factors use during the computation. The Bayes factors distribution depend on two parameters, however if the sample size (individuals) is big enough one of this parameter, can be set as zero safely. So the Bayes factors depend only on one parameter that can be computed as follow:

(it takes 2-3 minutes)
```{r , echo=FALSE}
lambda <- get_lambda1(Y,sigma_b = 0.2)
```


We then simulate the null distribution of the test statistics depending of the parameter lambda and the Wavelet Screaming settings:
```{r , echo=FALSE}
Sim_gam <- Simu_Lambda_null(nsimu=10000, lambda=lambda,lev_res = 6)
```


#### Pvalue estimation

We use the simulated observations of the null distribution to fit an Generalized Pareto Distribution that correctly approximate the high quantile of the null distribution of the test statistics.  We then use this fitted distribution to estimated the associated quantile of the observation and so the  p-value.
```{r , echo=FALSE}
#Should be preferred for smaller values of Lambda
library(fExtremes)
x <-  Sim_gam
#Selecting the threshold
thresh <- min(  max (c( 2.5*mean( Sim_Gamma) , quantile(Sim_Gamma,0.98)) ), quantile(Sim_Gamma,0.97))
z = gpdFit(x, u = thresh , type = "mle")
z
pval <- 1-fExtremes::pgpd(q=res["Lambda"], xi=z@fit$par.ests["xi"],
                          mu=z@parameter$u, beta=z@fit$par.ests["beta"])
pval
```

### Visualisation

You can visualize your result with the function *plot_WS*
It represents the Bayes factor for the different levels scales of the wavelets decomposition.
The size and the darkness of the points that represent the Bayes factor are scaled by the value of the Bayes factors. If a Bayes factor is greater than 1 then the region that represent the Bayes factor is filled up in order to give an orverview of the size and the origin of the genetic signal.


```{r , echo=FALSE}
bp <- c(min(my_bp),max(my_bp))
plot_WS(res=res,bp=bp,lev_res=6)
```




# Full Genome screening

You can find in the folder *examples* our script for full genome screening.

## Slicing

 As the Wavelet Screaming is made to analysis region of the genome, the first step is to divide your data set in *screenable regions*. The function *slice_definition* is made for that, the user has to specify the size of the slice that he/she wants to analyse (Loci_size) AND the maximum gap between two SNPs (thresh option) (our set up Loci_size=1e6,thresh=1e4).

Below a script that slice vcf file and store the slice in a directory. This slice will be analyze in a second time

```{r , echo=FALSE}
library(WaveletScreaming)

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




## Phenotype management

IMPORTANT: The Wavelet Screaminng function assumes that your covariate, your counfunder matrix and your genotype matrix are sorted in the same order.



## Genome Wide Run
Below a code to deal with genome wide data.
In the *example* folder you will find detailed code to use the Wavelet Screaming package through the terminal in combinaison with *bcftools*.
We choose to store the slice instead of working directly in memory, it might not be optimal however it is more handy when you want to perform various model: different type of wavelet coefficients, different counfounding structure and different phenotype.
The provided scripts are already parallelised.


```{r , echo=FALSE}
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

# read phenotype and covariate data
pheno <- read.table(paste0(workdir, "bw/bw_moms.txt"), h=T)
Yall <- pheno$VEKT #Our phenoytpe

Confounderall <- as.matrix(pheno[,3:ncol(pheno)])#Where we put our Counfunding matrix 


# identify slices to read
slices = read.table(paste0(workdir, "bw/slicedef_moms_autosomes.txt"), h=F, sep=" ")
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


res_coefc = res_coefd  = data.frame()

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


  ## Coeftype D
  # Run the WaveletScreaming function
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screaming(Y=Yall,loci=my_slice,bp=bp,confounder=Confounderall,
  						  lev_res=9, sigma_b=0.2, coeftype="d", para=FALSE)

  # Save the results into a result matrix after each run.
  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefc <- rbind(res_coefc, t(res))
  print("analysis with coeftype C, no AA46, completed successfully")
  print(proc.time()-ptm)

  ##  Coeftype C
  ptm <- proc.time()
  temp <- NULL
  temp <- Wavelet_screaming(Y=Yall,loci=my_slice,bp=bp,confounder=Confounderall,
  						  lev_res=9, sigma_b=0.2, coeftype="c", para=FALSE)
  # Save the results into a result matrix after each run
  res <- c("chr"=slices[i,1],"start"=slices[i,2],"end"=slices[i,3],temp)
  res_coefd <- rbind(res_coefd, t(res))
  print("analysis with coeftype D, no AA46, completed successfully")
  print(proc.time()-ptm)


}

write.table(res_coefc,
			paste0(workdir, "results/res_coefc_", start, "-", end, ".txt"),
			row.names=F, col.names=T, quote=F)
write.table(res_coefd,
			paste0(workdir, "results/res_coefd_", start, "-", end, ".txt"),
			row.names=F, col.names=T, quote=F)

```
