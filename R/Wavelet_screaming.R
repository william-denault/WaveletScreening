#'@title Main function to perform wavelet screaming
#'@description  Perform a wavelet screening of a locus for a given phenotype and a specified level of resolution
#'@param Y phenotype vector has to be numeric. For case-control  data, code it as 0 and 1. Multiple label phenotypes, e.g., ABO blood groups, will be implemented in the next version.
#'@param loci genotype matrix (either data.frame or numeric matrix). Lines=SNPs in increasing order in terms of base pair, columns=individuals. No missing values allowed.
#'@param bp vector of the base pairs positions. It has to be in the same order and length than the locus line order/length.
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included if missing will generate an intercept matrix.
#'@param lev_res the maximum level of resolution needed.
#'@param sigma_b the parameter of the NIG prior used for the Beta computation. We advised setting this value between 0.1 and 0.2
#'@param coeftype type of wavelet coefficient used for the screening (choice "c" or "d"). If missing set as "d"
#'@param base_shrink numeric, value used in the thresholding of the proportion of assocation, if non specificed set up as 1/sqrt(2*log(sample_size))
#'@param para logical parameter for parallelization, if not specified, set at FALSE by default.
#'@param BF logical parameter for obtainning the Bayes Factor of the wavelet regression. If not specified, set at FALSEby default.
#'@details The Wavelet_screaming function computes the likelihood ratio used for testing the significance of a genetic region. In addition, it computes
#'the proportion of wavelet coefficients associated by the level of resolution and the Beta used for this estimation. All the details
#'of the computation can be found in our paper, preliminarily titled "Wavelet Screaming: a novel look to GWAS data.".
#'@return A named vector. The first position contains the estimated value of the Lambda statistics. The next positions of the vector are the computed proportion of associations per level of resolution.
#'set.seed(1)
#'#########################################
#'#Generate a randomly sampled SNP from a locus of size=1Mb
#'#########################################
#'
#'#5000 Randomly choosen basepairs
#'my_bp <- sort(sample(1:1000000, size=5000,replace = FALSE))
#'#############################
#'#Three different bump signals
#'#############################
#'my_functions <-data.frame(f0 = c(rep(0,400000),rep(0,200000),rep(0,400000)),
#'                          f1 = c(rep(0,400000),rep(1,200000),rep(0,400000)) ,
#'                          f2=c(rep(0,400000),rep(2,200000),rep(0,400000)))
#'
#'
#'library(gridExtra)
#'###########################
#'#Minor allele frequency 30%
#'###########################
#'MAF=0.3
#'sampl_schem <- c((1-MAF)^2,2*MAF*(1-MAF),MAF^2)
#'#######################################
#'#Sampling at Hardy Weinberg equilibrium
#'#######################################
#'#Assigning class
#'
#'#sample size =4000
#'n_size=4000
#'type_fn <-sample(0:2,replace = TRUE,size=n_size,  prob=  sampl_schem  )
#'
#'
#'genotype <-  matrix(my_functions[my_bp,2 ], ncol=1 ) %*%t(matrix(type_fn,ncol=1))
#'genotype <- genotype+ rnorm(dim(genotype)[1]*dim(genotype)[2],sd=0.1)
#'#dim(genotype)= nSNP, nind
#'
#'###############################################################
#'#Generate a phenotype with variance explained by genotype =0.5%
#'###############################################################
#'varexp=0.005
#'var_noise <- (1-varexp)*var(sample(0:2,replace = TRUE,size=10000,
#'                                   prob=sampl_schem ))/varexp
#'Y <-  rnorm(n=n_size,sd=sqrt(var_noise)) +type_fn
#'df <- data.frame(y=Y,genotype =factor(type_fn))
#'P1 <- ggplot(df,aes(y=y,x=genotype))+
#'  geom_boxplot()+
#'  xlab("Type of genotype")+
#'  theme(axis.text=element_text(size=12),
#'        axis.title=element_text(size=14,face="bold"))+
#'  ylab("Simulated Phenotype")+
#'  theme_bw()+
#'  ggtitle("Variation of the phenotype\ndepending of the genotype, \nVariance explained =0.5%")
#'
#'df <- data.frame(bp= rep(my_bp,3),y=c(my_functions[my_bp,1],my_functions[my_bp,2],my_functions[my_bp,3]),
#'                 mycol = factor(c(rep("f0",length(my_bp)),rep("f1",length(my_bp)),rep("f2",length(my_bp))) ) )
#'
#'P2 <- ggplot(df,aes(y=y,x=bp,color=mycol))+
#'  geom_point(size=1)+
#'  xlab("Base pair")+
#'  ylab("Number of variants")+
#'  theme_bw()+
#'  theme(legend.title=element_blank())+
#'  ggtitle("Three different kind of genotype signal")
#'
#'grid.arrange(P1,P2,ncol=2)
#'
#'##################
#'#Wavelet screaming
#'##################
#'res <- Wavelet_screaming( Y,loci=genotype,bp=my_bp,
#'                          lev_res=6,sigma_b = 0.2)
#'res
#'#Value of the test statistics
#'res[c("L_h","min_ph_pv")]
#'#############
#'#Significance
#'#############
#'
#'#Simulate the null distribution using proxy covariance matrix
#'
#'Sim <- Simu_null_proxy(Y,lev_res = 6,sigma_b = 0.2,size=10000)
#'head(Sim)
#'#Calibration of the hyperparameter
#'lambda <- Search_lambda(Sim,plot=TRUE)
#'
#'Th <- Sim[,c("L_h")]+lambda*Sim[,c("min_ph_pv")]
#'muv <- median(Th,na.rm = TRUE)
#'sdv <- mad(Th,na.rm = TRUE)
#'####################################
#'#Test Value of the loci to be tested
#'####################################
#'th <-  res[c("L_h")]+lambda*res["min_ph_pv"]
#'#######
#'#P-value
#'#######
#'1-pnorm(th,mean=muv,sd=sdv)
#'
#'df <- data.frame(Th = Th,type = factor( c(rep("Null",length(Th)))) )
#'ggplot(df,aes(Th,fill=type))+
#'  xlim(c(min(c(Th,th)),max(Th,th)))+
#'  geom_density()+
#'  guides(fill=FALSE)+
#'  geom_point(aes(x=th, y=0), colour="red")+theme(legend.position="none")+
#'  geom_text(label="Value of the test statistics", x=th, y=0.001)+
#'  geom_text(label="Null distribution", x=mean(Th), y=0.001)+
#'  theme_bw()
#'##############
#'#Visualisation
#'##############
#'bp <- c(min(my_bp),max(my_bp))
#'plot_WS(res=res,bp=bp,lev_res=6)
#'
#'
#'}


Wavelet_screaming <- function(Y,
                              loci,
                              bp,
                              confounder,
                              lev_res,
                              sigma_b,
                              coeftype,
                              base_shrink,
                              para=FALSE,
                              BF=FALSE)
{



  #To ensure the length not to be 0
  Y <- as.vector(Y)
  sigma_b <- sigma_b


  # INPUT CHECKS
  if(missing(coeftype))
  {
    print( "missing coeftype set as d")
    coeftype <- "d"
  }
   message("Input dimensions:")
  if(!is.numeric(Y) || length(Y)==0){
    stop("ERROR: Y is not a numeric vector")
  } else {
   	message(sprintf("%i phenotypes detected", length(Y)))
  	if(all(Y %in% c(0,1))){
   		message("Binary phenotype detected")
  	} else if(!is.vector(Y)){
  		stop("ERROR: Y is not a vector. Multi-phenotype analysis not implemented yet.")
  	} else {
   		message("Continuous phenotype detected")
  	}
  }
  # Writing the design matrix
  if(missing(confounder)) {
   	message("no covariates provided, using intercept only")
  	confounder <- data.frame(confounding=rep(1,length(Y)) )
  } else if(nrow(confounder)!=length(Y)) {
    stop("ERROR: number of samples in Y and confounder does not match")
  } else {
   	message(sprintf("%i covariates for %i samples detected", ncol(confounder), nrow(confounder)))
	  confounder <- cbind(rep(1,length(Y)),confounder)
  }
  if(missing(BF)) {
   BF <- FALSE
  }


  # Check genotype matrix
  if(is.data.frame(loci)){
   	message("Converting genotype data to matrix")
  	loci <- as.matrix(loci)
  }
  if(missing(loci) || !is.numeric(loci)){
    stop("ERROR: genotype matrix missing or not numeric")
  } else if(ncol(loci)!=length(Y)){
    stop("ERROR: number of samples in Y and loci does not match")
  } else {
   	message(sprintf("%i SNPs for %i samples detected", nrow(loci), ncol(loci)))
  }

  # Check position vector
  if(!is.numeric(bp) || !is.vector(bp)){
    stop("ERROR: must provide numeric position vector")
  } else {
   	message(sprintf("positions for %i SNPs read", length(bp)))
  }

  # Clean missing samples from all inputs
  keepY <- complete.cases(Y)
  keepC <- complete.cases(confounder)
  keepGT <- complete.cases(t(loci))
  nonmissing_index <- which(keepGT & keepY & keepC)
  if(length(nonmissing_index) != length(Y)){
	warning(sprintf("Warning: %i individuals will be removed due to missingness",
  			length(Y) - length(nonmissing_index)))
  }

  Y <- Y[nonmissing_index]
  confounder <- confounder[nonmissing_index,]
  loci <- loci[,nonmissing_index]

   message(paste("N individuals analysed = ", dim(loci)[2],
   			", N SNPs analysed = ",dim(loci)[1]))

  # workaround for git issue #1 - mysteriously empty slices
  if(is.null(dim(loci)) || dim(loci)[1] < 2^lev_res || dim(loci)[2] < 2){
	warning("not enough genotypes remaining, returning empty output")

    # Naming the output
    names_Betas <- c("Beta_0_0")
    for(i in 1:lev_res){
      for (j in 1:(2^i)){
        names_Betas <- c(names_Betas,paste("Beta",i,j,sep = "_"))
      }
    }
    out = rep(NA, 1+1+length(names_Betas))
    names(out) <- c("L_h","min_ph_pv",names_Betas)
    return(out)
  }


  ####################################
  #Redefinition of the needed function
  ####################################
  n_coef_wc <- sum(2^(0:4))


  ###############
  #Paralelisation
  ###############
  if(para==TRUE)
  {
    cl <-makeCluster(detectCores(all.tests=TRUE)-1, type = "SOCK")
  }

  ###################
  #Wavelet processing
  ###################
   message("Wavelet processing")

  Time01 <- (bp- min(bp))/(max(bp)-min(bp))
  my_wavproc <- function(y)
  {
    #Kovac and Silvermann 2000
    mygrid <- wavethresh::makegrid(t=Time01,y=y)
    LDIRWD <- irregwd(mygrid,filter.number=1)
    class(LDIRWD) <- "wd"
    #Thresholding here
    LDIRWD <- threshold(LDIRWD,policy = "universal",type="hard",
                        dev = madmad,levels = 1:(LDIRWD$nlevels-1))

    res <- c()
    for(i in 0: lev_res){
      if(coeftype == "d"){
        res <- c(res, accessD( LDIRWD,lev = i) )
      } else if (coeftype == "c") {
        res <- c(res, accessC( LDIRWD,lev = i) )
      } else {
        stop(paste("ERROR: coeftype", coeftype, "not recognized"))
      }
    }

    return(res)
  }

  if(para==TRUE)
  {
    clusterExport(cl,"irregwd")
    clusterExport(cl,"threshold")
    clusterExport(cl,"madmad")
    clusterExport(cl,"accessD")
    clusterExport(cl,"accessC")
    clusterExport(cl,"my_wavproc")
    Gen_W_trans <- snow::parApply(cl,loci,2,my_wavproc)
  }
  else{
    Gen_W_trans <- apply(loci,2,my_wavproc)
  }

  #Quantile transform for non normal WCs for every scale location
  Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)

  ##########
  #Modeling
  ##########


  message("Computing Beta values")
  betas_f <- function(y)
  {

    confounder <- data.frame(confounder)
    pc <- dim(confounder)[2]
    Dmat <- cbind(confounder,Y)
    Dmat <- as.matrix(Dmat)

    res <- solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,dim(Dmat)[2])) %*% t(Dmat)%*% y
    index <- pc+1

    return(res[index,1])
  }


  if(para==TRUE)
  {

    clusterExport(cl,"betas_f")
    my_betas <- snow::parApply(Gen_W_trans, 2, betas_f )
  }
  else{
    my_betas <- apply(Gen_W_trans, 2, betas_f )
  }

  if(BF ==TRUE)
  {
    message("Computing Bayes Factors")
    W <- as.matrix(confounder, ncol=ncol(confounder))
    n = nrow(W)
    q = ncol(W)

    # L <- as.matrix(Y , ncol=ncol(Y)) #reversed regression
    L <- as.matrix(Y,ncol=1)

    p = 1
    PW = diag(n) - W %*% solve(t(W) %*% W) %*% t(W)
    X = PW %*% L
    HB = X %*% solve(t(X) %*% X + diag(1/sigma_b/sigma_b,p)) %*% t(X)
    delta = svd(X)$d
    lambda = delta^2 / (delta^2 + 1/sigma_b/sigma_b)
    log.T = sum(log(1-lambda))/2

    my_bf <- function( y ){
      y <-  as.matrix(y,ncol=1)
      log.R = -0.5*n*log(1 - (t(y) %*% HB %*% y) / (t(y) %*% PW %*% y ))

      bf = exp(log.T + log.R)
      return(c(bf))
    }
    my_bayes <- apply(Gen_W_trans, 2, my_bf )
  }

  ###########################
  #Computation test statstics
  ###########################
  message("Post-processing")

  Dmat <- cbind(confounder,Y)
  Dmat <- as.matrix(Dmat)

  null_sd <- sqrt(solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,dim(Dmat)[2]))["Y","Y"])
  alt_sd <- 100*null_sd
  #Shrinkage coefficient for the EM

  if(missing(base_shrink))
  {
    alp <-  1/sqrt(2*log(length(Y)))
  }
  else
  {
    alp <-  base_shrink
  }

  my_betas <- as.numeric(my_betas)
  rest <- max_EM_post_Beta(my_betas=my_betas, lev_res = lev_res, null_sd =  null_sd, alt_sd = alt_sd,alp = alp)

  test_stat <- rest[[1]]
  postH1 <- rest[[2]]


  if(BF ==FALSE)
  {
    out <- c(test_stat,my_betas,postH1)
  }
  else
  {
    out <-  c(test_stat,my_betas,postH1,my_bayes)
  }

  #Naming the output
  if(BF ==FALSE)
  {
    names_Betas <- c("Beta_0_0")
    names_postH1 <-  c("Pi_0_0")
    for(i in 1:lev_res)
    {
      for (j in 1:(2^i))
      {
        names_Betas <- c(names_Betas,paste("Beta",i,j,sep = "_"))
        names_postH1 <- c(names_postH1,paste("Pi",i,j,sep = "_"))
      }
    }

    names(out) <- c("L_h",
                    "min_ph_pv",
                    names_Betas,
                    names_postH1)
  }
  else
  {
    names_BF <- c("BF_0_0")
    for(i in 1:lev_res)
    {
      for (j in 1:(2^i))
      {
        names_BF <- c(names_BF,paste("BF",i,j,sep = "_"))
      }
    }
    names_Betas <- c("Beta_0_0")
    names_postH1 <-  c("Pi_0_0")
    for(i in 1:lev_res)
    {
      for (j in 1:(2^i))
      {
        names_Betas <- c(names_Betas,paste("Beta",i,j,sep = "_"))
        names_postH1 <- c(names_postH1,paste("Pi",i,j,sep = "_"))
      }
    }
    names(out) <- c("L_h",
                    "min_ph_pv",
                    names_Betas,
                    names_postH1,
                    names_BF)
  }




  if(para==TRUE)
  {
    stopCluster(cl)
  }
  return(out)
}
