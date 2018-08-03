#'@title Main function to perform wavelet screaming
#'@description  Perform a wavelet screening of a loci for a given phenotype and a specified level of resolution
#'@param Y phenotype vector, has to be numeric. For case control code it as 0 and 1. Multiple label phenotype will be implemented in the next version
#'@param loci genotype matrix (either data.frame or numeric matrix). Lines=SNPs in increasing order in term of base pair, columns=individuals. No missing values allowed.
#'@param bp vector of the base pairs positions. It has to be in the same order and length than the loci line order/length.
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param lev_res the maximum level of resolution needed
#'@param sigma_b the parameter of the NIG prior used for the Bayes Factor computation. We advised to set this value between 0.1 and 0.2
#'@param coeftype type of wavelet coefficient used for the screening. By default set as "d" difference, "c" can be used if you prefere to work in term of amount of variants instead in disrepency within sub loci.
#'@param para logical parameter for parallelisation, if not specified set at FALSE.
#'@details The Wavelet_screaming function computes the Likelihood ratio used for testing significance of a genetic region. In addition it computes
#'the porportion of wavelets coefficients associated by level of resolution, and the Bayes factor used for this estimation. All the details
#'of the computation can be found in our paper "Wavelet Screaming: a novel look to GWAS data"
#'@return A named vector. First position the estimated value of the Lambda statistics, then the proportion of association per level of resolution then the computed Bayes Factor per wavelet coefficient.
#'@examples \dontrun{
#'
#'
#'set.seed(666)
#'#########################################
#'#Generate a randomly sample loci size=1Mb
#'#########################################
#'
#'#5000 Randomly choosen bp
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
#'#sampling at Hardy Weinberg equilibrium
#'#######################################
#'#Assigning class
#'
#'#sample size =4000
#'n_size=4000
#'type_fn <-sample(0:2,replace = TRUE,size=n_size,  prob=  sampl_schem  )
#'
#'
#'genotype <-  matrix(my_functions[my_bp,2 ], ncol=1 ) %*%t(matrix(type_fn,ncol=1))
#'#dim(genotype)= nSNP, nind
#'
#'###############################################################
#'#Generate a phenotype with variance explained by genotype  0.5%
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
#'# or:
#'genotype_df <- as.data.frame(genotype)
#'res <- Wavelet_screaming( Y,loci=genotype_df,bp=my_bp,
#'                          lev_res=6,sigma_b = 0.2)
#'#############
#'#Significance
#'#############
#'
#'#Estimation of the Bayes factor lambda_1 distribution parameter
#'#Take a bit of time
#'lambda <- get_lambda1(Y,sigma_b = 0.2)
#'#Simulation of the test statistics under the null distribution of
#'# the Bayes factor
#'Sim_gam <- Simu_Lambda_null(nsimu=10000, lambda=lambda,lev_res = 6)
#'val <- res["Lambda"]
#'
#'#Via Simulation
#'
#'pval <-c(length(Sim_gam[which(Sim_gam>val)])+1)/(length(Sim_gam)+1)
#'pval
#'
#'#Via EVT
#'#Should be preferred for smaller values of Lambda
#'
#'library(fExtremes)
#'x <-  Sim_gam
#'z = gpdFit(x, u = min(x), type = "mle")
#'z
#'pval <- 1-fExtremes::pgpd(q=res["Lambda"], xi=z@fit$par.ests["xi"], mu=z@parameter$u, beta=z@fit$par.ests["beta"])
#'pval
#'
#'
#'##############
#'#Visualisation
#'##############
#'bp <- c(min(my_bp),max(my_bp))
#'plot_WS(res=res,bp=bp,lev_res=6)
#'
#'
#'}


Wavelet_screaming <- function(Y,loci,bp,confounder,lev_res,sigma_b,coeftype="d",para=FALSE)
{
  #Loci: genotype matrix, line=SNP order in increasing bp, column individual genoype
  #bp: position of the SNP in term of base pair
  #confounder: designed matrix of the confounding effect size = n,c
  #n= n ind, c= number of confounder
  #lev_res: lev of resolution for the wavelet filtering
  #sigma_b= Para of prior, should be <1 advised 0.2


  #To ensure the length not to be 0
  Y <- as.vector(Y)
  sigma_b <- sigma_b


  # INPUT CHECKS
  print("Input dimensions:")
  if(!is.numeric(Y) || length(Y)==0){
    stop("ERROR: Y is not a numeric vector")
  } else {
  	print(sprintf("%i phenotypes detected", length(Y)))
  	if(all(Y %in% c(0,1))){
  		print("Binary phenotype detected")
  	} else if(!is.vector(Y)){
  		stop("ERROR: Y is not a vector. Multi-phenotype analysis not implemented yet.")
  	} else {
  		print("Continuous phenotype detected")
  	}
  }

  # Writing the design matrix
  if(missing(confounder)) {
  	print("no covariates provided, using intercept only")
  	confounder <- data.frame(confounding=rep(1,length(Y)) )
  } else if(nrow(confounder)!=length(Y)) {
  	stop("ERROR: number of samples in Y and confounder does not match")
  } else {
  	print(sprintf("%i covariates for %i samples detected", ncol(confounder), nrow(confounder)))
	  confounder <- cbind(rep(1,length(Y)),confounder)
  }


  # Check genotype matrix
  if(is.data.frame(loci)){
  	print("Converting genotype data to matrix")
  	loci <- as.matrix(loci)
  }
  if(missing(loci) || !is.numeric(loci)){
  	stop("ERROR: genotype matrix missing or not numeric")
  } else if(ncol(loci)!=length(Y)){
  	stop("ERROR: number of samples in Y and loci does not match")
  } else {
  	print(sprintf("%i SNPs for %i samples detected", nrow(loci), ncol(loci)))
  }

  # Check position vector
  if(!is.numeric(bp) || !is.vector(bp)){
  	stop("ERROR: must provide numeric position vector")
  } else {
  	print(sprintf("positions for %i SNPs read", length(bp)))
  }

  # Clean missing samples from all inputs
  keepY <- complete.cases(Y)
  keepC <- complete.cases(confounder)
  keepGT <- complete.cases(t(loci))
  nonmissing_index <- which(keepGT & keepY & keepC)
  if(length(nonmissing_index) != length(Y)){
  	print(sprintf("Warning: %i individuals will be removed due to missingness",
  			length(Y) - length(nonmissing_index)))
  }

  Y <- Y[nonmissing_index]
  confounder <- confounder[nonmissing_index,]
  loci <- loci[,nonmissing_index]

  print(paste("N individuals analysed = ", dim(loci)[2],
  			", N SNPs analysed = ",dim(loci)[1]))

  # workaround for git issue #1 - mysteriously empty slices
  if(is.null(dim(loci)) || dim(loci)[1] < 2^lev_res || dim(loci)[2] < 2){
  	print("Warning: not enough genotypes remaining, returning empty output")

    # Naming the output
    names_BF <- c("BF_0_0")
    for(i in 1:lev_res){
  	  for (j in 1:(2^i)){
  		names_BF <- c(names_BF,paste("BF",i,j,sep = "_"))
  	  }
    }
    out = rep(NA, 1+lev_res+1+length(names_BF))
    names(out) <- c("Lambda", paste("pi",0:lev_res, sep = "_"), names_BF)
    return(out)
  }

  ####################################
  #Redefinition of the needed function
  ####################################
  n_coef_wc <- function(lev_res)
  {
    temp <- c()
    for(i in 0:lev_res)
    {
      temp <- c(temp,2^i)
    }
    sum(temp)
  }

  #Quantile transform to prevent for non normaliy distrib WCs
  Quantile_transform  <- function(x)
  {
    x.rank = rank(x, ties.method="random")
    return(qqnorm(x.rank,plot.it = F)$x)
  }

  #Estimation of Lambda
  Lambda_stat <- function (my_pi, my_bayes)
  {
    # vector: pi1 pi2 pi2 pi3 pi3 pi3 pi3...
    my_pi_vec = rep(my_pi, 2^(1:length(my_pi)-1))
    coefs = 1-my_pi_vec + my_pi_vec * my_bayes[1:(2^length(my_pi)-1)]
    prod(coefs)
  }

  sumlog <- function (A1, A2)
  {
    if(A1 > A2){
      res = A1 + log(1 + exp(A2 - A1))
    }else{
      res = A2 + log(exp(A1 - A2) + 1)
    }

    return(res)
  }

  max_EM_Lambda <- function(my_bayes)
  {
    niter=10000
    epsilon <- 10^-4
    p_vec <- c()
    for(gi in 0: lev_res)
    {
      # EM algorithm for each group separately

      N_obllikli = 0
      logpi = log(0.5)
      pi <- 0.5
      log1pi = logpi

      pp = 0
      logPiBF = log(my_bayes[(2^gi):(2^(gi+1)-1)]) + logpi
      logden <- c()
      for (i in 1:length(logPiBF))
      {
        logden[i] <- sumlog(logPiBF[i],log1pi)
      }
      pp = pp+sum(exp(logPiBF - logden))
      N_obllikli = sum(logden)
      O_obllikli = N_obllikli


      for(iter in  0:niter){
        pi = pp/(2^(gi))
        logpi  = log(pi)
        log1pi = log(1-pi)
        logPiBF =   log(my_bayes[(2^gi):(2^(gi+1)-1)]) + logpi
        logden <- c()
        for (i in 1:length(logPiBF))
        {
          logden[i] <- sumlog(logPiBF[i],log1pi)
        }
        pp=0
        pp = pp+sum(exp(logPiBF - logden))
        N_obllikli = sum(logden)
        diff = abs(N_obllikli - O_obllikli)

        if(diff < epsilon){
          break
        }else{
          O_obllikli = N_obllikli
        }
      }
      p_vec <-c(p_vec,pi)
    }
    return(p_vec)
  }


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
  print("Wavelet processing")

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
  print("Computing Bayes Factors")
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


  if(para==TRUE)
  {
    clusterExport(cl,"log.T")
    clusterExport(cl,"sigma_b")
    clusterExport(cl,"my_bf")
    my_bayes <- snow::parApply(cl,Gen_W_trans, 2, my_bf )
  }
  else{
    my_bayes <- apply(Gen_W_trans, 2, my_bf )
  }


  #################
  #Estimation Lambda
  #################
  print("Post-processing")
  my_pis <- max_EM_Lambda(my_bayes = my_bayes)
  trueLambda <- Lambda_stat(my_pi = my_pis,my_bayes = my_bayes)

  out <- c(trueLambda,my_pis,my_bayes)

  #Naming the output
  names_BF <- c("BF_0_0")
  for(i in 1:lev_res)
  {
    for (j in 1:(2^i))
    {
      names_BF <- c(names_BF,paste("BF",i,j,sep = "_"))
    }
  }

  names(out) <- c("Lambda",
                  paste("pi",0:lev_res, sep = "_"),
                  names_BF)
  if(para==TRUE)
  {
    stopCluster(cl)
  }
  return(out)
}

