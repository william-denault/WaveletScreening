#'@title Main function to perform wavelet screaming
#'@description  Perform a wavelet screening of a loci for a given phenotype and a specified level of resolution
#'@param Y phenotype vector, has to be numeric. For case control code it as 0 and 1. Multiple label phenotype will be implemented in the next version
#'@param loci genotype matrix (either data.frame or numeric matrix). Lines=SNPs in increasing order in term of base pair, columns=individuals. No missing values allowed.
#'@param bp vector of the base pairs positions. It has to be in the same order and length than the loci line order/length.
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param lev_res the maximum level of resolution needed
#'@param sigma_b the parameter of the NIG prior used for the betas computation. We advised to set this value between 0.1 and 0.2
#'@param coeftype type of wavelet coefficient used for the screening. By default set as "d" difference, "c" can be used if you prefere to work in term of amount of variants instead in disrepency within sub loci.
#'@param para logical parameter for parallelisation, if not specified set at FALSE.
#'@param BF logical parameter for getting Bayes Factor of wavelet regression, if not specified set at FALSE.
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
#'res
#'#value of the test statistics
#'res[c("L_h","min_ph_pv")]
#'#############
#'#Significance
#'#############
#'
#'#Simulate the null distribution using proxy covariance matrix
#'Sim <- Simu_null(Y,lev_res = 6,sigma_b = 0.2,size=10000)
#'head(Sim)
#'#Calibration of the hyperparameter
#'lambda <- Search_lambda(Sim,plot=TRUE)
#'
#'Th <- Sim[,c("L_h")]+lambda*Sim[,c("min_ph_pv")]
#'mu <- median(Th,na.rm = TRUE)
#'sdv <- mad(Th,na.rm = TRUE)
#'####################################
#'#Test Value of the loci to be tested
#'####################################
#'th <-  res[c("L_h")]+lambda*res["min_ph_pv"]
#'#######
#'#Pvalue
#'#######
#'1-pnorm(th,mean=muv,sd=sdv)
#'
#'hist(Th,nclass=1000,xlim=c(min(c(Th,th)),max(Th,th)))
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


Wavelet_screaming <- function(Y,loci,bp,confounder,lev_res,sigma_b,coeftype="d",para=FALSE,BF=FALSE)
{

  #Quantile transform to prevent for non normaliy distrib WCs
  Quantile_transform  <- function(x)
  {
    .ex.seed <- exists(".Random.seed")
    if(.ex.seed) .oldseed <- .Random.seed
    set.seed(666)
    if(.ex.seed) on.exit(.Random.seed <<- .oldseed)


    x.rank = rank(x, ties.method="random")
    #x.rank = rank(x, ties.method="average")
    return(qqnorm(x.rank,plot.it = F)$x)
  }

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
  	  Y <-   Quantile_transform(Y)
  	}
  }
  if(missing(BF)) {
    BF <- FALSE
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
    names_Betas <- c("Beta_0_0")
    for(i in 1:lev_res){
      for (j in 1:(2^i)){
        names_BF <- c(names_BF,paste("Beta",i,j,sep = "_"))
      }
    }
    out = rep(NA, 1+1+length(names_BF))
    names(out) <- c("L_h","min_ph_pv", names_BF)
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
    .ex.seed <- exists(".Random.seed")
    if(.ex.seed) .oldseed <- .Random.seed
    set.seed(666)
    if(.ex.seed) on.exit(.Random.seed <<- .oldseed)


    x.rank = rank(x, ties.method="random")
    #x.rank = rank(x, ties.method="average")
    return(qqnorm(x.rank,plot.it = F)$x)
  }



  max_EM_post_Beta <- function(my_betas, lev_res,null_sd,alt_sd,alp) {
    niter = 100
    epsilon <- 10^-4
    p_vec <- c()
    sd_vec <- c()
    erreur<-1+epsilon
    eps <-10^-10#slight correction in case of non identifiable mixture
    #prevent from having division by 0 in update parameter sum(temp)
    betasub = my_betas
    m0.hat<-0
    m1.hat<-0
    sigma0.hat<-sqrt(null_sd)
    sigma1.hat<-alt_sd
    #Prevent from label swapping
    if(sigma1.hat < sigma0.hat){
      sigma1.hat <- 3*sigma0.hat+sigma1.hat
    }

    p.hat<-0.5
    new.params<-c(m0.hat,m1.hat,sigma0.hat,sigma1.hat,p.hat)
    erreur<-1+epsilon
    iter <- 1

    while((erreur>epsilon)&(iter<=niter))
    {
      old.log.lik<- sum(log(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat)))
      old.params<-new.params
      #vecteur des pi_{i1}^{t}
      temp<-p.hat*dnorm( betasub ,m1.hat,sigma1.hat)/(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat))
      #Update parameter
      p.hat<-mean(temp)
      m1.hat<-sum(temp* betasub)/(sum(temp)+eps)
      m0.hat<-sum((1-temp)* betasub)/(sum(1-temp)+eps)
      sigma1.hat<-sqrt( sum(temp*( betasub-m1.hat)^2)/(sum(temp)+eps) )+alt_sd
      sigma0.hat<-sqrt( sum((1-temp)*( betasub-m0.hat)^2)/(sum(1-temp)+eps) )
      #limit the decrease of sigma0.hat in case of non identifiable mixture
      if(sigma0.hat < 0.5*sqrt(null_sd) ){
        sigma0.hat <- 0.5*sqrt(null_sd)
      }
      new.params<-c(m0.hat,m1.hat,sigma0.hat,sigma1.hat,p.hat)
      #Check end
      new.log.lik<- sum(log(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat)))
      #epsilon <- abs( new.log.lik -old.log.lik)
      iter<-iter+1

    }

    #Proba Belong belong to the alternative:
    pos.prob <- rep(NA,length(my_betas))
    for (i in 1:length(my_betas))
    {
      pos.prob[i] <- p.hat*dnorm( my_betas[i] ,m1.hat,sigma1.hat)/(p.hat*dnorm( my_betas[i] ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( my_betas[i] ,m0.hat,sigma0.hat))

    }

    lambcom <- rep(NA,(lev_res+1))
    p_vec   <- rep(NA,(lev_res+1))
    for (gi in 0:lev_res) {

      ####################################
      #Soft Thresholding
      ####################################
      temp <- pos.prob[(2^gi):(2^(gi + 1) - 1)]- alp*sqrt(1/2^(gi -1))
      pos.prob[(2^gi):(2^(gi + 1) - 1)] <-ifelse(temp<0,0, temp)
      ####################################
      #proportion of association per level
      ####################################
      p_vec[(gi+1)]    <-  mean(pos.prob[(2^gi):(2^(gi + 1) - 1)])
      lambcom[(gi+1)]  <-  mean(pos.prob[(2^gi):(2^(gi + 1) - 1)]*dnorm( betasub[(2^gi):(2^(gi + 1) - 1)] ,m1.hat,sigma1.hat)-(1-pos.prob[(2^gi):(2^(gi + 1) - 1)])*dnorm( betasub[(2^gi):(2^(gi + 1) - 1)] ,m0.hat,sigma0.hat))
    }

    porth   <- rep(NA,(lev_res))
    start <- 2^(1:lev_res)
    end  <-2^((1+1):(lev_res+1))-1
    for( gi in 0:(lev_res-1))
    {
      tempstart <- round(start + (gi)*(end-start)/lev_res)
      tempend <-round(start + (gi+1)*(end-start)/lev_res)
      ind <- c()
      for ( j in 1:lev_res)
      {
        temp <- cbind(tempstart,tempend)
        p1 <- temp[j,1]
        p2 <- temp[j,2]
        ind <- c(p1:p2 )
      }
      porth[gi+1] <-  mean(pos.prob[ind])

    }
    ph <- sum(p_vec)
    pv <- sum(porth)
    min_ph_pv <- min( ph,pv)
    L_h <- sum(lambcom)
    return(c(L_h, min_ph_pv))
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


  print("Computing Beta values")
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
    my_bayes <- apply(Gen_W_trans, 2, my_bf )
  }

  ###########################
  #Computation test statstics
  ###########################
  print("Post-processing")

  Dmat <- cbind(confounder,Y)
  Dmat <- as.matrix(Dmat)

  index <- dim(confounder)[2]
  resM <- (1/sigma_b/sigma_b)*solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,dim(Dmat)[2]))
  #Starting position for the EM
  null_sd <-sigma_b*as.numeric(resM["Y","Y"])^2
  alt_sd <- sigma_b
  #Shrinkage coefficient for the EM
  alp <-  1/sqrt(2*log(length(Y)))
  my_betas <- as.numeric(my_betas)
  test_stat <- max_EM_post_Beta(my_betas=my_betas, lev_res = lev_res, null_sd =  null_sd, alt_sd = alt_sd,alp = alp)





  if(BF ==FALSE)
  {
    out <- c(test_stat,my_betas)
  }
  else
  {
    out <-  c(test_stat,my_betas,my_bayes)
  }

  #Naming the output
  if(BF ==FALSE)
  {
   names_Betas <- c("Beta_0_0")
   for(i in 1:lev_res)
   {
     for (j in 1:(2^i))
     {
       names_Betas <- c(names_Betas,paste("BF",i,j,sep = "_"))
     }
   }

   names(out) <- c("L_h",
                   "min_ph_pv",
                   names_Betas)
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
    for(i in 1:lev_res)
    {
      for (j in 1:(2^i))
      {
        names_Betas <- c(names_Betas,paste("Beta",i,j,sep = "_"))
      }
    }
    names(out) <- c("L_h",
                    "min_ph_pv",
                    names_Betas,names_BF)
  }




  if(para==TRUE)
  {
    stopCluster(cl)
  }
  return(out)
}

