#'@title Simulation of the null statsitics
#'@description  Simulation of the null statsitics
#'@param Y a vector of numeric values used in in the wavelet screaming function for association (phenotype or simulated phenotype with same distribution).
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param lev_res the level of resolution in the wavelet transform
#'@param size number of simulation to be performed
#'@param sigma_b the parameter of the NIG prior used for the Betas computation.
#'@param print logical parameter set as TRUE, if TRUE send message when 10\% of the simulations has been completed.
#'@return The simulation under the null of the two test statistics used to build the final test (i.e L_h and min(ph,pv))
#'@examples \dontrun{
#'Y <- rnorm(4000)
#'Sim <- Simu_null_proxy(Y=Y,lev_res = 6,sigma_b = 0.2,size=10000)
#'}
Simu_null_proxy <- function(Y,confounder,lev_res,size,sigma_b,print=TRUE)
{
  smp_size= length(Y)
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
  #####################################
  #Size of the multivaraite to simulate
  #####################################
  nbeta <- sum(2^(0:lev_res))
  ######################################
  #compute theoretical variance of Betas
  ######################################
  if(missing(confounder))
  {
    print("no covariates provided, using intercept only")
    Dmat <- cbind(rep(1,length(Y)),Y)
    Dmat <- as.matrix(Dmat)
  }
  else if(nrow(confounder)!=length(Y)) {
    stop("ERROR: number of samples in Y and confounder does not match")
  } else {
    print(sprintf("%i covariates for %i samples detected", ncol(confounder), nrow(confounder)))
    confounder <- cbind(rep(1,length(Y)),confounder)
    Dmat <- cbind(confounder,Y)
    Dmat <- as.matrix(Dmat)
  }


  sigma_b <- sigma_b
  beta_0 <- c()
  for( i in 1:1000)
  {
    y <- rnorm(length(Y),sd=1)
    temp <- solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,dim(Dmat)[2])) %*% t(Dmat)%*% y
    beta_0 <- c(beta_0,temp[2])

  }


  null_sd <-var(beta_0)

  ##################################
  #Computing proxy covariance matrix
  ##################################
  print("Empirical Covariance missing computing proxy covariance for simulation")
  lev_res <-lev_res
  coeftype="c"
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

  N=length(Y)
  #Generate random signal
  SNP=lev_res*500
  loci<- matrix(runif(N*SNP,min=0,2),ncol=N)
  bp= 1:SNP
  #Preparing for the wavelet transform
  Time01 <- (bp- min(bp))/(max(bp)-min(bp))
  Gen_W_trans <- apply(loci,2,my_wavproc)
  Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)
  #Compute proxy for empirical covariance matrix
  emp_cov <- (cov(Gen_W_trans))*null_sd
  print("Proxy Covariance computed")




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
      if(sigma0.hat < 0.1*sqrt(null_sd) ){
        sigma0.hat <- 0.1*sqrt(null_sd)
      }
      new.params<-c(m0.hat,m1.hat,sigma0.hat,sigma1.hat,p.hat)
      #Check end
      new.log.lik<- sum(log(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat)))
      epsilon <- abs( new.log.lik -old.log.lik)
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
  ######################
  #Set up for simulation
  ######################
  null_sd <- mean(diag(emp_cov))
  Pi_nt <- list()
  alt_sd <- sigma_b
  alp <- 1/sqrt(2*log(smp_size))
  print("Simulation of test statistics")
  temp <- seq(from=1,to=size,by=size/10)[-1]-1
  out <- list()
  my_f <- function(y)
  {
    pis <- max_EM_post_Beta(y, lev_res = 9, null_sd = null_sd ,alt_sd = alt_sd,alp=alp)
    return(pis)
  }
  if(print==TRUE)
  {

    for (i in 1:10)
    {
      y <- rmvnorm(floor(size/10),mean=rep(0,nbeta),sigma =emp_cov )
      out[[i]] <-t( apply(y,1,my_f))
      print(paste(i*floor(size/10), "simulations performed"))
    }

  }
  if(print==FALSE)
  {
    for (i in 1:10)
    {
      y <- rmvnorm(floor(size/10),mean=rep(0,nbeta),sigma =emp_cov )
      out[[i]] <-t( apply(y,1,my_f))
    }
  }



  out <- do.call(rbind,out)
  colnames(out) <- c("L_h","min_ph_pv")
  return(out)
}
