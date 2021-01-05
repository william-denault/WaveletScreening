#'@title Simulation of the null statistics
#'@description  Simulation of the null statistics
#'@param Y a vector of numeric values used in the wavelet screening function for association (a phenotype or simulated phenotype with the same distribution).
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included if missing will generate an intercept matrix.
#'@param lev_res the level of resolution in the wavelet transform
#'@param size number of simulation to be performed. If not specified set at 10000
#'@param base_shrink numeric, value used in the thresholding of the proportion of assocation, if non specificed set up as 1/sqrt(2*log(sample_size)
#'@param sigma_b the parameter of the NIG prior used for the Betas computation.
#'@param coeftype type of wavelet coefficient used for the screening (choice "c" or "d"). If missing set as "c"
#'@param verbose logical parameter, set as TRUE by default. ID
#'@return The simulation under the null of the two test statistics used to build the final test (i.e., L_h and min(ph,pv))
#'@examples \dontrun{
#'Y <- rnorm(4000)
#'Sim <- Simu_null_proxy(Y=Y,lev_res = 6,sigma_b = 0.2,size=10000)
#'}
Simu_null_proxy <- function(Y,
                            confounder,
                            lev_res,
                            size=10000,
                            sigma_b= NA,
                            base_shrink,
                            coeftype= "c",
                            verbose= TRUE)
{
  smp_size= length(Y)
  if(missing(coeftype))
  {
    print( "missing coeftype set as c")
    coeftype <- "c"
  }

  if( is.na(sigma_b))
  {
    if(verbose)
    {
      message("No prior size provided, using frequentist modeling")
    }
    analysis_type <- "Frequentist"
  }
  if( !is.na(sigma_b))
  {
    if(verbose)
    {
      message("Using Bayesian modeling")
    }
    analysis_type <- "Bayesian"
  }
  Quantile_transform  <- function(x)
  {

    x.rank = rank(x, ties.method="random")
    #x.rank = rank(x, ties.method="average")
    return(qqnorm(x.rank,plot.it = F)$x)
  }
  # INPUT CHECKS
  if(verbose)
  {
   print("Input dimensions:")
  }
  if(!is.numeric(Y) || length(Y)==0){
    stop("ERROR: Y is not a numeric vector")
  } else {
    if(verbose)
    {
      print(sprintf("%i phenotypes detected", length(Y)))
    }

    if(all(Y %in% c(0,1))){
      if(verbose)
      {
      print("Binary phenotype detected")
      }
    } else if(!is.vector(Y)){
      stop("ERROR: Y is not a vector. Multi-phenotype analysis not implemented yet.")
    } else {
      if(verbose)
      {
        print("Continuous phenotype detected")
      }
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
    if(verbose)
    {
     print("no covariates provided, using intercept only")
    }
    Dmat <- cbind(rep(1,length(Y)),Y)
    Dmat <- as.matrix(Dmat)
  }
  else
  {
    if(nrow(confounder)!=length(Y)) {
      stop("ERROR: number of samples in Y and confounder does not match")
    } else {
      if(verbose)
      {
       print(sprintf("%i covariates for %i samples detected", ncol(confounder), nrow(confounder)))
      }
      confounder <- cbind(rep(1,length(Y)),confounder)
      Dmat <- cbind(confounder,Y)
      Dmat <- as.matrix(Dmat)
    }

  }



  sigma_b <- sigma_b

  Dmat <- as.matrix(Dmat)

  if( analysis_type == "Bayesian" )
  {
    null_sd <- sqrt(solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,dim(Dmat)[2]))["Y","Y"])
  }
  if( analysis_type == "Frequentist" )
  {
    null_sd <- sqrt(solve(t(Dmat) %*% Dmat )["Y","Y"])
  }


  ##################################
  #Computing proxy covariance matrix
  ##################################
  if(verbose)
  {
   print("Empirical Covariance missing computing proxy covariance for simulation")
  }
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

  N=1000 #number of simulation used to compute the covariance matrix
  #Generate random signal
  SNP=lev_res*500
  loci<- matrix(runif(N*SNP,min=0,2),ncol=N)
  bp= 1:SNP
  #Preparing for the wavelet transform
  Time01 <- (bp- min(bp))/(max(bp)-min(bp))
  Gen_W_trans <- apply(loci,2,my_wavproc)
  Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)
  #Compute proxy for empirical covariance matrix
  emp_cov <- (cov(Gen_W_trans))*(null_sd^2)
  if(verbose)
  {
    print("Proxy Covariance computed")
  }



  ######################
  #Set up for simulation
  ######################

  Pi_nt <- list()
  alt_sd <- 100*null_sd
  if(missing(base_shrink))
  {
    alp <-  1/sqrt(2*log(length(Y)))
  }
  else
  {
    alp <-  base_shrink
  }
  if(verbose)
  {
   print("Simulation of test statistics")
  }
  temp <- seq(from=1,to=size,by=size/10)[-1]-1
  out <- list()
  my_f <- function(y)
  {
    pis <- max_EM_post_Beta(y, lev_res = lev_res , null_sd = null_sd ,alt_sd = alt_sd,alp=alp)[[1]]
    return(pis)
  }
  if(verbose)
  {

    for (i in 1:10)
    {
      y <- rmvnorm(floor(size/10),mean=rep(0,nbeta),sigma =emp_cov )
      out[[i]] <-t( apply(y,1,my_f))
      print(paste(i*floor(size/10), "simulations performed"))
    }

  }
  if(verbose==FALSE)
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
