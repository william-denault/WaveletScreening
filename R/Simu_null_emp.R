#'@title Simulation of the null statistics
#'@description  Simulation of the null statistics
#'@param res an output of Wavelet_screening function. The user can also provide a matrix of results of Wavelet_screening (from the same analysis), where the results have been concatenated by row (rbind).
#'@param smp_size Sample size from the main run of Wavelet screening
#'@param lev_res the level of resolution in the wavelet transform
#'@param coeftype type of wavelet coefficient used for the screening (choice "c" or "d"). If missing set as "c"
#'@param size number of simulation to be performed. If not specified set at 10000
#'@param base_shrink numeric, value used in the thresholding of the proportion of assocation, if non specificed set up as 1/sqrt(2*log(sample_size)
#'@param print logical parameter set as TRUE, if TRUE sends a message when 10\% of the simulations have been completed.
#'@return The simulation under the null of the two test statistics used to build the final test (i.e., L_h and min(ph,pv))
Simu_null_emp <- function(res,
                          smp_size,
                          lev_res,
                          coeftype,
                          size= 10000,
                          base_shrink,
                          print= TRUE)
{

  #####################################
  #Size of the multivaraite to simulate
  #####################################
  nbeta <- sum(2^(0:lev_res))
  ######################################
  #compute theoretical variance of Betas
  ######################################

  if(missing(coeftype))
  {
    print( "missing coeftype set as d")
    coeftype <- "c"
  }

  #computing the structure of the covariance for the simulations
  my_wavproc <- function(y)
  {
    #Kovac and Silvermann 2000
    mygrid <- wavethresh::makegrid(t=Time01,y=y)
    LDIRWD <- irregwd(mygrid,filter.number=1)
    class(LDIRWD) <- "wd"

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

  N=1000
  #Generate random signal
  SNP=1:2^(lev_res+3)
  bp= SNP
  #Preparing for the wavelet transform
  Time01 <- (bp- min(bp))/(max(bp)-min(bp))
  loci<- matrix(runif(N*2^(lev_res+3),min=0,2),ncol=N)
  for ( i in 1:N)
  {
    loci[,i]  <- rnorm(n = 2^(lev_res+3))


  }

  bp=1:2^(lev_res)
  #Preparing for the wavelet transform
  Gen_W_trans <- apply(loci,2,my_wavproc)
  Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)
  #Compute proxy for empirical covariance matrix
  emp_cov <- (cov(Gen_W_trans))

  ################################################
  #Computing robustified variance for the diagonal
  ################################################
  if( length(dim(res) )==0) {
    var_sim <-  mad( res[ grep( pattern = "Beta", names(res))])^2
  }   else   {
    subres <- res [, grep( pattern = "Beta*", colnames(res))]
    temp <-  c(subres)
    var_sim <- mad(temp)^2
  }

  emp_cov <-  var_sim * (emp_cov/(max(emp_cov)))

  ######################
  #Set up for simulation
  ######################

  null_sd <-  mean(res[, grep(pattern = "null_sd_start_EM",colnames(res))])
  Pi_nt <- list()
  alt_sd <- 100*null_sd
  if(missing(base_shrink))
  {
    alp <-  1/sqrt(2*log(smp_size))
  }
  else
  {
    alp <-  base_shrink
  }

  print("Simulation of test statistics")
  temp <- seq(from=1,to=size,by=size/10)[-1]-1
  out <- list()
  my_f <- function(y)
  {
    pis <- max_EM_post_Beta(y, lev_res = lev_res , null_sd = null_sd ,alt_sd = alt_sd,alp=alp)[[1]]
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
