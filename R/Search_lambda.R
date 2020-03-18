#'@title Calibration of the hyperparameter
#'@description  Perform a grid search for the tuning hyperparameter
#'@param Sim a data frame with two columns named L_h and min_ph_pv, normally generated from the Sim_null function. At least 10000 simulation
#'@param lev_res the level of resolution in the wavelet transform
#'@param emp_cov Empirical covariance matrix of the beta values. It can be computed using several results of the wavelet screening using the betas values for different loci. If missing the function computes an approximation of the covariance matrix, this leads to a loss of power and a more conservative test statistics.
#'@param size number of simulation to be performed
#'@param sigma_b the parameter of the NIG prior used for the Betas computation.
#'@return The simulation under the null of the two test statistics used to build the final test (i.e., L_h and min(ph,pv))
#'@examples \dontrun{
#'Y <- rnorm(4000)
#'Sim <- Simu_null(Y=Y,lev_res = 6,sigma_b = 0.2,size=100000)
#'lambda <-Search_lambda(Sim,plot=TRUE)
#'
#' par(mfrow=c(2,1))
#'pen <- i*0
#'t1 <- lh+pen*(mph_pv)
#'muv <- median(t1,na.rm = TRUE)
#'sdv <- mad(t1,na.rm = TRUE)
#' hist(t1,nclass= 100,
#'          main=latex2exp("Histogramm of null using L_h"),xlab = "L_h")
#'pen <- j
#'t1 <- lh+pen*(mph_pv)
#'muv <- median(t1,na.rm = TRUE)
#'sdv <- mad(t1,na.rm = TRUE)
#'h <- hist(t1,nclass= 100,
#'          main=latex2exp("Histogramm of null using T_{S_{lambda}}"),xlab = " T_{S_{lambda}}")
#'par(mfrow=c(1,1))
#'}




Search_lambda <- function(Sim,plot=FALSE)
{
 if(dim(Sim)[1]<10000)
 {
   stop("Not enough simulations, please provide at least 10,000 simulations")
 }
  lh <- Sim[,"L_h"]
  mph_pv <- Sim[,"min_ph_pv"]
  #Starting rank
  muv <- median( lh,na.rm = TRUE)
  sdv <- mad( lh,na.rm = TRUE)



  if(dim(Sim)[1]<50000)
  {
    n_class <- 100
  }
  else
  {
    n_class <-1000
  }
  h <- hist(1-pnorm(lh,mean=muv,sd=sdv),nclass= n_class,plot = FALSE)
  x <- h$counts
  Rank0 <- rank(x)[1]


  hyp <- c()
  pos_check <- c()
  i <-0
  temp  <-0
  print("Rought search")
  while( temp < Rank0+1)
  {
    i <-i +1
    pen <- i*1000
    t1 <- lh+pen*(mph_pv)
    muv <- median(t1,na.rm = TRUE)
    sdv <- mad(t1,na.rm = TRUE)
    h <- hist(1-pnorm(t1,mean=muv,sd=sdv),nclass= n_class,plot = FALSE)
    x <- h$counts
    temp <- rank(x)[1]
    hyp <- c( hyp,temp )
    pos_check <- c(pos_check,i)

    if(i >10000)
    {
      print("Penalization parameter over 10 millions, provide bigger simulation set (suggested size 50k simulation)")
      break
    }
  }
  print("Dichotomised search")
  upper_lim <- i*1000
  lower_lim <- (i-1)*1000
  for(j in lower_lim:upper_lim)
  {

    pen <- j
    t1 <- lh+pen*(mph_pv)
    muv <- median(t1,na.rm = TRUE)
    sdv <- mad(t1,na.rm = TRUE)
    h <- hist(1-pnorm(t1,mean=muv,sd=sdv),nclass= n_class,plot = FALSE)
    x <- h$counts
    x <- rank(x)[1]
    if(x>(Rank0+1))
    {
      break
    }

  }

  if(plot==TRUE)
  {

    par(mfrow=c(2,1))
    pen <- i*0
    t1 <- lh+pen*(mph_pv)
    muv <- median(t1,na.rm = TRUE)
    sdv <- mad(t1,na.rm = TRUE)
    h <- hist(1-pnorm(t1,mean=muv,sd=sdv),nclass= n_class,
              main=c("Histogram of null p-values using L_h"),xlab = "p value",xlim = c(0,0.1))
    pen <- j
    t1 <- lh+pen*(mph_pv)
    muv <- median(t1,na.rm = TRUE)
    sdv <- mad(t1,na.rm = TRUE)
    h <- hist(1-pnorm(t1,mean=muv,sd=sdv),nclass= n_class,
              main=paste("Histogram of null p-values using T_",j),xlab = "p value",xlim = c(0,0.1))
    par(mfrow=c(1,1))
  }

out <- j-1
return(out)

}
