#'@title Calibration of the hyperparameter
#'@description  Perform a grid search for the tuning hyperparameter
#'@param Sim a data frame with two columns named L_h and min_ph_pv, normaly generated from the Sim_null function. At least 10000 simulation
#'@param lev_res the level of resolution in the wavelet transform
#'@param emp_cov Emprical covariance matrix of the beta values. Can be computed using several results of the wavelet screaming using the betas values for different loci. If missing the function compute an approximation of the covariance matrix, this leads to a lost of power and a more conservative test statitics.
#'@param size number of simulation to be performed
#'@param sigma_b the parameter of the NIG prior used for the Betas computation.
#'@retrun The simulation under the null of the two test statistics used to build the final test (i.e L_h and min(ph,pv))





Search_lambda(Sim)
{
 if(dim(Sim)[1]<10000)
 {
   print("Not enough simulation, please provide at least 10,000 simulations")
   break
 }
  lh <- Sim[,"L_h"]
  mph_pv <- Sim[,"min_ph_pv"]
  #Starting rank
  muv <- median( lh,na.rm = TRUE)
  sdv <- mad( lh,na.rm = TRUE)
  n_class <-1000
  h <- hist(1-pnorm(df0$t1[c(1000:dim(df0)[1])],mean=muv,sd=sdv),nclass= n_class,plot = FALSE)
  x <- h$counts
  Rank0 <- rank(x)[1]



  for ( i in (length(hyp)+1):1000)
  {
    pen <- i*100
    df0$t1 <- df0$lh+pen*(mph_pv)
    muv <- median(df0$t1[c(1000:dim(df0)[1])],na.rm = TRUE)
    sdv <- mad(df0$t1[c(1000:dim(df0)[1])],na.rm = TRUE)
    n_class <-1000
    h <- hist(1-pnorm(df0$t1[c(1000:dim(df0)[1])],mean=muv,sd=sdv),nclass= n_class,plot = FALSE)
    #h <- hist(1-pnorm(df0$t1[c(1000:dim(df0)[1])],mean=muv,sd=sdv),nclass= n_class)

    x <- h$counts
    chec <- rank(x)[1]

    hyp <- c( hyp, chec)
    print(i)
  }
  hyp <- c()
  pos_check <- c()
  i <-0
  print("Rought search")
  while( temp < Rank0+1)
  {
    i <-i +1
    pen <- i*100
    t1 <- lh+pen*(mph_pv)
    muv <- median(t1,na.rm = TRUE)
    sdv <- mad(t1,na.rm = TRUE)
    n_class <-1000
    h <- hist(1-pnorm(t1,mean=muv,sd=sdv),nclass= n_class,plot = FALSE)
    x <- h$counts
    chec <- rank(x)[1]
    hyp <- c( hyp, chec)
    pos_check <- c(pos_check,i)

    if(i >10000)
    {
      break
    }
  }
  upper_lim <- i*100
  lower_lim <- (i-1)*100




}
