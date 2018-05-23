#'@title Simulation of Lambda null distribution
#'@description  Simulation of Lambda null distribution.
#'@param nsimu number of simulation to perform.
#'@param lambda the lambda parameter of the null distribution of the Bayes Factor.
#'@param lev_res the maximum level of resolution needed.
#'@param ncp the lambda parameter of the null distribution of the Bayes Factor, if not specified set at 0.
#'@details Use the theoretical null distribution of the generated Bayes factor to perform simulation of null distriution of the test statistics of the Wavelet screaming procedure.
#'@references Quan Zhou and Yongtao Guan, On the Null Distribution of Bayes Factors in linear Regression, Journal of the American Statistical Association, 518, 2017.
#'@seealso \code{\link{Wavelet_screaming}}
#'@export
#'@examples \dontrun{
#'set.seed(1)
#'Simu_Lambda_null(100,0.7,2)
#'}














Simu_Lambda_null <- function(nsimu,lambda,lev_res,ncp)

{
  if(missing(ncp)){
    ncp=0
  }


  sumlog <- function (A1, A2)
  {

    if(A1 > A2){
      res = A1 + log(1 + exp(A2 - A1))
    }else{
      res = A2 + log(exp(A1 - A2) + 1)
    }

    return (res)

  }
  n_coef_wc <- function(lev_res)
  {
    temp <- c()
    for( i in 0:lev_res)
    {
      temp <- c(temp,2^i)
    }
    sum(temp)
  }
  max_EM_Lambda <- function(my_bayes)
  {
    niter=10000
    epsilon <- 10^-4
    p_vec <- c()
    for(gi  in  0: lev_res)
    {
      # EM algorithm for each group separately


      N_obllikli = 0
      logpi = log(0.5)
      pi <-0.5
      log1pi = logpi


      pp = 0
      logPiBF =   log(my_bayes[(2^gi):(2^(gi+1)-1)]) + logpi
      logden <- c()
      for (i in 1:length(logPiBF))
      {
        logden[i] <- sumlog(logPiBF[i],log1pi)
      }
      pp = pp+sum(exp(logPiBF - logden))
      N_obllikli = sum(logden)
      O_obllikli = N_obllikli


      for( iter in  0:niter){

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
  nBF <- n_coef_wc(lev_res)
  BF <- matrix(NA,ncol=nBF,nrow=nsimu)


  print("Simulation Bayes Factors")
    for (i in  1:nsimu)

    {

      for (j in 1: nBF)
      {
        BF[i,j] <- exp( (lambda*rchisq( n=1 ,ncp=ncp, df=1)+log(1-lambda) )/2 )
      }
    }
  Lambda_stat <- function (my_pi, my_bayes)
  {
    # vector: pi1 pi2 pi2 pi3 pi3 pi3 pi3...
    my_pi_vec = rep(my_pi, 2^(1:length(my_pi)-1))
    coefs = 1-my_pi_vec + my_pi_vec * my_bayes[1:(2^length(my_pi)-1)]
    prod(coefs)
  }

  print("Lambda maximisation")
    my_pis<- t(apply(BF,1,max_EM_Lambda))

    Simu_Lambda <- c()
    for (i in 1:dim(my_pis)[1])
    {
      Simu_Lambda <-c( Simu_Lambda , Lambda_stat(my_bayes=BF[i,],my_pi=my_pis[i,]))

    }


    return(Simu_Lambda)

}

