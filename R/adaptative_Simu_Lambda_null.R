#'@title Simulation of Lambda null distribution on a sub tree
#'@title Simulation of Lambda null distribution on a sub tree
#'@description  Simulation of Lambda null distribution.
#'@param nsimu number of simulation to perform.
#'@param lambda the lambda parameter of the null distribution of the Bayes Factor.
#'@param lev_res the maximum level of resolution needed.
#'@param ncp the lambda parameter of the null distribution of the Bayes Factor, if not specified set at 0.
#'@details Use the theoretical null distribution of the generated Bayes factor to perform simulation of null distriution of the test statistics of the Wavelet screaming procedure.
#'@references Quan Zhou and Yongtao Guan, On the Null Distribution of Bayes Factors in linear Regression, Journal of the American Statistical Association, 518, 2017.










adaptative_Simu_Lambda_null <- function(nsimu,lambda,lev_res,ncp,sub,thresh)

{

  if(missing(thresh))
  {
    thresh <-1
  }
  if(length( which( as.numeric(sub)>thresh ) )==1 )
  {
    print( "Only one Bayes Factor above the choosen threshold(set as one per deault), analytical formula available. Please use the function analytical_p")
    break
  }


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


  #Adapted version that match the sub structure
  loc_EM_Lambda <- function(my_bayes)
  {
    niter=10000
    epsilon <- 10^-4
    p_vec <- c()

    BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))
    for(gi in unique(BF_class))
    {
      # EM algorithm for each group separately

      N_obllikli = 0
      logpi = log(0.5)
      pi <- 0.5
      log1pi = logpi

      pp = 0
      logPiBF = log(my_bayes[which(BF_class==gi)]) + logpi
      logden <- c()
      for (i in 1:length(logPiBF))
      {
        logden[i] <- sumlog(logPiBF[i],log1pi)
      }
      pp = pp+sum(exp(logPiBF - logden))
      N_obllikli = sum(logden)
      O_obllikli = N_obllikli


      for(iter in  0:niter){
        pi = pp/(length(which(BF_class ==gi)))
        logpi  = log(pi)
        log1pi = log(1-pi)
        logPiBF =   log(my_bayes[which(BF_class==gi)]) + logpi
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

  nBF <- length(as.numeric(sub))
  BF <- matrix(NA,ncol=nBF,nrow=nsimu)


  print("Simulation Bayes Factors")
  for (i in  1:nsimu)

  {

    for (j in 1: nBF)
    {
      BF[i,j] <- exp( (lambda*rchisq( n=1 ,ncp=ncp, df=1)+log(1-lambda) )/2 )
    }
  }
  loc_Lambda_stat <- function (my_pi, my_bayes,sub)
  {
    BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))

    my_pi_vec <- c()
    for ( i in 1: length(unique(BF_class)))
    {
      temp <- rep(my_pi[i], c( length( which(BF_class == unique(BF_class)[i] ) ) ) )
      my_pi_vec <- c(my_pi_vec, temp )
    }
    coefs = 1-my_pi_vec + my_pi_vec * my_bayes

    prod(coefs)
  }

  print("Lambda maximisation")
  my_pis<- t(apply(BF,1,loc_EM_Lambda))

  Simu_Lambda <- c()
  for (i in 1:dim(my_pis)[1])
  {
    Simu_Lambda <-c( Simu_Lambda , Lambda_stat(my_bayes=BF[i,],my_pi=my_pis[i,]))

  }


  return(Simu_Lambda)





}





