#'@title Internal EM procedure
#'@description  Internal EM procedure
#'@param my_bayes a vector of Bayes Factors of length 2^lev_res
#'@param lev_res the level of resolution in the wavelet transform




max_EM_Lambda <- function(my_bayes,lev_res)
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
