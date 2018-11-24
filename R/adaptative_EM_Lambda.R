#'@title  EM procedure for sub tree
#'@description  EM procedure for sub tree
#'@param sub a Output of extract_tree
#'@return Estimated proportion for the sub tree
#'@examples \dontrun{
#'
#'#using res for the Wavelet_screaming exemple
#'
#'
#'sub_analysis <- function(res, lev_res )
#'{
#'  sub <- extract_tree(res,lev_res=lev_res)
#'  my_pi <- adaptative_EM_Lambda(sub)
#'  out <-  adaptative_Lambda (my_pi, sub)
#'  return(out)
#'}
#'
#'
#'sub_analysis(res, 6)
#'
#'}




adaptative_EM_Lambda <- function(sub)
{
  niter=10000
  epsilon <- 10^-4
  p_vec <- c()

  my_bayes <- as.numeric(sub)
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

