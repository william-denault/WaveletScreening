#'@title Compute posterior variance Beta value from Bayesian linear regression
#'@description  Compute posterior variance Beta value from Bayesian linear regression
#'@param y phenotype vector/variable of interest, has to be numeric.
#'@param confounder the confounding matrix with the same sample order as Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param sigma_b the parameter of the NIG prior used for the Bayes Factor computation. We advised to set this value between 0.1 and 0.2
#'@param all logical, if set as TRUE return all the Beta value (including the ones form the confounding factors). If set as FALSE only return the estimate for x, set as FALSE if missing.
#'@details The Wavelet_screaming performed reverse regression so variance for all posterior distribution are equal.
#'@return A matrix variance covariance matrix
#'@examples \dontrun{
#'
#'x <- rnorm(1000)
#'y <- x+3
#'sigma_b <- 0.2
#'Cov_Beta(y=y,sigma_b = sigma_b,all=TRUE)
#'}




Cov_Beta <- function(y,confounder, sigma_b,all)
{
  if(missing(confounder) )
  {
    confounder <- data.frame( confounder=rep(1,length(y)) )
  }
  else{
    confounder <- data.frame(confounder)
  }
  if(missing(all))
  {
    all=FALSE
  }

  pc <- dim(confounder)[2]
  Dmat <- cbind(confounder,x)
  Dmat <- as.matrix(Dmat)
  index <- pc+1
  resM <- resCI <- (1/sigma_b/sigma_b)*solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,
                                                                      dim(Dmat)[2]))
  #postvar <- diag( solve(t(Dmat) %*% Dmat + diag(1/sigma_b/sigma_b,2)))*1/sigma_b/sigma_b
  if(all==FALSE)
  {
    index <- pc + 1
    return(resM[index, index])
  }
  if(all==TRUE)
  {
    return(resM)
  }

}

