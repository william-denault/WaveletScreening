#'@title Compute lambda1 for BF null distribution
#'@description Compute the lambda parameter of the Bayes factors null distribution.
#'@param Y a vector of the variable of interest.
#'@param counfounder the counfounding matrix with same a number of line equal to the length of Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param sigma_b value of the prior used in the Wavelet screaming.
#'@export
#'@references Quan Zhou and Yongtao Guan, On the Null Distribution of Bayes Factors in linear Regression, Journal of the American Statistical Association, 518, 2017.
#'@details Compute the lambda parameter of the Bayes factors distribution using its closed form which is provided in the paper of Quan Zhou and Yongtao Guan.
#'@seealso \code{\link{get_ncp}}
#'@examples \dontrun{
#'Y <- rnorm(n=1000)
#'sigma <-0.2
#'get_lambda1(Y=Y,sigma_b=sigma)
#'}







get_lambda1 <- function(Y,counfounder,sigma_b)

{

  if(missing(counfounder))
  {
    counfounder <- data.frame(counfounding =rep(1,length(Y)) )
  } else
  {
    counfounder <- rbind(rep(1,length(Y)),counfounder)
  }
  W <- as.matrix(counfounder, ncol=ncol(counfounder))
  L <- as.matrix(Y , ncol=ncol(Y)) #reversed regression
  n = nrow(W)
  q = ncol(W)
  L <- as.matrix(L,ncol=1)
  PW = diag(n) - W %*% solve(t(W) %*% W) %*% t(W)
  p=1 #1df test
  X= PW %*% L
  HB = X %*% solve(t(X) %*% X + diag(1/sigma_b/sigma_b,p)) %*% t(X)

  temp <- svd(HB)
  lambda <- temp$d[1]

  return(lambda)
}








