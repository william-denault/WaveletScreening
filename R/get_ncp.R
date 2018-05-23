#'@title Compute ncp for BF null distribution
#'@description  Compute non central parameter of the Bayes factors distribution.
#'@param Y a vector of the variable of interest.
#'@param Counfounder the counfounding matrix with same a number of line equal to the length of Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param WCs the vector of a specific wavelet coefficient for each inidivual.
#'@param sigma_b value of the prior used in the Wavelet screaming.
#'@export
#'@references Quan Zhou and Yongtao Guan, On the Null Distribution of Bayes Factors in linear Regression, Journal of the American Statistical Association, 518, 2017.
#'@details Compute non central parameter of the Bayes factors distribution using its closed form which is provided in the paper of Quan Zhou and Yongtao Guan.
#'@seealso \code{\link{get_lambda1}}
#'@examples \dontrun{
#'Y <- rnorm(n=1000)
#'WCs <- rnorm(n=1000)
#'sigma <-0.2
#'get_ncp(Y=Y,WCs=WCs,sigma_b=sigma)
#'}




get_ncp <- function( Y,Counfounder, WCs,sigma_b ){

  #Y phenotype
  #Counfounder= Counfounding matrix
  #WCs= the wavelets coefficients

  if(missing(Counfounder))
  {
    Counfounder <- data.frame(Counfounding =rep(1,length(Y)) )
  } else
  {
    Counfounder <- rbind(rep(1,length(Y)),Counfounder)
  }
  my_Y <- Y
  W <- as.matrix(Counfounder, ncol=ncol(Counfounder))
  L <- as.matrix(my_Y , ncol=ncol(my_Y)) #reversed regression
  n = nrow(W)
  q = ncol(W)
  L <- as.matrix(L,ncol=1)
  y <-  as.matrix(WCs,ncol=1)

  p = 1
  PW = diag(n) - W %*% solve(t(W) %*% W) %*% t(W)
  X= PW %*% L
  HB = X %*% solve(t(X) %*% X + diag(1/sigma_b/sigma_b,p)) %*% t(X)
  log.R = -0.5*n*log(1 - (t(y) %*% HB %*% y) / (t(y) %*% PW %*% y ))
  delta = svd(X)$d
  lambda = delta^2 / (delta^2 + 1/sigma_b/sigma_b)
  log.T = sum(log(1-lambda))/2
  # one can check:
  # log.T = -0.5*log(det( t(X) %*% X * sigma_b * sigma_b + diag(p)))
  bf = exp(log.T + log.R)

  #Computation null distribution para
  temp <- svd(HB)
  lambda <- temp$d[1]
  u <- temp$u[,1]
  Beta <- solve(t(L) %*% L) %*% t(L)%*%y
  my_ncp = ((t(u)%*%L*Beta)^2)/n
  #print(c(bf,my_ncp,lambda))

  return(my_ncp)
}
