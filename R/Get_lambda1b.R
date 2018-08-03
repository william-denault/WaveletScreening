
#'@title Alternative computation of lambda1 for BF null distribution based on the Spectra C++ library
#'@description Compute the lambda parameter of the Bayes factors null distribution. Using the rARPACK SVD R routine.
#'@param Y a vector of the variable of interest.
#'@param confounder the confounding matrix with same a number of line equal to the length of Y. The intercept should not be included, if missing will generate a intercept matrix.
#'@param sigma_b value of the prior used in the Wavelet screaming.
#'@export
#'@references Quan Zhou and Yongtao Guan, On the Null Distribution of Bayes Factors in linear Regression, Journal of the American Statistical Association, 518, 2017.
#'@details Compute the lambda parameter of the Bayes factors distribution using its closed form which is provided in the paper of Quan Zhou and Yongtao Guan.
#'@seealso \code{\link{get_ncp}}
#'@examples \dontrun{
#'library(rARPACK)
#'Y <- rnorm(n=1000)
#'sigma <-0.2
#'get_lambda1B(Y=Y,sigma_b=sigma)
#'
#'}
get_lambda1B <- function(Y, confounder, sigma_b)
{
  Y <- as.vector(Y)

  # INPUT CHECKS
  print("Input dimensions:")
  if(!is.numeric(Y) || length(Y)==0){
    stop("ERROR: Y is not a numeric vector")
  } else {
    print(sprintf("%i phenotypes detected", length(Y)))
    if(all(Y %in% c(0, 1, NA))){
      print("Binary phenotype detected")
    } else if(!is.vector(Y)){
      stop("ERROR: Y is not a vector. Multi-phenotype analysis not implemented yet.")
    } else {
      print("Continuous phenotype detected")
    }
  }
  if(missing(confounder)) {
    print("no covariates provided, using intercept only")
    confounder <- data.frame(confounding=rep(1,length(Y)) )
  } else if(nrow(confounder)!=length(Y)) {
    stop("ERROR: number of samples in Y and confounder does not match")
  } else {
    print(sprintf("%i covariates for %i samples detected", ncol(confounder), nrow(confounder)))
    confounder <- cbind(rep(1,length(Y)),confounder)
  }

  # Clean missing samples from all inputs
  keepY <- complete.cases(Y)
  keepC <- complete.cases(confounder)
  nonmissing_index <- which(keepY & keepC)
  if(length(nonmissing_index) != length(Y)){
    print(sprintf("Warning: %i individuals will be removed due to missingness",
                  length(Y) - length(nonmissing_index)))
  }

  Y <- Y[nonmissing_index]
  confounder <- confounder[nonmissing_index,]

  W <- as.matrix(confounder, ncol=ncol(confounder))
  L <- as.matrix(Y , ncol=ncol(Y)) #reversed regression
  n = nrow(W)
  q = ncol(W)
  L <- as.matrix(L,ncol=1)
  PW = diag(n) - W %*% solve(t(W) %*% W) %*% t(W)
  p=1 #1df test
  X= PW %*% L
  HB = X %*% solve(t(X) %*% X + diag(1/sigma_b/sigma_b,p)) %*% t(X)

  lambda <- rARPACK::svds(HB, k=1)$d[1]

  return(lambda)
}
