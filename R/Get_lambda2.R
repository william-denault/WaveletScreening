get_lambda2 <- function(Y, confounder, sigma_b)
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
