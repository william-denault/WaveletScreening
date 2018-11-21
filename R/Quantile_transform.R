#'@title Perform a quantile transform
#'@description  Perform a quantile transform with ties taken at random
#'@param x a vector to be quantile transform
#'@export
#'@examples \dontrun{
#'x <- rnorm(n=1000)
#'Quantile_transform(x)
#'}



Quantile_transform <- function(x)
{
  x.rank = rank(x, ties.method="average")
  return(qqnorm(x.rank,plot.it = F)$x)
}
