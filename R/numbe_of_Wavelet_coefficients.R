#'@title Number of wavelet coefficients for a given level of resolution
#'@description  Compute the number of wavelet coefficients up to a certain levl of resolution
#'@param lev_res the maximum level of resolution needed
#'@examples \dontrun{
#'n_coef_wc(7)
#'}




n_coef_wc <- function(lev_res)
{
  temp <- c()
  for( i in 0:lev_res)
  {
    temp <- c(temp,2^i)
  }
  sum(temp)
}
