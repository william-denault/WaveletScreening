#'@title Internal likelyhood ratio
#'@description   internal function objective function.
#'@param my_pi a vector of the proportion of association per level of resolution.
#'@param my_bayes a vector of bayes factor for all the wavelet transform component.




Lambda_stat <- function (my_pi, my_bayes)
{
  # vector: pi1 pi2 pi2 pi3 pi3 pi3 pi3...
  my_pi_vec = rep(my_pi, 2^(1:length(my_pi)-1))
  coefs = 1-my_pi_vec + my_pi_vec * my_bayes[1:(2^length(my_pi)-1)]
  prod(coefs)
}
