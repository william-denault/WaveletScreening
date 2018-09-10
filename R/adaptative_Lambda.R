#'@title Likelihood ratio for sub tree
#'@description   internal function objective function for sub trees.
#'@param my_pi a vector of the proportion of association per level of resolution.
#'@param sub Output of extract_treet.
#'@return Value of the likelihood on a sub tree.




adaptative_Lambda <- function(my_pi,sub)
{

  my_bayes <- as.numeric(sub)
  BF_class <-  as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))

  my_pi_vec <- c()
  for ( i in 1: length(unique(BF_class)))
  {
    temp <- rep(my_pi[i], c( length( which(BF_class == unique(BF_class)[i] ) ) ) )
    my_pi_vec <- c(my_pi_vec, temp )
  }
  coefs = 1-my_pi_vec + my_pi_vec * my_bayes

  prod(coefs)
}

