#'@title Compute the relative size of a sub region.
#'@description  Compute the relative size of a sub region.
#'@param sub Output of extract_tree.
#'@return A proportion, used for multiple testing correction for the region






prop_sub <- function(sub)
{
  my_s <- as.numeric(gsub("^[^_]*_|_[^_]*$", "", colnames(sub )))
  out <- length(which(my_s ==min(my_s)))/   2^(min(my_s))
  return(out)

}
