#'@title Extract sub result from result of the wavelet_screaming
#'@description  Function to performed zoomed analysis of the wavelet screaming output
#'@param res Output of Wavelet_screaming.
#'@param lev_res the maximum level of resolution needed, has to be less or equal to the request level of resolution in the Wavelet_screaming.
#'@param thresh Minimal value of the Bayes Factor to  defined the a sub region, if missing set as 1.
#'@return A vector correpsonding of the sub tree for the zoomed analysis.
#'@examples \dontrun{
#'
#'#using res for the Wavelet_screaming exemple
#'
#'
#'sub_analysis <- function(res, lev_res )
#'{
#'  sub <- extract_tree(res,lev_res=lev_res)
#'  my_pi <- adaptative_EM_Lambda(sub)
#'  out <-  adaptative_Lambda (my_pi, sub)
#'  return(out)
#'}
#'
#'
#'sub_analysis(res, 6)
#'
#'}



extract_tree <- function(res,lev_res,thresh)
{

  if(missing(thresh))
  {
    thresh <- 1
  }
  res <- res[-c(1:(lev_res+2))]
  temp <- colnames(res)[which(as.numeric(res)>thresh)]

  #get the Scale
  my_s <- as.numeric(gsub("^[^_]*_|_[^_]*$", "", temp))

  cor <- c()
  for (i in 1:length(my_s))
  {
    temp <- sum(2^(0:(my_s[i]-1)))
    cor <- c(cor, temp)
  }
  #Get the location of the BF over thresh
  my_l <-   which(as.numeric(res)>thresh)-cor

  #check if the lower scale are nested into the upper scale

  temp <- which(my_s >min(my_s))

  #Case where ther is only one Bayes Factor above 1/ all at the same scale
  if(length(temp)>0)
  {
    temps <- my_s[temp] - min(my_s)

    #my_l essemble of the location needed at the min scale
    cors <- c()
    for( i in 1: length(temp))
    {
      temp1 <- my_l[temp[i]]

      #divided enought time to get the corresponding scale location
      for(j in 1: temps[i])
      {
        temp1 <- ceiling(temp1/2)
      }
      cors <- c(cors,temp1)
    }
    #Corresponding WC at min scale
    locations <- c( my_l[which(my_s==min(my_s))],  cors)
  }
  else{
    locations <- c( my_l[which(my_s==min(my_s))])
  }

  #Index to select
  my_index <- c()

  span <- min(locations):max(locations)
  for ( i in min(my_s):lev_res)
  {
    temp <- sum(2^( 0:( min(i) -1) ) )+ span
    span <- (2*min(span)-1):(2*max(span))
    my_index <- c( my_index,temp)

  }

  return(res[my_index])

}

