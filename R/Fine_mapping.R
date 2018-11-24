
#'@title Defintion of the regions with Bayes factor over a threshold.
#'@description  Fine mapping tool for output of the Wavelet Screaming function.
#'@param res Output of Wavelet_screaming, without Betas.
#'@param lev_res the maximum level of resolution of the previous analysis.
#'@param thresh numeric, Bayes factor threshold to defined the fine mapping. If missing set as 1.
#'@param start numeric, start in base pair of the analyzed regions .
#'@param chr numeric, end in base pair of the analyzed regions .
#'@details return a list of chr, start, end position, that correspond of the sub regions defined by the dyadic decomposition of wavelet that are associated with a Bayes factor over the defined threshold.


fine_map <- function(res,lev_res,thresh,start,end,chr)
{
  if(missing(thresh))
  {
    thresh=1

  }
  k <-1
  temp <- lev_res+2
  l1 <- which(res[- c(1:temp )] >thresh )
  fmap <- list()
  for (i in 1: lev_res)
  {
    lt <- l1[which(l1 < sum(2^(0:i)))]
    lt <-  lt[which(lt> sum(2^(0:(i-1 ) )  ) ) ]



    if(length(lt) >0 )
    {
      ltt <- ( lt- sum(2^(0:(i-1 ) ) ) )
      for ( j in 1: length(ltt))
      {
        startpos <- as.numeric(start)-1 +
          (ltt[j] -1)*(1/(2^i)) *(as.numeric(end) -as.numeric(start) )
        endpos <- as.numeric(start)-1 +
          (ltt[j] )*(1/(2^i)) *(as.numeric(end) -as.numeric(start)  )

        chr <- chr

        fmap[[k]] <-  c( chr, startpos,endpos)

        k <- k+1
      }



    }



  }
  return(fmap)
}
