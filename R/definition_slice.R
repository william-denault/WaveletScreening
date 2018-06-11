#'@title Define chromosomal regions for screening
#'@description  Define overlapping loci starting and ending positions for Wavelet screaming analysis
#'@param bp vector of the observed based pair positions in a chromosome
#'@param Loci_size size of the defined loci, limited by thresh size gaps on ends. Slices smaller than Loci_size will be skipped.
#'@param thresh maximal distance between two SNP within a loci. E.g 10000
#'@param Chr the chromosome's number where the slicing is made. By default set as NA
#'@export
#'@examples \dontrun{
#'Loci_size=1000000
#'thresh=10000
#'temp <- runif(n = 50000,min=1,max=10050)
#'for (i in 2:length(temp))
#'{
#'  temp[i] <- temp[i]+temp[i-1]
#'}
#'bp <- temp
#'df <-slice_definition(bp=bp,Loci_size=Loci_size,thresh = thresh ,Chr=5)
#'head(df)
#'}

slice_definition <- function(bp,Loci_size=1e6,thresh=1e4,Chr=NA)
{
  if(!is.vector(bp) | !is.numeric(bp)){
  	stop("ERROR: bp was not a numeric vector")
  }
  bp <- sort(bp)
  
  #First SNp at distance 0 of itself
  espacement <- c(0,diff(bp))
  
  ##################################
  ##location where spacing is to big
  ##################################

  my_index <- c(1,which(espacement > thresh) )
  if(length(my_index)==1)
  {
    my_index <- c(1, length(bp))
  }
  
  #Check distance between problematic spacing
  Width_loci <- rep(0,length(my_index)-1)


  #######################################
  #Simple window sliding half overlapping
  #######################################

  ##Strategy: if region to run analysis superior of 1.5 loci definition then
  #run half overlapping analysis

  #if lower than 1.5 than loci_size then run one analysis from the start one from the end

  #df output to define the extraction
  df <- data.frame(Chr= numeric(),posStart= numeric(),posEnd= numeric())

  for(i in 1:(length(my_index)-1))
  {
    my_diff <- bp[my_index][i+1]-bp[my_index][i]
 
    if( my_diff >= Loci_size ) #True means ok to run a wavelet analysis between my_index[i] and my_index[i+1]
    {
      Width_loci[i] <- my_diff

      if(my_diff -1>= 1.5*Loci_size)
      {
        temp1 <-0
        while(temp1 + Loci_size < my_diff)
        {
	  # add one slice if there's at least one SNP inside region
          # (-1 +1 just to ensure that we keep edge SNPs)
          my_loci <- data.frame(Chr, posStart=bp[my_index][i]+temp1-1, posEnd=bp[my_index][i]+temp1+Loci_size+1)
          if(any(bp >= my_loci[,2] & bp <= my_loci[,3])){
            df <- rbind(df,my_loci)
	  }

          temp1 <- temp1 +Loci_size/2
        }

        #to get the last part which is the rest of width_loci[i]/1.5*loci_size
        my_loci <- data.frame(Chr, posStart=bp[my_index][i+1]-Loci_size-1, posEnd=bp[my_index][i+1]+1 )
        if(any(bp >= my_loci[,2] & bp <= my_loci[,3])){
	  df <- rbind(df,my_loci)
	}
      }
      else
      {
        my_loci <- data.frame(Chr, posStart=bp[my_index][i]-1, posEnd=bp[my_index][i] + Loci_size + 1 )
        if(any(bp >= my_loci[,2] & bp <= my_loci[,3])){
	  df <- rbind(df,my_loci)
	}
        my_loci <- data.frame(Chr, posStart=bp[my_index][i+1] - Loci_size-1, posEnd=bp[my_index][i+1] + 1 )
        if(any(bp >= my_loci[,2] & bp <= my_loci[,3])){
	  df <- rbind(df,my_loci)
	}
      }
    }
  }

  # percentage of chromosome covered
  # note: does not necessarily correspond to coverage of genotyped positions
  coverage_of_analysis <- 100*sum(Width_loci)/(max(bp)-min(bp))
  
  colnames(df) <- c("Chr","posStart","posEnd")
  print(paste("number of slices defined: ", nrow(df)))
  print(paste("percentage of chromosome covered after slice definition: ",
              coverage_of_analysis,"%"))
  return(df)
}
