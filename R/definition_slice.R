#'@title defined the chromosal regions for screening
#'@description  Define overlapping loci starting and ending positions for Wavelet screaming analysis
#'@param bp vector of the observed based pair in a Chromosome
#'@param Loci_size size of the defined loci
#'@param thresh maximal distance between two SNP within a loci. E.g 10000
#'@param Chr the chromosome's number  where the slicing is made. By default set as NA
#'@export
#'@examples \dontrun{
#'Loci_size =1000000
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

slice_definition <- function(bp,Loci_size,thresh,Chr)
{
  if(missing(Chr))
  {
    Chr <-NA
  }

  #First SNp at distance 0 of itself
  espacement <- c(0,rep(NA,length(bp)[1]-1))


  for(i in 1:(length(espacement)-1 ) )
  {
    espacement[i+1] <-bp[i+1] -bp[i]
  }

  ##################################
  ##location where spacing is to big
  ##################################

  my_index <- c(1,which(espacement > thresh) )
  if(length(my_index)==1)
  {
    my_index <- c(1, length(bp))
  }
  #Check distance between problematic spacing

  possible_loci <- rep(FALSE,length(bp[my_index]) )
  Width_loci <- rep(FALSE,length(bp[my_index]) )

  for( i in 1: (length(possible_loci)-1) )
  {
    my_diff <- bp[my_index][i+1]-bp[my_index][i]
    if( (  Loci_size <= my_diff) )
    {
      possible_loci[i] <- TRUE #Truemean ok to run a wavelet analysis between my_index[i] and my_index[i+1]
      Width_loci[i] <- my_diff
    }
  }


  #percentage of included genotyped data
  coverage_of_analysis <- 100*sum(Width_loci)/(max(bp)-min(bp))


  #######################################
  #Simple window sliding half overlapping
  #######################################

  ##Strategy: if region to run analysis superior of 1.5 loci definition then
  #run half overlapping analysis

  #if lower than 1.5 than loci_size then run one analysis from the start one from the end

  #df output to define the extraction
  df <- data.frame(Chr= numeric(),poStart= numeric(),posEnd= numeric())

  for(i in 1:(length(possible_loci)-1))
  {

    if(possible_loci[i] ==TRUE)
    {
      if(Width_loci[i] >= 1.5*Loci_size)
      {


        temp1 <-0
        while(temp1  + Loci_size < Width_loci[i])
        {

          #definition of one slide with a "dense" enoguh region (-1 +1 just to insure that we keep the SNPs)
          my_loci <- c(Chr,bp[my_index][i] +temp1-1 ,bp[my_index][i]+ temp1 +Loci_size +1 )
          df <- rbind(df,my_loci)

          temp1 <- temp1 +Loci_size/2

        }

        #to get the last part whihc is the rest of width_loci[i]/1.5*loci_size
        my_loci <- c(Chr,bp[my_index][i+1] -Loci_size-1 ,bp[my_index][i+1]+1 )
        df <- rbind(df,my_loci)

      }
      else
      {
        my_loci <- c(Chr,bp[my_index][i]-1 ,bp[my_index][i] + Loci_size + 1 )
        df <- rbind(df,my_loci)
        my_loci <- c(Chr,bp[my_index][i+1] - Loci_size-1 ,bp[my_index][i+1] + 1 )
        df <- rbind(df,my_loci)

      }






    }


  }

  colnames(df) <- c("Chr","poStart","posEnd")
  print(paste("percentage of included of the SNPs included after slice definition: ",
              coverage_of_analysis,"%"))
  return(df)
}
