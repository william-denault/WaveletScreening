#'@title Data visualization for the output of Wavelet Screaming
#'@description  Data visualization of Wavelet Screaming results.
#'@param res Output of Wavelet_screaming, without Bayes Factor.
#'@param bp a vector of the basepair position of the locus. You can provide only the starting and the endpoint of the locus. If missing ,set as 0, 1.
#'@param lev_res the maximum level of resolution needed. This has to be less or equal to the requested level of resolution in the Wavelet_screaming.
#'@param fill logical. If not, provide set as TRUE.
#'@param dg numerical. This is the number of digits displayed on the x-axis. If missing, set at 3.
#'@return return a ggplot
#' @details The function generates a ggplot from the Wavelet Screaming output. It represents the Betas for the different scales of the wavelet decomposition.
#'The size and the darkness of the points that represent the Betas are scaled by the value of the Betas.
#'If a Beta is not thresholded, then the region that represents the Beta is highlighted in order to give an overview of the size and the origin of the genetic signal.
#'@seealso \code{\link{Wavelet_screaming}}
#'
plot_WS <- function(res,bp,lev_res,fill,dg)
{

  if(missing(fill))
  {
    fill=TRUE
  }
  if(missing(dg))
  {
    dg=3
  }
  if(missing(bp))
  {
    bp=c(0,1)
  }
  sel <- names(res)
  pos <- grep("Pi",sel)
  res <- as.numeric(res)
  prt <- res[pos]#posterior probabilities of H1

  ypos <- seq(from = 1, to = 0, length.out=lev_res+1 )
  y<- c()
  x <-c()
  levres <- c()
  for( i in 0:lev_res)
  {

    xpos <-  seq(1,by=2,len=2^i)/2^(i+1)

    x <- c(x,xpos)

    for (j in (2^i):(2^(i+1)-1))
    {
      y[j] <- ypos[i+1]
      levres <- c(levres,i)
    }

  }
  point_size <- res[grep("Beta",sel)]#Betas for size of point
  point_size[which(prt== 0 )]<-NA

  #############
  #Option here
  #############
  disp <- c(point_size,point_size)

  point_size<-res[grep("Beta",sel)]
  df <- data.frame(x=x,y=y,levres=levres, ps =point_size)


  xstart <- df$x -1/(2^(df$levres+1))
  xend <- df$x +1/(2^(df$levres+1))

  group <- as.factor(1:dim(df)[1])

  xfil <- c(xstart,xend)
  ########################
  #avoid to display BF<1
  ########################
  xfil[which(is.na(disp))] <-NA

  gl <- c(group,group)
  my=c(df$y,df$y)+0.1

  mx <- c(x,x)
  mycol <- c(df$levres,df$levres)
  ps <-abs(c(point_size ,point_size ))
  ps2 <- c(point_size ,point_size )
  df_fill <- data.frame(xfil=xfil ,y=my,group=gl,col =mycol,x=mx, ps =abs(ps))

  lim <- c( min(-max(point_size),min(point_size)),
            max(max(point_size),-min(point_size)))
  lim


  if(fill==TRUE)
  {

    P1 <- ggplot(df_fill, aes(x=xfil,y=y, group=group,fill=col))+
      geom_area(position="identity")+
      geom_point(aes(x=x,y=my,size=ps+1,col=ps2))+

      #scale_color_gradient(low = 'grey80', high ='grey20',guide='none' )+
      scale_color_gradient2(low = 'blue', high ='red', mid="white",guide='none',limits=lim )+


      scale_fill_gradient( high= "green", low = "purple" ,guide='none')+
      guides( size = FALSE,col=FALSE)+
      ylab("Level of resolution")+
      xlab("Base pair position")+
      scale_y_continuous(breaks=unique(df_fill$y ), labels =0:lev_res)+
      scale_x_continuous(breaks=seq(0,1 ,by=0.125), labels =format(seq(min(bp),max(bp),length.out = 9), digits =dg, scientific = TRUE))+
      theme_bw()


  }
  if(fill==FALSE)
  {
    P1 <-     ggplot(df_fill, aes(x=x,y=my,size=abs(ps)))+
      guides( size = FALSE,col=FALSE)+
      geom_point(aes(x=x,y=my,size=ps+1,col=ps))+
      scale_color_gradient(low = 'grey80', high ='grey20',guide='none' )+
      scale_y_continuous(breaks=unique(df_fill$y ), labels =0:lev_res)+
      scale_x_continuous(breaks=seq(0,1 ,by=0.125), labels =seq(min(bp),max(bp),length.out = 9))+
      theme_bw()

  }

  return(P1)

}


