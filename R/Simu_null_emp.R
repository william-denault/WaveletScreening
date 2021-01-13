
#'@title Simulation of the null statistics
#'@description  Simulation of the null statistics
#'@param res an output of Wavelet_screening function. The user can also provide a matrix of results of Wavelet_screening (from the same analysis), where the results have been concatenated by row (rbind).
#'@param smp_size Sample size from the main run of Wavelet screening
#'@param lev_res the level of resolution in the wavelet transform
#'@param coeftype type of wavelet coefficient used for the screening (choice "c" or "d"). If missing set as "c"
#'@param size number of simulation to be performed. If not specified set at 10000
#'@param base_shrink numeric, value used in the thresholding of the proportion of assocation, if non specificed set up as 1/sqrt(2*log(sample_size)
#'@param print logical parameter set as TRUE, if TRUE sends a message when 10\% of the simulations have been completed.
#'@return The simulation under the null of the two test statistics used to build the final test (i.e., L_h and min(ph,pv))

Simu_null_emp <- function(res,
                          smp_size,
                          lev_res,
                          coeftype,
                          size= 10000,
                          base_shrink,
                          print= TRUE)
{

  #####################################
  #Size of the multivaraite to simulate
  #####################################
  nbeta <- sum(2^(0:lev_res))
  ######################################
  #compute theoretical variance of Betas
  ######################################

  if(missing(coeftype))
  {
    print( "missing coeftype set as c")
    coeftype <- "c"
  }

  #computing the structure of the covariance for the simulations
  my_wavproc <- function(y)
  {
    #Kovac and Silvermann 2000
    mygrid <- wavethresh::makegrid(t=Time01,y=y)
    LDIRWD <- irregwd(mygrid,filter.number=1)
    class(LDIRWD) <- "wd"

    res <- c()
    for(i in 0: lev_res){
      if(coeftype == "d"){
        res <- c(res, accessD( LDIRWD,lev = i) )
      } else if (coeftype == "c") {
        res <- c(res, accessC( LDIRWD,lev = i) )
      } else {
        stop(paste("ERROR: coeftype", coeftype, "not recognized"))
      }
    }

    return(res)
  }



  max_EM_post_Beta <- function(my_betas, lev_res,null_sd,alt_sd,alp) {
    niter = 100
    epsilon <- 10^-4
    p_vec <- c()
    sd_vec <- c()
    erreur<-1+epsilon
    eps <-10^-10#slight correction in case of non identifiable mixture
    #prevent from having division by 0 in update parameter sum(temp)
    betasub = my_betas
    m0.hat<-0
    m1.hat<-0
    sigma0.hat<- null_sd
    sigma1.hat<-alt_sd
    #Prevent from label swapping
    if(sigma1.hat < sigma0.hat){
      sigma1.hat <- 3*sigma0.hat+sigma1.hat
    }

    p.hat<-0.25
    new.params<-c(m0.hat,m1.hat,sigma0.hat,sigma1.hat,p.hat)
    erreur<-1+epsilon
    iter <- 1

    while((erreur>epsilon)&(iter<=niter))
    {
      old.log.lik<- sum(log(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat)))
      old.params<-new.params
      #vecteur des pi_{i1}^{t}
      temp<-p.hat*dnorm( betasub ,m1.hat,sigma1.hat)/(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat))
      #Update parameter
      p.hat<-mean(temp)
      m1.hat<-sum(temp* betasub)/(sum(temp)+eps)
      #m0.hat<-sum((1-temp)* betasub)/(sum(1-temp)+eps)
      sigma1.hat<-sqrt( sum(temp*( betasub-m1.hat)^2)/(sum(temp)+eps) )+alt_sd
      sigma0.hat<-sqrt( sum((1-temp)*( betasub-m0.hat)^2)/(sum(1-temp)+eps) )
      #limit the decrease of sigma0.hat in case of non identifiable mixture
      if(sigma0.hat < 0.1*null_sd ){
        sigma0.hat <- 0.1*null_sd
      }
      new.params<-c(m0.hat,m1.hat,sigma0.hat,sigma1.hat,p.hat)
      #Check end
      new.log.lik<- sum(log(p.hat*dnorm( betasub ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( betasub ,m0.hat,sigma0.hat)))
      epsilon <- abs( new.log.lik -old.log.lik)
      iter<-iter+1

    }

    #Proba Belong belong to the alternative:
    pos.prob <- rep(NA,length(my_betas))
    for (i in 1:length(my_betas))
    {
      pos.prob[i] <- p.hat*dnorm( my_betas[i] ,m1.hat,sigma1.hat)/(p.hat*dnorm( my_betas[i] ,m1.hat,sigma1.hat)+(1-p.hat)*dnorm( my_betas[i] ,m0.hat,sigma0.hat))

    }

    lambcom <- rep(NA,(lev_res+1))
    p_vec   <- rep(NA,(lev_res+1))
    for (gi in 0:lev_res) {

      ####################################
      #Soft Thresholding
      ####################################
      temp <- pos.prob[(2^gi):(2^(gi + 1) - 1)]- alp*((1/2^(gi -1))^(1-1/7))
      pos.prob[(2^gi):(2^(gi + 1) - 1)] <-ifelse(temp<0,0, temp)
      ####################################
      #proportion of association per level
      ####################################
      p_vec[(gi+1)]    <-  mean(pos.prob[(2^gi):(2^(gi + 1) - 1)])
      lambcom[(gi+1)]  <-  mean(pos.prob[(2^gi):(2^(gi + 1) - 1)]*dnorm( betasub[(2^gi):(2^(gi + 1) - 1)] ,m1.hat,sigma1.hat)-(1-pos.prob[(2^gi):(2^(gi + 1) - 1)])*dnorm( betasub[(2^gi):(2^(gi + 1) - 1)] ,m0.hat,sigma0.hat))
    }
    porth   <- rep(NA,(lev_res))
    start <- 2^(1:lev_res)
    end  <-2^((1+1):(lev_res+1))-1
    for( gi in 0:(lev_res-1))
    {
      tempstart <- round(start + (gi)*(end-start)/lev_res)
      tempend <-round(start + (gi+1)*(end-start)/lev_res)
      ind <- c()
      for ( j in 1:lev_res)
      {
        temp <- cbind(tempstart,tempend)
        p1 <- temp[j,1]
        p2 <- temp[j,2]
        ind <- c(ind, p1:p2 )
      }
      porth[gi+1] <-  mean(pos.prob[ind])

    }
    ph <- sum(p_vec)
    pv <- sum(porth)
    min_ph_pv <- min( ph,pv)
    L_h <- sum(lambcom)
    return(c(L_h, min_ph_pv))
  }


  ################################################
  #Computing robustified variance for the diagonal
  ################################################
  if( length(dim(res) )==0) {
    N=1000
    #Generate random signal
    SNP=1:2^(lev_res+3)
    bp= SNP
    #Preparing for the wavelet transform
    Time01 <- (bp- min(bp))/(max(bp)-min(bp))
    loci<- matrix(runif(N*2^(lev_res+3),min=0,2),ncol=N)
    for ( i in 1:N)
    {
      loci[,i]  <- rnorm(n = 2^(lev_res+3))


    }

    bp=1:2^(lev_res)
    #Preparing for the wavelet transform
    Gen_W_trans <- apply(loci,2,my_wavproc)
    Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)
    #Compute proxy for empirical covariance matrix
    null_sd <-  mean(res[, grep(pattern = "null_sd_start_EM",names(res))])
    var_sim <-  null_sd^2
    emp_cov <- (cov(Gen_W_trans))
    emp_cov <-  var_sim * (emp_cov/(max(emp_cov)))

  }   else   {
    N=1000
    #Generate random signal
    SNP=1:2^(lev_res+3)
    bp= SNP
    #Preparing for the wavelet transform
    Time01 <- (bp- min(bp))/(max(bp)-min(bp))
    loci<- matrix(runif(N*2^(lev_res+3),min=0,2),ncol=N)
    for ( i in 1:N)
    {
      loci[,i]  <- rnorm(n = 2^(lev_res+3))


    }

    bp=1:2^(lev_res)
    #Preparing for the wavelet transform
    Gen_W_trans <- apply(loci,2,my_wavproc)
    Gen_W_trans = apply(Gen_W_trans, 1, Quantile_transform)
    emp_cov <- cov(Gen_W_trans)
    emp_cov <-  (mad(c(res [, grep( pattern = "Beta*", colnames(res))]))^2)*emp_cov/(max(emp_cov))


    var_sim <- (mad(c(res [, grep( pattern = "Beta*", colnames(res))]))^2)
    null_sd <-sqrt(var_sim)
  }



  ######################
  #Set up for simulation
  ######################

  #null_sd <-  mean(res[, grep(pattern = "null_sd_start_EM",colnames(res))])
  Pi_nt <- list()
  alt_sd <- 100*null_sd
  if(missing(base_shrink))
  {
    alp <-  1/sqrt(2*log(smp_size))
  }
  else
  {
    alp <-  base_shrink
  }

  print("Simulation of test statistics")
  temp <- seq(from=1,to=size,by=size/10)[-1]-1
  out <- list()
  my_f <- function(y)
  {
    pis <- max_EM_post_Beta(y, lev_res = lev_res , null_sd = null_sd ,alt_sd = alt_sd,alp=alp)
    return(pis)
  }
  if(print==TRUE)
  {

    for (i in 1:10)
    {
      y <- rmvnorm(floor(size/10),mean=rep(0,nbeta),sigma =emp_cov )
      out[[i]] <-t( apply(y,1,my_f))
      print(paste(i*floor(size/10), "simulations performed"))
    }

  }
  if(print==FALSE)
  {
    for (i in 1:10)
    {
      y <- rmvnorm(floor(size/10),mean=rep(0,nbeta),sigma =emp_cov )
      out[[i]] <-t( apply(y,1,my_f))
    }
  }



  out <- do.call(rbind,out)
  colnames(out) <- c("L_h","min_ph_pv")
  return(out)
}

