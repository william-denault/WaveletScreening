#'@title Simulation of the null statsitics
#'@description  Simulation of the null statsitics
#'@param emp_cov Emprical covariance matrix of the beta values. Can be computed using several results of the wavelet screaming using the betas values for different loci. If missing the function compute an approximation of the covariance matrix, this leads to a lost of power and a more conservative test statitics.
#'@param smp_size Sample size from the main run of Wavelet Screaming
#'@param lev_res the level of resolution in the wavelet transform
#'@param size number of simulation to be performed
#'@param sigma_b the parameter of the NIG prior used for the Betas computation.
#'@param print logical parameter set as TRUE, if TRUE send message when 10\% of the simulations has been completed.
#'@return The simulation under the null of the two test statistics used to build the final test (i.e L_h and min(ph,pv))

Simu_null_emp <- function(emp_cov,smp_size,lev_res,size,sigma_b,print=TRUE)
{

  #####################################
  #Size of the multivaraite to simulate
  #####################################
  nbeta <- sum(2^(0:lev_res))
  ######################################
  #compute theoretical variance of Betas
  ######################################

  if(!(dim(emp_cov)[1]==nbeta &dim(emp_cov)[2]==nbeta )){
    print("Warning: Empirical Covariance Matrix not of the good size")
    break
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
    sigma0.hat<-sqrt(null_sd)
    sigma1.hat<-alt_sd
    #Prevent from label swapping
    if(sigma1.hat < sigma0.hat){
      sigma1.hat <- 3*sigma0.hat+sigma1.hat
    }

    p.hat<-0.5
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
      m0.hat<-sum((1-temp)* betasub)/(sum(1-temp)+eps)
      sigma1.hat<-sqrt( sum(temp*( betasub-m1.hat)^2)/(sum(temp)+eps) )+alt_sd
      sigma0.hat<-sqrt( sum((1-temp)*( betasub-m0.hat)^2)/(sum(1-temp)+eps) )
      #limit the decrease of sigma0.hat in case of non identifiable mixture
      if(sigma0.hat < 0.1*sqrt(null_sd) ){
        sigma0.hat <- 0.1*sqrt(null_sd)
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
      temp <- pos.prob[(2^gi):(2^(gi + 1) - 1)]- alp*sqrt(1/2^(gi -1))
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
        ind <- c(p1:p2 )
      }
      porth[gi+1] <-  mean(pos.prob[ind])

    }
    ph <- sum(p_vec)
    pv <- sum(porth)
    min_ph_pv <- min( ph,pv)
    L_h <- sum(lambcom)
    return(c(L_h, min_ph_pv))
  }
  ######################
  #Set up for simulation
  ######################
  null_sd <- mean(diag(emp_cov))
  Pi_nt <- list()
  alt_sd <- sigma_b
  alp <- 1/sqrt(2*log(smp_size))
  print("Simulation of test statistics")
  temp <- seq(from=1,to=size,by=size/10)[-1]-1
  out <- list()
  my_f <- function(y)
  {
    pis <- max_EM_post_Beta(y, lev_res = 9, null_sd = null_sd ,alt_sd = alt_sd,alp=alp)
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
