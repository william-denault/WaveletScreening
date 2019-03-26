#'@title Internal EM procedure
#'@description  Internal EM procedure
#'@param my_betas a vector of beta of length 2^lev_res
#'@param lev_res the level of resolution in the wavelet transform
#'@param null_sd starting point for sd of the null distribution in the EM procedure
#'@param alt_sd starting point for sd of the alternative distribution in the EM procedure
#'@param alp shrinkage parameter for computation of the test statsitcs
#'@retrun The two test statistics used to build the final test (i.e L_h and min(ph,pv))



max_EM_post_Beta <- function(my_betas, lev_res,null_sd,alt_sd,alp) {
  niter = 1000
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
