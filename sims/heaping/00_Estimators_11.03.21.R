#Program: 00_Estimators_01.18.20.R
#Developed by: BES
#Purpose: This program contains functions for MSM, PG, and DR estimation when data are not heaped
#Created: 01.18.20
#Last Updated: 11.03.21

library(numDeriv)
library(boot)
library(geex)

##########################################################################################################
#### MSM estimation functions  ####
##########################################################################################################

##wts known
geex_msm_fixed <- function(data, CM0.IV, CM1.IV, logRR.IV, RR.IV){

  m_estimate(
    estFUN = msm_estFUN_fixed, 
    data   = data, 
    roots = c(CM1.IV, CM0.IV, logRR.IV, RR.IV), 
    compute_roots = FALSE)
}

msm_estFUN_fixed <- function(data){
  A <- data$exp
  Y <- data$out
  W <- data$wt
 
  function(theta){
    cm1num <- (A*W*Y)
    cm1den <- (A*W)
    cm0num <-((1-A)*W*Y)
    cm0den <-((1-A)*W)
    p<-length(theta)
    
    c(cm1num-cm1den*theta[p-3], 
      cm0num-cm0den*theta[p-2],
      theta[p-3]-(exp(theta[p-1]) * theta[p-2]),
      exp(theta[p-1])-theta[p])
    
  }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

MSM_est_fixed <- function(data,exposure,outcome,weight,CM0.IV,CM1.IV,logRR.IV,RR.IV) {
  data$exp<-exposure
  data$out<-outcome
  data$wt<-weight
  geex_results_fixed <- geex_msm_fixed(data,CM0.IV,CM1.IV,logRR.IV,RR.IV)
  
  RR.MSM.fixed<-geex_results_fixed@estimates[length(geex_results_fixed@estimates)]
  seRR.MSM.fixed <- sqrt(geex_results_fixed@vcov[length(geex_results_fixed@estimates),length(geex_results_fixed@estimates)])
  MSM_ests_fixed<-cbind(RR.MSM.fixed,seRR.MSM.fixed)
  
  return(MSM_ests_fixed)
}


##wts estimated
geex_msm <- function(data, propensity_formula, CM0.IV, CM1.IV, logRR.IV, RR.IV){
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  models <- list(e = e_model)

  m_estimate(
    estFUN = msm_estFUN, 
    data   = data, 
    roots = c(coef(e_model), CM1.IV, CM0.IV, logRR.IV, RR.IV), 
    compute_roots = FALSE,
    outer_args = list(models = models))
}

msm_estFUN <- function(data,models){
  A <- data$exp
  Y <- data$out
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 1:ncol(Xe)
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    e  <- plogis(Xe %*% theta[e_pos])
    cm1num <- (A*Y/e)
    cm1den <- (A/e)
    cm0num <-(((1-A)*Y)/(1 - e))
    cm0den <-((1-A)/(1 - e))
    p<-length(theta)
    
    c(e_scores(theta[e_pos]), 
      cm1num-cm1den*theta[p-3], 
      cm0num-cm0den*theta[p-2],
      theta[p-3]-(exp(theta[p-1]) * theta[p-2]),
      exp(theta[p-1])-theta[p])
  }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

MSM_est <- function(data,exposure,outcome,propensity_formula,CM0.IV,CM1.IV,logRR.IV,RR.IV) {
data$exp<-exposure
data$out<-outcome
geex_results <- geex_msm(data,propensity_formula,CM0.IV,CM1.IV,logRR.IV,RR.IV)

RR.MSM<-geex_results@estimates[length(geex_results@estimates)]
seRR.MSM <- sqrt(geex_results@vcov[length(geex_results@estimates),length(geex_results@estimates)])
MSM_ests<-cbind(RR.MSM,seRR.MSM)

return(MSM_ests)
}


##########################################################################################################
#### Score Equations for the four distributions
##########################################################################################################

ScorePois<- function(A,L1,L2,L3,Y,lam){
  c(Y-lam,
    A*Y-A*lam,
    L1*Y-L1*lam,
    L2*Y-L2*lam, 
    L3*Y-L3*lam)
}

ScorePoisCig<- function(A,L1,Y,lam){
  c(Y-lam,
    A*Y-A*lam,
    L1*Y-L1*lam)
}

ScorePoisMisSpec<- function(A,L2,L3,Y,lam){
  c(Y-lam,
    A*Y-A*lam,
    L2*Y-L2*lam, 
    L3*Y-L3*lam)
}

calc_sum<-function(X,odisp){
  sum<-0 
  for (j in 0:(X-1)){
    sumpart<-(X>0)*(j+odisp^(-1))^-1
    sum<-sum+sumpart
  }
  return(sum)}

ScoreNBin<- function(A,L1,L2,L3,odisp,Y,lam){
  c(1*(Y-lam)/(1+odisp*lam),
    A*(Y-lam)/(1+odisp*lam),
    L1*(Y-lam)/(1+odisp*lam),
    L2*(Y-lam)/(1+odisp*lam),
    L3*(Y-lam)/(1+odisp*lam),
    odisp^(-2)*(log(1+odisp*lam)-calc_sum(Y,odisp))+(Y-lam)/(odisp*(1+odisp*lam)))
}

ScoreNBinMisSpec<- function(A,L2,L3,odisp,Y,lam){
  c(1*(Y-lam)/(1+odisp*lam),
    A*(Y-lam)/(1+odisp*lam),
    L2*(Y-lam)/(1+odisp*lam),
    L3*(Y-lam)/(1+odisp*lam),
    odisp^(-2)*(log(1+odisp*lam)-calc_sum(Y,odisp))+(Y-lam)/(odisp*(1+odisp*lam)))
}


ScoreZIP<- function(A,L1,L2,L3,Y,lam,nu){
c(as.numeric(Y==0)*(-1*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*1),
  as.numeric(Y==0)*(-A*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*A),
  as.numeric(Y==0)*(-L1*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*L1),
  as.numeric(Y==0)*(-L2*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*L2),
  as.numeric(Y==0)*(-L3*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*L3),
  as.numeric(Y==0)*(1*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*1,
  as.numeric(Y==0)*(L1*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*L1,
  as.numeric(Y==0)*(L2*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*L2,
  as.numeric(Y==0)*(L3*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*L3
)
}
    
ScoreZIPMisSpec<- function(A,L2,L3,Y,lam,nu){
  c(as.numeric(Y==0)*(-1*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*1),
    as.numeric(Y==0)*(-A*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*A),
    as.numeric(Y==0)*(-L2*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*L2),
    as.numeric(Y==0)*(-L3*lam/(nu*exp(lam)+1))+as.numeric(Y>0)*((Y-lam)*L3),
    as.numeric(Y==0)*(1*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*1,
    as.numeric(Y==0)*(L2*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*L2,
    as.numeric(Y==0)*(L3*nu*exp(lam)/(nu*exp(lam)+1))-(nu/(1+nu))*L3
  )}

calc_sum3<-function(X,odisp){
  sum<-0 
  for (j in 0:(X-1)){
    sumpart<-(X>0)*(-1/(j*odisp^(2)+odisp))
    sum<-sum+sumpart
  }
  return(sum)}

ScoreZINB<- function(A,L1,L2,L3,odisp,Y,lam,nu){
  c(as.numeric(Y==0)*1*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*1*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*A*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*A*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*L1*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*L1*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*L2*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*L2*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*L3*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*L2*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*1*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-1*nu/(1+nu),
    as.numeric(Y==0)*L1*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-L1*nu/(1+nu),
    as.numeric(Y==0)*L2*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-L2*nu/(1+nu),
    as.numeric(Y==0)*L3*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-L3*nu/(1+nu),
    as.numeric(Y==0)*(((1+odisp*lam)*log(1+odisp*lam)-odisp*lam)/(odisp^2*(1+odisp*lam)*(nu*(1+odisp*lam)^(odisp^(-1))+1)))+as.numeric(Y>0)*(log(1+odisp*lam)*odisp^(-2)+(Y-lam)/(odisp*(1+odisp*lam))+calc_sum3(Y,odisp))
  )
    }

ScoreZINBMisSpec<- function(A,L2,L3,odisp,Y,lam,nu){
  c(as.numeric(Y==0)*1*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*1*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*A*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*A*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*L2*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*L2*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*L3*((-lam*(1+odisp*lam)^(-1-odisp^(-1)))/(nu+(1+odisp*lam)^(-odisp^(-1))))+as.numeric(Y>0)*L2*((Y-lam)/(1+odisp*lam)),
    as.numeric(Y==0)*1*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-1*nu/(1+nu),
    as.numeric(Y==0)*L2*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-L2*nu/(1+nu),
    as.numeric(Y==0)*L3*(nu/(nu+(1+odisp*lam)^(-odisp^(-1))))-L3*nu/(1+nu),
    as.numeric(Y==0)*(((1+odisp*lam)*log(1+odisp*lam)-odisp*lam)/(odisp^2*(1+odisp*lam)*(nu*(1+odisp*lam)^(odisp^(-1))+1)))+as.numeric(Y>0)*(log(1+odisp*lam)*odisp^(-2)+(Y-lam)/(odisp*(1+odisp*lam))+calc_sum3(Y,odisp))
  )}

##########################################################################################################
#### Max LL for initial values, NB and ZI models
##########################################################################################################

nLLNBin<- function(par,data){
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  lam=exp(par[1]+par[2]*A+par[3]*L1+par[4]*L2+par[5]*L3)
  logodisp<-par[6]
  LL=sum(Y*log(exp(logodisp)*lam/(1+exp(logodisp)*lam))-(exp(logodisp))^(-1)*log(1+exp(logodisp)*lam)+lgamma(Y+(exp(logodisp))^(-1))-lgamma(Y+1)-lgamma((exp(logodisp))^(-1))) 
  return(-LL)
}

ML.NB<-function(data,init.NB){
  MLEs<-optim(init.NB, nLLNBin, data=data, method = "BFGS")
  MLEs2<-c(MLEs$par[1:5],exp(MLEs$par[6]))
  return(MLEs2)
}

nLLZIP<- function(par,data){
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  lam=exp(par[1]+par[2]*A+par[3]*L1+par[4]*L2+par[5]*L3)
  nu=exp(par[6]+par[7]*L1+par[8]*L2+par[9]*L3)
  LL=sum(as.numeric(Y==0)*log(nu+exp(-lam)))+sum(as.numeric(Y>0)*(Y*log(lam)-lam-lfactorial(Y)))-sum(log(1+nu)) 
  return(-LL)
}

ML.ZIP<-function(data){
  MLEs<-optim(init.ZIP, nLLZIP, data=simdat2, method="BFGS")
  return(MLEs$par)
}

calc_sum2<-function(X,odisp){
  sum<-0 
  for (j in 0:(X-1)){
    sumpart<-(X>0)*log(j+odisp^(-1))
    sum<-sum+sumpart
  }
  return(sum)}

nLLZINB<- function(par,data){
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  lam=exp(par[1]+par[2]*A+par[3]*L1+par[4]*L2+par[5]*L3)
  nu=exp(par[6]+par[7]*L1+par[8]*L2+par[9]*L3)
  logodisp<-par[10]
  LL<-0
  for (i in 1:num.part){
    partsum<-as.numeric(Y[i]==0)*log(nu[i]+(1+exp(logodisp)*lam[i])^(-exp(logodisp)^(-1)))+as.numeric(Y[i]>0)*(-lfactorial(Y[i])-(Y[i]+exp(logodisp)^(-1))*log(1+exp(logodisp)*lam[i])+Y[i]*log(exp(logodisp))+Y[i]*log(lam[i]))-log(1+nu[i])+as.numeric(Y[i]>0)*calc_sum2(Y[i],exp(logodisp))
    LL<-LL+partsum}
  return(-LL)
}

ML.ZINB<-function(data,init.ZINB){
  MLEs<-optim(init.ZINB, nLLZINB, data=simdat2, method="BFGS")
  MLEs2<-c(MLEs$par[1:9],exp(MLEs$par[10]))
  return(MLEs2)
}


nLLNBin.mo<- function(par,data){
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  lam=exp(par[1]+par[2]*A+par[3]*L2+par[4]*L3)
  logodisp<-par[5]
  LL=sum(Y*log(exp(logodisp)*lam/(1+exp(logodisp)*lam))-(exp(logodisp))^(-1)*log(1+exp(logodisp)*lam)+lgamma(Y+(exp(logodisp))^(-1))-lgamma(Y+1)-lgamma((exp(logodisp))^(-1))) 
  return(-LL)
}

ML.NB.mo<-function(data,init.NB){
  MLEs<-optim(init.NB, nLLNBin.mo, data=data, method = "BFGS")
  MLEs2<-c(MLEs$par[1:4],exp(MLEs$par[5]))
  return(MLEs2)
}

nLLZIPms<- function(par,data){
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  lam=exp(par[1]+par[2]*A+par[3]*L2+par[4]*L3)
  nu=exp(par[5]+par[6]*L2+par[7]*L3)
  LL=sum(as.numeric(Y==0)*log(nu+exp(-lam)))+sum(as.numeric(Y>0)*(Y*log(lam)-lam-lfactorial(Y)))-sum(log(1+nu)) 
  return(-LL)
}

ML.ZIP.ms<-function(data){
  MLEs<-optim(init.ZIP.ms, nLLZIPms, data=simdat2, method="BFGS")
  return(MLEs$par)
}

nLLZINBms<- function(par,data){
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  lam=exp(par[1]+par[2]*A+par[3]*L2+par[4]*L3)
  nu=exp(par[5]+par[6]*L2+par[7]*L3)
  logodisp<-par[8]
  LL<-0
  for (i in 1:num.part){
    partsum<-as.numeric(Y[i]==0)*log(nu[i]+(1+exp(logodisp)*lam[i])^(-exp(logodisp)^(-1)))+as.numeric(Y[i]>0)*(-lfactorial(Y[i])-(Y[i]+exp(logodisp)^(-1))*log(1+exp(logodisp)*lam[i])+Y[i]*log(exp(logodisp))+Y[i]*log(lam[i]))-log(1+nu[i])+as.numeric(Y[i]>0)*calc_sum2(Y[i],exp(logodisp))
    LL<-LL+partsum}
  return(-LL)
}

ML.ZINB.ms<-function(data){
  MLEs<-optim(init.ZINB.ms, nLLZINBms, data=simdat2, method="BFGS")
  MLEs2<-c(MLEs$par[1:7],exp(MLEs$par[8]))
  return(MLEs2)
}


##########################################################################################################
#### Parametric g-formula ####
##########################################################################################################

###Correct Specification
geex_PG <- function(data, estfun, MOD.IV, CM0.IV, CM1.IV, RR.IV){
 
  m_estimate(
    estFUN = estfun, 
    data   = data, 
    roots = c(MOD.IV,CM1.IV,CM0.IV,RR.IV),
    compute_roots = FALSE)
}

estfun_gf_Pois <- function(data){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    CM0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    CM1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    c(ScorePois(A,L1,L2,L3,Y,lam),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

estfun_gf_NBin <- function(data){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    CM0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    CM1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    odisp<-theta[6]
    c(ScoreNBin(A,L1,L2,L3,odisp,Y,lam),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

estfun_gf_ZIP <- function(data){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    nu=exp(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3)
    CM0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    CM1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    c(ScoreZIP(A,L1,L2,L3,Y,lam,nu),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

estfun_gf_ZINB <- function(data){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    nu=exp(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3)
    CM0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    CM1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    odisp<-theta[10]
    c(ScoreZINB(A,L1,L2,L3,odisp,Y,lam,nu),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}


estfun_gf_Pois_cig <- function(data){
  
  A<-data$inc
  L1<-data$income
  Y<-data$CIG.heap
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*L1+theta[3]*A)
    CM0<-exp(theta[1]+theta[2]*L1)
    CM1<-exp(theta[1]+theta[2]*L1+theta[3])
    c(ScorePoisCig(A,L1,Y,lam),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

###Incorrect Specification
estfun_gf_inc_Pois <- function(data){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    CM0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)
    CM1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)
    c(ScorePoisMisSpec(A,L2,L3,Y,lam),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

estfun_gf_inc_NBin <- function(data){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    CM0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)
    CM1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)
    odisp<-theta[5]
    c(ScoreNBinMisSpec(A,L2,L3,odisp,Y,lam),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

estfun_gf_inc_ZIP <- function(data){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    nu=exp(theta[5]+theta[6]*L2+theta[7]*L3)
    CM0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    CM1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    c(ScoreZIPMisSpec(A,L2,L3,Y,lam,nu),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}

estfun_gf_inc_ZINB <- function(data){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  function(theta){
    p<-length(theta)
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    nu=exp(theta[5]+theta[6]*L2+theta[7]*L3)
    CM0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    CM1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    odisp<-theta[8]
    c(ScoreZINBMisSpec(A,L2,L3,odisp,Y,lam,nu),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]-theta[p]*theta[p-1])
  }
}


PG_est <- function(data, estfun, MOD.IV, CM0.IV, CM1.IV, RR.IV) {
  geex_resultspg <- geex_PG(data, estfun, MOD.IV, CM0.IV, CM1.IV, RR.IV)
  
  RR.PG<-geex_resultspg@estimates[length(geex_resultspg@estimates)]
  seRR.PG <- sqrt(geex_resultspg@vcov[length(geex_resultspg@estimates),length(geex_resultspg@estimates)])
  PG_ests<-cbind(RR.PG,seRR.PG)
  
  return(PG_ests)
}


##########################################################################################################
#### DR Estimation ####
##########################################################################################################

#correct outcome model specification
geex_DR <- function(data, estfun, propensity_formula, MOD.IV, CM0.IV, CM1.IV, RR.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  models <- list(e = e_model)
 
  m_estimate(
    estFUN = estfun, 
    data   = data, 
    roots = c(MOD.IV,coef(e_model), CM1.IV, CM0.IV, RR.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
}



estfun_dr_Pois <- function(data,models){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 6:(5+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    m0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    m1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScorePois(A,L1,L2,L3,Y,lam),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_NBin <- function(data,models){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 7:(6+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    m0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    m1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    odisp<-theta[6]
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScoreNBin(A,L1,L2,L3,odisp,Y,lam),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_ZIP <- function(data,models){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 10:(9+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    nu=exp(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3)
    m0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    m1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScoreZIP(A,L1,L2,L3,Y,lam,nu),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_ZINB <- function(data,models){
  
  A<-data$inc
  L1<-data$druguse
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 11:(10+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L1+theta[4]*L2+theta[5]*L3)
    nu=exp(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3)
    m0<-exp(theta[1]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    m1<-exp(theta[1]+theta[2]+theta[3]*L1+theta[4]*L2+theta[5]*L3)*(1-inv.logit(theta[6]+theta[7]*L1+theta[8]*L2+theta[9]*L3))
    odisp<-theta[10]
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScoreZINB(A,L1,L2,L3,odisp,Y,lam,nu),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_Pois_Cig <- function(data,models){
  
  A<-data$inc
  L1<-data$income
  Y<-data$CIG.heap
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 4:(3+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*L1+theta[3]*A)
    m0<-exp(theta[1]+theta[2]*L1)
    m1<-exp(theta[1]+theta[3]+theta[2]*L1)
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScorePoisCig(A,L1,Y,lam),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

#Incorrect outcome model specification
estfun_dr_inc_Pois <- function(data,models){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 5:(4+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    m0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)
    m1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScorePoisMisSpec(A,L2,L3,Y,lam),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_inc_NBin <- function(data,models){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 6:(5+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    m0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)
    m1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)
    odisp<-theta[5]
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScoreNBinMisSpec(A,L2,L3,odisp,Y,lam),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_inc_ZIP <- function(data,models){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 8:(7+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    nu=exp(theta[5]+theta[6]*L2+theta[7]*L3)
    m0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    m1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScoreZIPMisSpec(A,L2,L3,Y,lam,nu),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

estfun_dr_inc_ZINB <- function(data,models){
  
  A<-data$inc
  L2<-data$sxex
  L3<-data$agevis10
  Y<-data$SP
  
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  e_pos  <- 9:(8+ncol(Xe))
  e_scores  <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    lam=exp(theta[1]+theta[2]*A+theta[3]*L2+theta[4]*L3)
    nu=exp(theta[5]+theta[6]*L2+theta[7]*L3)
    m0<-exp(theta[1]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    m1<-exp(theta[1]+theta[2]+theta[3]*L2+theta[4]*L3)*(1-inv.logit(theta[5]+theta[6]*L2+theta[7]*L3))
    odisp<-theta[8]
    CM1<-((A*Y - (A-e)*m1)/e)
    CM0<-(((1-A)*Y + (A-e)*m0)/(1-e)) 
    c(ScoreZINBMisSpec(A,L2,L3,odisp,Y,lam,nu),
      e_scores(theta[e_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p]*theta[p-1]- theta[p-2])
  }
}

DR_est <- function(data,estfun,propensity_formula, MOD.IV, CM0.IV, CM1.IV, RR.IV) {

  geex_resultsdr <- geex_DR(data,estfun,propensity_formula, MOD.IV, CM0.IV, CM1.IV, RR.IV)
  
  RR.DR<-geex_resultsdr@estimates[length(geex_resultsdr@estimates)]
  seRR.DR <- sqrt(geex_resultsdr@vcov[length(geex_resultsdr@estimates),length(geex_resultsdr@estimates)])
  DR_ests<-cbind(RR.DR,seRR.DR)
  
  return(DR_ests)
}







