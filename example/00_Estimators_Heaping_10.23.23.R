#Program: 00_Estimators_Heaping_06.28.23.R
#Developed by: BES
#Purpose: This program contains functions for MSM, PG, and DR estimation with data heaping

library(numDeriv)
library(boot)
library(geex)
library(rootSolve)

##########################################################################################################
#### LL and Score Equations
##########################################################################################################

#Poisson (for PG and DR estimators)

LLHeapPois<-function(Y,lam,logitpi1){
  sumpart<-log(ifelse(Y%%10>0,inv.logit(logitpi1)*dpois(Y,lam),inv.logit(logitpi1)*dpois(Y,lam)+(1-inv.logit(logitpi1))*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam))))
  neglogl<-(-1)*sum(sumpart)
  return(neglogl)
}

LLHeapPois.sumpart<-function(Y,lam,logitpi1){
  sumpart<-log(ifelse(Y%%10>0,inv.logit(logitpi1)*dpois(Y,lam),inv.logit(logitpi1)*dpois(Y,lam)+(1-inv.logit(logitpi1))*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam))))
  return(sumpart)
}

LLHeapPois.IH<-function(Y,heap.int,lam,pi1,pi2){
  sumpart<-as.numeric(heap.int==1)*LLHeapPois.sumpart(Y,lam,pi1)+as.numeric(heap.int==2)*LLHeapPois.sumpart(Y,lam,pi2)
  neglogl<-(-1)*sum(sumpart)
  return(neglogl)
}

#HCAR
ScoreHeapPois1cov<-function(Y,lam,COV,pi1){
  sumpart<-ifelse(Y%%10>0,((pi1*dpois(Y,lam))^(-1))*(pi1*COV*dpois(Y,lam)*(Y-lam)),((pi1*dpois(Y,lam)+(1-pi1)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(pi1*COV*dpois(Y,lam)*(Y-lam)+(1-pi1)*(COV*(Y-5-lam)*dpois(Y-5,lam)+COV*(Y-4-lam)*dpois(Y-4,lam)+COV*(Y-3-lam)*dpois(Y-3,lam)+COV*(Y-2-lam)*dpois(Y-2,lam)+COV*(Y-1-lam)*dpois(Y-1,lam)+COV*(Y-lam)*dpois(Y,lam)+COV*(Y+1-lam)*dpois(Y+1,lam)+COV*(Y+2-lam)*dpois(Y+2,lam)+COV*(Y+3-lam)*dpois(Y+3,lam)+COV*(Y+4-lam)*dpois(Y+4,lam))))
  return(sumpart)
}

ScoreHeapPoispi<-function(Y,lam,pi1){
  sumpart<-ifelse(Y%%10>0,((pi1*dpois(Y,lam))^(-1))*(dpois(Y,lam)),((pi1*dpois(Y,lam)+(1-pi1)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(dpois(Y,lam)-(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam))))
  return(sumpart)
}

# IH
ScoreHeapPois1cov.IH<-function(Y,lam,COV,pi1,pi2,int1){
  sumpart<-ifelse(Y%%10>0,int1*((pi1*dpois(Y,lam))^(-1))*(pi1*COV*dpois(Y,lam)*(Y-lam))+(1-int1)*((pi2*dpois(Y,lam))^(-1))*(pi2*COV*dpois(Y,lam)*(Y-lam)),int1*((pi1*dpois(Y,lam)+(1-pi1)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(pi1*COV*dpois(Y,lam)*(Y-lam)+(1-pi1)*(COV*(Y-5-lam)*dpois(Y-5,lam)+COV*(Y-4-lam)*dpois(Y-4,lam)+COV*(Y-3-lam)*dpois(Y-3,lam)+COV*(Y-2-lam)*dpois(Y-2,lam)+COV*(Y-1-lam)*dpois(Y-1,lam)+COV*(Y-lam)*dpois(Y,lam)+COV*(Y+1-lam)*dpois(Y+1,lam)+COV*(Y+2-lam)*dpois(Y+2,lam)+COV*(Y+3-lam)*dpois(Y+3,lam)+COV*(Y+4-lam)*dpois(Y+4,lam)))+(1-int1)*((pi2*dpois(Y,lam)+(1-pi2)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(pi2*COV*dpois(Y,lam)*(Y-lam)+(1-pi2)*(COV*(Y-5-lam)*dpois(Y-5,lam)+COV*(Y-4-lam)*dpois(Y-4,lam)+COV*(Y-3-lam)*dpois(Y-3,lam)+COV*(Y-2-lam)*dpois(Y-2,lam)+COV*(Y-1-lam)*dpois(Y-1,lam)+COV*(Y-lam)*dpois(Y,lam)+COV*(Y+1-lam)*dpois(Y+1,lam)+COV*(Y+2-lam)*dpois(Y+2,lam)+COV*(Y+3-lam)*dpois(Y+3,lam)+COV*(Y+4-lam)*dpois(Y+4,lam)))) 
  return(sumpart)
}

ScoreHeapPoispi.IH<-function(Y,lam,pi,int){
  sumpart<-ifelse(Y%%10>0,int*((pi*dpois(Y,lam))^(-1))*(dpois(Y,lam)),int*((pi*dpois(Y,lam)+(1-pi)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(dpois(Y,lam)-(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam))))  
  return(sumpart)
}

#Negative Binomial (for marginal models)

#HCAR
LLHeapNB<-function(Y,lam,pi,theta){
  sumpart<-log(ifelse(Y%%10>0,inv.logit((pi))*dnbinom(Y, mu=lam, size=1/theta),inv.logit(pi)*dnbinom(Y, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y-5, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y-4, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y-3, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y-2, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y-1, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y+1, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y+2, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y+3, mu=lam, size=1/theta)+(1-inv.logit(pi))*dnbinom(Y+4, mu=lam, size=1/theta)))
  neglogl<-(-1)*sum(sumpart)
  return(neglogl)
}

#deriv of NBin pdf with respect to beta
derivnbin<-function(odisp,Z,lam,COV){
  deriv<-(odisp^(Z)*gamma(Z+odisp^(-1)))/(factorial(Z)*gamma(odisp^(-1)))*(Z*COV*(lam^(Z))*(1+odisp*lam)^(-1*(odisp^(-1)+Z))-odisp*lam^(Z+1)*COV*(Z+odisp^(-1))*(1+odisp*lam)^(-1*(Z+odisp^(-1))-1))
  return(deriv)  
}

#pdf of NBin
Nbinpdf<-function(theta,Z){
  Y<-Z
  function(theta){
    (gamma(Y+theta[1]^(-1))/(factorial(Y)*gamma(theta[1]^(-1))))*((theta[1]*theta[2]/(1+theta[1]*theta[2]))^(Y))*(1/(1+theta[1]*theta[2]))^(theta[1]^(-1))
  }}

#approx deriv of NBin pdf wrt overdispersion parameter
derivnbinodisp<-function(theta,Z){
  grad(Nbinpdf(theta,Z),x=theta)[1]
}

#components of score eqn
HeapNBinScorepart<-function(Y,lam,odisp,pi){
  (1-pi)*dnbinom(Y-5, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y-4, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y-3, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y-2, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y-1, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+1, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+2, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+3, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+4, mu=lam, size=1/odisp)
}

HeapNBinScorepartpos<-function(Y,lam,odisp,pi){
  (1-pi)*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+1, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+2, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+3, mu=lam, size=1/odisp)+(1-pi)*dnbinom(Y+4, mu=lam, size=1/odisp)
}

#score for a covariate
ScoreHeapNB1cov<-function(Y,lam,COV,pi,odisp){
  sumpart<-ifelse(Y%%10>0,((pi*dnbinom(Y, mu=lam, size=1/odisp))^(-1))*(pi*derivnbin(odisp,Y,lam,COV)),ifelse(Y==0,((pi*dnbinom(Y, mu=lam, size=1/odisp)+HeapNBinScorepartpos(Y,lam,odisp,pi)))^(-1)*(pi*derivnbin(odisp,Y,lam,COV)+(1-pi)*derivnbin(odisp,Y,lam,COV)+(1-pi)*derivnbin(odisp,Y+1,lam,COV)+(1-pi)*derivnbin(odisp,Y+2,lam,COV)+(1-pi)*derivnbin(odisp,Y+3,lam,COV)+(1-pi)*derivnbin(odisp,Y+4,lam,COV)),((pi*dnbinom(Y, mu=lam, size=1/odisp)+ HeapNBinScorepart(Y,lam,odisp,pi)))^(-1)*(pi*derivnbin(odisp,Y,lam,COV)+(1-pi)*derivnbin(odisp,Y-5,lam,COV)+(1-pi)*derivnbin(odisp,Y-4,lam,COV)+(1-pi)*derivnbin(odisp,Y-3,lam,COV)+(1-pi)*derivnbin(odisp,Y-2,lam,COV)+(1-pi)*derivnbin(odisp,Y-1,lam,COV)+(1-pi)*derivnbin(odisp,Y,lam,COV)+(1-pi)*derivnbin(odisp,Y+1,lam,COV)+(1-pi)*derivnbin(odisp,Y+2,lam,COV)+(1-pi)*derivnbin(odisp,Y+3,lam,COV)+(1-pi)*derivnbin(odisp,Y+4,lam,COV))))
  return(sumpart)
}

#scores for heaping prob model 
ScoreHeapNBinpi<-function (Y,lam,pi,odisp){
  sumpart<-ifelse(Y%%10>0,((pi*dnbinom(Y, mu=lam, size=1/odisp))^(-1))*dnbinom(Y, mu=lam, size=1/odisp), ifelse(Y==0, (((pi*dnbinom(Y, mu=lam, size=1/odisp)+ HeapNBinScorepartpos(Y,lam,odisp,pi))^(-1))*(dnbinom(Y, mu=lam, size=1/odisp)-(dnbinom(Y, mu=lam, size=1/odisp)+dnbinom(Y+1, mu=lam, size=1/odisp)+dnbinom(Y+2, mu=lam, size=1/odisp)+dnbinom(Y+3, mu=lam, size=1/odisp)+dnbinom(Y+4, mu=lam, size=1/odisp)))), (((pi*dnbinom(Y, mu=lam, size=1/odisp)+ HeapNBinScorepart(Y,lam,odisp,pi))^(-1))*(dnbinom(Y, mu=lam, size=1/odisp)-(dnbinom(Y-5, mu=lam, size=1/odisp)+dnbinom(Y-4, mu=lam, size=1/odisp)+dnbinom(Y-3, mu=lam, size=1/odisp)+dnbinom(Y-2, mu=lam, size=1/odisp)+dnbinom(Y-1, mu=lam, size=1/odisp)+dnbinom(Y, mu=lam, size=1/odisp)+dnbinom(Y+1, mu=lam, size=1/odisp)+dnbinom(Y+2, mu=lam, size=1/odisp)+dnbinom(Y+3, mu=lam, size=1/odisp)+dnbinom(Y+4, mu=lam, size=1/odisp))))))
  return(sumpart)
}

#score for odisp parameter
ScoreHeapNBodisp<-function(Y,lam,pi,odisp){
  sumpart<-ifelse(Y%%10>0,
                  ((pi*dnbinom(Y, mu=lam, size=1/odisp))^(-1))*(pi*derivnbinodisp(theta=c(odisp,lam),Z=Y)),
                  ifelse(Y==0, ((pi*dnbinom(Y, mu=lam, size=1/odisp)+HeapNBinScorepartpos(Y,lam,odisp,pi)))^(-1)*(pi*derivnbinodisp(theta=c(odisp,lam),Z=Y)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+1)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+2)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+3)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+4)),
                         ((pi*dnbinom(Y, mu=lam, size=1/odisp)+HeapNBinScorepart(Y,lam,odisp,pi)))^(-1)*(pi*derivnbinodisp(theta=c(odisp,lam),Z=Y)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y-5)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y-4)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y-3)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y-2)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y-1)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+1)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+2)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+3)+(1-pi)*derivnbinodisp(theta=c(odisp,lam),Z=Y+4))
                  ))
  return(sumpart)
}


##########################################################################################################
#### Stacked Estimating Equations
##########################################################################################################
#NOTE: covariates (L#) can be added/removed from these functions, as appropriate, to align with the number of covariates. 

###for IPTW
ScoreHeapNB.HCAR<- function(Y,lam,pi,odisp){
  c(ScoreHeapNB1cov(Y,lam,1,pi,odisp),
    ScoreHeapNBinpi(Y,lam,pi,odisp),
    ScoreHeapNBodisp(Y,lam,pi,odisp))}


###for PG

#HCAR
ScorePGHeapPois.HCAR<- function(A,L1,Y,lam,pi){
  c(ScoreHeapPois1cov(Y,lam,1,pi),
    ScoreHeapPois1cov(Y,lam,L1,pi),
    ScoreHeapPois1cov(Y,lam,A,pi),
    ScoreHeapPoispi(Y,lam,pi))}

#IH
ScorePGHeapPois.IH<- function(A,L1,Y,int1,lam,pi1,pi2){
  c(ScoreHeapPois1cov.IH(Y,lam,1,pi1,pi2,int1),
    ScoreHeapPois1cov.IH(Y,lam,L1,pi1,pi2,int1),
    ScoreHeapPois1cov.IH(Y,lam,A,pi1,pi2,int1),
    ScoreHeapPoispi.IH(Y,lam,pi1,int1),
    ScoreHeapPoispi.IH(Y,lam,pi2,(1-int1)))}


###for DR
ScoreDRHeapPois.HCAR<- function(A,L1,Y,lam.IPTW,odisp,pi.IPTW,lam.PG,pi.PG){
  c(ScoreHeapNB1cov(Y,lam.IPTW,1,pi.IPTW,odisp),
    ScoreHeapNBinpi(Y,lam.IPTW,pi.IPTW,odisp),
    ScoreHeapNBodisp(Y,lam.IPTW,pi.IPTW,odisp),
    ScoreHeapPois1cov(Y,lam.PG,1,pi.PG),
    ScoreHeapPois1cov(Y,lam.PG,L1,pi.PG),
    ScoreHeapPois1cov(Y,lam.PG,A,pi.PG),
    ScoreHeapPoispi(Y,lam.PG,pi.PG))}


##########################################################################################################
#### Find MLEs, Heaping Models
##########################################################################################################
#NOTE: covariates (L#) can be added/removed from these functions, as appropriate, to align with the number of covariates in the outcome model. 

#marginal heaping models, to estimate pi (IPTW) - HCAR
maxLLNB.heap.HCAR<-function(par,data){
  Y<-data$CIG.heap
  pi=par[1]
  theta=par[2]
  lam=exp(par[3])
  negLL=(LLHeapNB(Y,lam,pi,theta))
  return(negLL)
}

ML.Heap.marginal.HCAR<-function(data,init.NB){
  MLEs<-optim(init.NB, maxLLNB.heap.HCAR, data=data, method = "L-BFGS-B", lower = c(-Inf, 0.0001, -Inf), upper = c(Inf, Inf, Inf))
  MLEs2<-c(inv.logit(MLEs$par[1]),MLEs$par[2:3])
  return(MLEs2)
}

#PG models 
PGmaxLL.HCAR<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  L1<-data$income
  pi=par[1]
  lam=exp(par[2]+par[3]*L1+par[4]*A)
  negLL=(LLHeapPois(Y,lam,pi))
  return(negLL)
}

ML.HeapPoisPG.HCAR<-function(data,init.Pois){
  MLEs<-optim(init.Pois, PGmaxLL.HCAR, data=data, method="BFGS")
  MLEs2<-c(inv.logit(MLEs$par[1]),(MLEs$par[2:4]))
  return(MLEs2)
} 

PGmaxLL.IH<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  L1<-data$income
  heap.int<-data$interval.col
  pi1=par[1]
  pi2=par[2]
  lam=exp(par[3]+par[4]*L1+par[5]*A)
  negLL=(LLHeapPois.IH(Y,heap.int,lam,pi1,pi2))  
  return(negLL)
}

ML.HeapPoisPG.IH<-function(data,init.Pois){
  MLEs<-optim(init.Pois, PGmaxLL.IH, data=data, method="BFGS")
  MLEs2<-c(inv.logit(MLEs$par[1]),inv.logit(MLEs$par[2]),MLEs$par[3:5])
  return(MLEs2)
} 


##########################################################################################################
#### IPTW estimation functions  ####
##########################################################################################################

###HCAR

##wts known
geex_IPTW_fixed_heap_NB.HCAR <- function(data, IV){
  m_estimate(
    estFUN = IPTW_estFUN_fixed_heap_NB.HCAR, 
    data   = data, 
       roots = c(IV), 
        compute_roots = FALSE)
}

IPTW_estFUN_fixed_heap_NB.HCAR <- function(data){
  A <- data$exp
  Y <- data$out
  Y_Hall<-round(data$out+0.00001,-1)
  W <- data$wt
  #num and denom of estimators of E(Y_h^a)
  Yhcm1num <- (W*A*Y)
  Yhcm1den <- (W*A)
  Yhcm0num <-(W*(1-A)*Y)
  Yhcm0den <-(W*(1-A))
  #nums of estimators of E(h(Y^a)) (Note: denominators are the same as for E(Y_h^a)))
  Yhallcm1num <- (W*A*Y_Hall)
  Yhallcm0num <-(W*(1-A)*Y_Hall)

  function(theta){
    p<-length(theta)
    pi=theta[1]
    odisp<-theta[2]
    lam=exp(theta[3])
    c(ScoreHeapNB.HCAR(Y,lam,pi,odisp),
      Yhcm1num-Yhcm1den*theta[4],   
      Yhcm0num-Yhcm0den*theta[5],
      Yhallcm1num-Yhcm1den*theta[6], 
      Yhallcm0num-Yhcm0den*theta[7],
      ((theta[1]^(-1))*(theta[4]-(1-theta[1])*theta[6]))-theta[8],
      ((theta[1]^(-1))*(theta[5]-(1-theta[1])*theta[7]))-theta[9],
      theta[8]-theta[10]*theta[9])}}
  

IPTW_est_fixed_heap_NB.HCAR <- function(data,exposure,outcome,weight,IV) {
    data$exp<-exposure
    data$out<-outcome
    data$wt<-weight
    geex_results_fixed <- geex_IPTW_fixed_heap_NB.HCAR(data,IV)
    RR.IPTW.fixed<-geex_results_fixed@estimates[length(geex_results_fixed@estimates)]
    seRR.IPTW.fixed <- sqrt(geex_results_fixed@vcov[length(geex_results_fixed@estimates),length(geex_results_fixed@estimates)])
    IPTW_ests_fixed<-cbind(RR.IPTW.fixed,seRR.IPTW.fixed)
    return(IPTW_ests_fixed)
  }
  

  ##wts estimated
   geex_IPTW_heap.HCAR <- function(data, propensity_formula, IV){
     e_model  <- glm(propensity_formula, data = data, family =binomial)
     models <- list(e = e_model)
    
     m_estimate(
       estFUN = IPTW_estFUN_heap_NB.HCAR, 
       data   = data, 
       roots = c(coef(e_model), IV), 
       compute_roots = FALSE,
       outer_args = list(models = models))
   }
   
   IPTW_estFUN_heap_NB.HCAR <- function(data,models){
     A <- data$exp
     Y <- data$out
     Y_Hall<-round(data$out+0.00001,-1)
     Xe <- grab_design_matrix(data = data,
                              rhs_formula = grab_fixed_formula(models$e))
     e_pos  <- 1:ncol(Xe)
     e_scores  <- grab_psiFUN(models$e, data)
     lam_start<-ncol(Xe)+1
     
     function(theta){
       e <- plogis(Xe %*% theta[e_pos])
       W <-A*e^(-1)+(1-A)*(1-e)^(-1)
       pi=theta[lam_start]
       odisp=theta[lam_start+1]
       lam=exp(theta[lam_start+2])
       p<-length(theta)
       
       c(e_scores(theta[e_pos]), 
         ScoreHeapNB.HCAR(Y,lam,pi,odisp),
         W*A*Y-(W*A)*theta[lam_start+3],   
         (W*(1-A)*Y)-(W*(1-A))*theta[lam_start+4],
         (W*A*Y_Hall)-(W*A)*theta[lam_start+5], 
         (W*(1-A)*Y_Hall)-(W*(1-A))*theta[lam_start+6],
         ((theta[lam_start]^(-1))*(theta[lam_start+3]-(1-theta[lam_start])*theta[lam_start+5]))-theta[lam_start+7],
         ((theta[lam_start]^(-1))*(theta[lam_start+4]-(1-theta[lam_start])*theta[lam_start+6]))-theta[lam_start+8],
         theta[lam_start+7]-theta[lam_start+9]*theta[lam_start+8])
     }}
  
   IPTW_est_heap_NB.HCAR <- function(data,exposure,outcome,propensity_formula,IV) {
     data$exp<-exposure
     data$out<-outcome
     geex_results <- geex_IPTW_heap.HCAR(data,propensity_formula,IV)
     
     RR.IPTW<-geex_results@estimates[length(geex_results@estimates)]
     seRR.IPTW <- sqrt(geex_results@vcov[length(geex_results@estimates),length(geex_results@estimates)])
     IPTW_ests<-cbind(RR.IPTW,seRR.IPTW)
     
     return(IPTW_ests)
   }
   
 
   
   ##wts estimated
   geex_IPTW_heap.HCAR <- function(data, propensity_formula, IV){
     e_model  <- glm(propensity_formula, data = data, family =binomial)
     models <- list(e = e_model)
     
     m_estimate(
       estFUN = IPTW_estFUN_heap_NB.HCAR, 
       data   = data, 
       roots = c(coef(e_model), IV), 
       compute_roots = FALSE,
       outer_args = list(models = models))
   }
   
   IPTW_estFUN_heap_NB.HCAR <- function(data,models){
     A <- data$exp
     Y <- data$out
     Y_Hall<-round(data$out+0.00001,-1)
     Xe <- grab_design_matrix(data = data,
                              rhs_formula = grab_fixed_formula(models$e))
     e_pos  <- 1:ncol(Xe)
     e_scores  <- grab_psiFUN(models$e, data)
     lam_start<-ncol(Xe)+1
     
     function(theta){
       e <- plogis(Xe %*% theta[e_pos])
       W <-A*e^(-1)+(1-A)*(1-e)^(-1)
       pi=theta[lam_start]
       odisp=theta[lam_start+1]
       lam=exp(theta[lam_start+2])
       p<-length(theta)
       
       c(e_scores(theta[e_pos]), 
         ScoreHeapNB.HCAR(Y,lam,pi,odisp),
         W*A*Y-(W*A)*theta[lam_start+3],   
         (W*(1-A)*Y)-(W*(1-A))*theta[lam_start+4],
         (W*A*Y_Hall)-(W*A)*theta[lam_start+5], 
         (W*(1-A)*Y_Hall)-(W*(1-A))*theta[lam_start+6],
         ((theta[lam_start]^(-1))*(theta[lam_start+3]-(1-theta[lam_start])*theta[lam_start+5]))-theta[lam_start+7],
         ((theta[lam_start]^(-1))*(theta[lam_start+4]-(1-theta[lam_start])*theta[lam_start+6]))-theta[lam_start+8],
         theta[lam_start+7]-theta[lam_start+9]*theta[lam_start+8])
     }}
   
   IPTW_est_heap_NB.HCAR <- function(data,exposure,outcome,propensity_formula,IV) {
     data$exp<-exposure
     data$out<-outcome
     geex_results <- geex_IPTW_heap.HCAR(data,propensity_formula,IV)
     
     RR.IPTW<-geex_results@estimates[length(geex_results@estimates)]
     seRR.IPTW <- sqrt(geex_results@vcov[length(geex_results@estimates),length(geex_results@estimates)])
     IPTW_ests<-cbind(RR.IPTW,seRR.IPTW)
     
     return(IPTW_ests)
   }
  
  
  
  ##########################################################################################################
  #### Parametric g-formula ####
  ##########################################################################################################
  
   #NOTE: covariates (L#) can be added/removed from these functions, as appropriate, to align with the number of covariates in the outcome model. 
   
   #### HCAR ####
   
   
   geex_PG_heap.HCAR <- function(data, estfun, IV){
     
     m_estimate(
       estFUN = estfun, 
       data   = data, 
       roots = c(IV),
       compute_roots = FALSE)
   }
   
   estfun_gf_Pois_heap.HCAR <- function(data){
     
     A<-data$inc
     L1<-data$income
     Y<-data$CIG.heap
     
     
     function(theta){
       p<-length(theta)
       pi=theta[1]
       lam=exp(theta[2]+theta[3]*L1+theta[4]*A)
       CM0<-exp(theta[2]+theta[3]*L1)
       CM1<-exp(theta[2]+theta[3]*L1+theta[4])
       c(ScorePGHeapPois.HCAR(A,L1,Y,lam,pi),
         CM1-theta[p-2], 
         CM0-theta[p-1],
         theta[p-2]-theta[p]*theta[p-1])
     }
   }
   
   
   PG_est_heap.HCAR <- function(data, estfun, IV) {
     geex_resultspg <- geex_PG_heap.HCAR(data, estfun, IV)
     
     RR.PG<-geex_resultspg@estimates[length(geex_resultspg@estimates)]
     seRR.PG <- sqrt(geex_resultspg@vcov[length(geex_resultspg@estimates),length(geex_resultspg@estimates)])
     PG_ests<-cbind(RR.PG,seRR.PG)
     
     return(PG_ests)
   }
   
   
   
   #### IH ####
   
   
   geex_PG_heap.IH <- function(data, estfun, IV){
     
     m_estimate(
       estFUN = estfun, 
       data   = data, 
       roots = c(IV),
       compute_roots = FALSE)
   }
   
   estfun_gf_Pois_heap.IH <- function(data){
     
     A<-data$inc
     L1<-data$income
     Y<-data$CIG.heap
     int1<-as.numeric(data$interval.col==1)
     
     
     function(theta){
       p<-length(theta)
       pi1=theta[1]
       pi2=theta[2]
       lam=exp(theta[3]+theta[4]*L1+theta[5]*A)
       CM0<-exp(theta[3]+theta[4]*L1)
       CM1<-exp(theta[3]+theta[4]*L1+theta[5])
       c(ScorePGHeapPois.IH(A,L1,Y,int1,lam,pi1,pi2), 
         CM1-theta[p-2], 
         CM0-theta[p-1],
         theta[p-2]-theta[p]*theta[p-1])
     }
   }
   
   
   PG_est_heap.IH <- function(data, estfun, IV) {
     geex_resultspg <- geex_PG_heap.IH(data, estfun, IV)
     
     RR.PG<-geex_resultspg@estimates[length(geex_resultspg@estimates)]
     seRR.PG <- sqrt(geex_resultspg@vcov[length(geex_resultspg@estimates),length(geex_resultspg@estimates)])
     PG_ests<-cbind(RR.PG,seRR.PG)
     
     return(PG_ests)
   }
   
  
  
  
  ##########################################################################################################
  #### DR Estimation ####
  ##########################################################################################################
   #NOTE: covariates (L#) can be added/removed from these functions, as appropriate, to align with the number of covariates in the outcome model. 

   geex_DR_heap.HCAR <- function(data, estfun, propensity_formula, IV){
     
     e_model  <- glm(propensity_formula, data = data, family =binomial)
     models <- list(e = e_model)

     m_estimate(
       estFUN = estfun, 
       data   = data, 
       roots = c(coef(e_model), IV),
       compute_roots = FALSE,
       outer_args = list(models = models))
   }
   
   estfun_dr_heap.HCAR <- function(data,models){
     
     A<-data$inc
     L1<-data$income
     Y<-data$CIG.heap
     Y_Hall<-round(Y+0.00001,-1)
     
     Xe <- grab_design_matrix(data = data,
                              rhs_formula = grab_fixed_formula(models$e))
     e_pos  <- 1:(0+ncol(Xe))
     e_scores  <- grab_psiFUN(models$e, data)
     
     function(theta){
       p<-length(theta)
       #IPTW parameters
       e <- plogis(Xe %*% theta[e_pos]) 
       pi.IPTW=theta[ncol(Xe)+1]
       odisp=theta[ncol(Xe)+2]
       lam.IPTW=exp(theta[ncol(Xe)+3])
       #PG parameters
       pi.PG=theta[ncol(Xe)+4]
       lam.PG=exp(theta[ncol(Xe)+5]+theta[ncol(Xe)+6]*L1+theta[ncol(Xe)+7]*A)
       #components of DR estimator
       y0hat<-exp(theta[ncol(Xe)+5]+theta[ncol(Xe)+6]*L1)
       y1hat<-exp(theta[ncol(Xe)+5]+theta[ncol(Xe)+6]*L1+theta[ncol(Xe)+7])
       hy0hat<-round(y0hat+0.00001,-1)
       hy1hat<-round(y1hat+0.00001,-1)
       hatYh0<-pi.PG*y0hat+(1-pi.PG)*hy0hat
       hatYh1<-pi.PG*y1hat+(1-pi.PG)*hy1hat
              
       c(ScoreDRHeapPois.HCAR(A,L1,Y,lam.IPTW,odisp,pi.IPTW,lam.PG,pi.PG),  
         e_scores(theta[e_pos]),
         ((e^(-1))*(A*Y-(A-e)*hatYh1)) - theta[p-6], #E(Y_h^1)
         (((1-e)^(-1))*((1-A)*Y+(A-e)*hatYh0)) - theta[p-5], #E(Y_h^0)
         ((e^(-1))*(A*Y_Hall-(A-e)*hy1hat)) - theta[p-4], #E(h(Y^1))
         (((1-e)^(-1))*((1-A)*Y_Hall+(A-e)*hy0hat)) - theta[p-3], #E(h(Y^0))
         ((pi.IPTW^(-1))*(theta[p-6]-(1-pi.IPTW)*theta[p-4])) - theta[p-2], #CM1
         ((pi.IPTW^(-1))*(theta[p-5]-(1-pi.IPTW)*theta[p-3])) - theta[p-1], #CM0
         theta[p]*theta[p-1]- theta[p-2])}
   }
   
  
DR_est_heap.HCAR <- function(data,estfun,propensity_formula, IV) {
     
     geex_resultsdr <- geex_DR_heap.HCAR(data,estfun,propensity_formula, IV)
     
     RR.DR<-geex_resultsdr@estimates[length(geex_resultsdr@estimates)]
     seRR.DR <- sqrt(geex_resultsdr@vcov[length(geex_resultsdr@estimates),length(geex_resultsdr@estimates)])

     DR_ests<-cbind(RR.DR,seRR.DR)
     
    return(DR_ests)
   }
   
  
  
  
  
  
  