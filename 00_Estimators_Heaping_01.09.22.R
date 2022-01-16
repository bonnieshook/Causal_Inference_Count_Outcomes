#Program: 00_Estimators_Heaping_04.08.21.R
#Developed by: BES
#Purpose: This program contains functions for MSM, PG, and DR estimation with data heaping

library(numDeriv)
library(boot)
library(geex)
library(rootSolve)

##########################################################################################################
#### LL and Score Equationds
##########################################################################################################

#Poisson (for PG and DR estimators)

LLHeapPois<-function(Y,lam,logitpi1){
  sumpart<-log(ifelse(Y%%10>0,inv.logit(logitpi1)*dpois(Y,lam),inv.logit(logitpi1)*dpois(Y,lam)+(1-inv.logit(logitpi1))*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam))))
  neglogl<-(-1)*sum(sumpart)
  return(neglogl)
}

ScoreHeapPois1cov<-function(Y,lam,COV,pi1){
  sumpart<-ifelse(Y%%10>0,((pi1*dpois(Y,lam))^(-1))*(pi1*COV*dpois(Y,lam)*(Y-lam)),((pi1*dpois(Y,lam)+(1-pi1)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(pi1*COV*dpois(Y,lam)*(Y-lam)+(1-pi1)*(COV*(Y-5-lam)*dpois(Y-5,lam)+COV*(Y-4-lam)*dpois(Y-4,lam)+COV*(Y-3-lam)*dpois(Y-3,lam)+COV*(Y-2-lam)*dpois(Y-2,lam)+COV*(Y-1-lam)*dpois(Y-1,lam)+COV*(Y-lam)*dpois(Y,lam)+COV*(Y+1-lam)*dpois(Y+1,lam)+COV*(Y+2-lam)*dpois(Y+2,lam)+COV*(Y+3-lam)*dpois(Y+3,lam)+COV*(Y+4-lam)*dpois(Y+4,lam))))
  return(sumpart)
}

ScoreHeapPois1cov.pi1<-function(Y,lam,pi1){
  sumpart<-ifelse(Y%%10>0,((pi1*dpois(Y,lam))^(-1))*(dpois(Y,lam)),((pi1*dpois(Y,lam)+(1-pi1)*(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam)))^(-1))*(dpois(Y,lam)-(dpois(Y-5,lam)+dpois(Y-4,lam)+dpois(Y-3,lam)+dpois(Y-2,lam)+dpois(Y-1,lam)+dpois(Y,lam)+dpois(Y+1,lam)+dpois(Y+2,lam)+dpois(Y+3,lam)+dpois(Y+4,lam))))
  return(sumpart)
}

#Negative Binomial
sumpdfheap<-function(Y,lam,odisp){
  dnbinom(Y-5, mu=lam, size=1/odisp)+dnbinom(Y-4, mu=lam, size=1/odisp)+dnbinom(Y-3, mu=lam, size=1/odisp)+dnbinom(Y-2, mu=lam, size=1/odisp)+dnbinom(Y-1, mu=lam, size=1/odisp)+dnbinom(Y, mu=lam, size=1/odisp)+dnbinom(Y+1, mu=lam, size=1/odisp)+dnbinom(Y+2, mu=lam, size=1/odisp)+dnbinom(Y+3, mu=lam, size=1/odisp)+dnbinom(Y+4, mu=lam, size=1/odisp)  
}

LLHeapNBWtd<-function(Y,lam,logitpi1,theta,W){
  sumpart<-W*log(ifelse(Y%%10>0,inv.logit(logitpi1)*dnbinom(Y, mu=lam, size=1/theta), inv.logit(logitpi1)*dnbinom(Y, mu=lam, size=1/theta)+(1-inv.logit(logitpi1))*(sumpdfheap(Y,lam,theta))))
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

ScoreHeapNB1covWtd<-function(Y,W,lam,COV,pi1,odisp){
  sumpart<-ifelse(Y%%10>0,((pi1*dnbinom(Y, mu=lam, size=1/odisp))^(-1))*(pi1*derivnbin(odisp,Y,lam,COV)),ifelse(Y==0, ((pi1*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi1)*(dnbinom(Y, mu=lam, size=1/odisp)+dnbinom(Y+1, mu=lam, size=1/odisp)+dnbinom(Y+2, mu=lam, size=1/odisp)+dnbinom(Y+3, mu=lam, size=1/odisp)+dnbinom(Y+4, mu=lam, size=1/odisp)))^(-1))*(pi1*derivnbin(odisp,Y,lam,COV)+(1-pi1)*(derivnbin(odisp,Y,lam,COV)+derivnbin(odisp,Y+1,lam,COV)+derivnbin(odisp,Y+2,lam,COV)+derivnbin(odisp,Y+3,lam,COV)+derivnbin(odisp,Y+4,lam,COV))),((pi1*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi1)*(sumpdfheap(Y,lam,odisp)))^(-1))*(pi1*derivnbin(odisp,Y,lam,COV)+(1-pi1)*(derivnbin(odisp,Y-5,lam,COV)+derivnbin(odisp,Y-4,lam,COV)+derivnbin(odisp,Y-3,lam,COV)+derivnbin(odisp,Y-2,lam,COV)+derivnbin(odisp,Y-1,lam,COV)+derivnbin(odisp,Y,lam,COV)+derivnbin(odisp,Y+1,lam,COV)+derivnbin(odisp,Y+2,lam,COV)+derivnbin(odisp,Y+3,lam,COV)+derivnbin(odisp,Y+4,lam,COV)))))
  return(W*sumpart)
}

ScoreHeapNB1cov.pi1Wtd<-function(Y,W,lam,pi1,odisp){
  sumpart<-ifelse(Y%%10>0,((pi1*dnbinom(Y, mu=lam, size=1/odisp))^(-1))*(dnbinom(Y, mu=lam, size=1/odisp)),((pi1*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi1)*(sumpdfheap(Y,lam,odisp)))^(-1))*(dnbinom(Y, mu=lam, size=1/odisp)-(sumpdfheap(Y,lam,odisp))))
  return(W*sumpart)
}

ScoreHeapNBodisp<-function(Y,W,lam,pi1,odisp){
  sumpart<-ifelse(Y%%10>0,((pi1*dnbinom(Y, mu=lam, size=1/odisp))^(-1))*(pi1*derivnbinodisp(theta=c(odisp,lam),Z=Y)),ifelse(Y==0, ((pi1*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi1)*(dnbinom(Y, mu=lam, size=1/odisp)+dnbinom(Y+1, mu=lam, size=1/odisp)+dnbinom(Y+2, mu=lam, size=1/odisp)+dnbinom(Y+3, mu=lam, size=1/odisp)+dnbinom(Y+4, mu=lam, size=1/odisp)))^(-1))*(pi1*derivnbinodisp(theta=c(odisp,lam),Z=Y)+(1-pi1)*(derivnbinodisp(theta=c(odisp,lam),Z=Y)+derivnbinodisp(theta=c(odisp,lam),Z=Y+1)+derivnbinodisp(theta=c(odisp,lam),Z=Y+2)+derivnbinodisp(theta=c(odisp,lam),Z=Y+3)+derivnbinodisp(theta=c(odisp,lam),Z=Y+4))),((pi1*dnbinom(Y, mu=lam, size=1/odisp)+(1-pi1)*(sumpdfheap(Y,lam,odisp)))^(-1))*(pi1*derivnbinodisp(theta=c(odisp,lam),Z=Y)+(1-pi1)*(derivnbinodisp(theta=c(odisp,lam),Z=Y-5)+derivnbinodisp(theta=c(odisp,lam),Z=Y-4)+derivnbinodisp(theta=c(odisp,lam),Z=Y-3)+derivnbinodisp(theta=c(odisp,lam),Z=Y-2)+derivnbinodisp(theta=c(odisp,lam),Z=Y-1)+derivnbinodisp(theta=c(odisp,lam),Z=Y)+derivnbinodisp(theta=c(odisp,lam),Z=Y+1)+derivnbinodisp(theta=c(odisp,lam),Z=Y+2)+derivnbinodisp(theta=c(odisp,lam),Z=Y+3)+derivnbinodisp(theta=c(odisp,lam),Z=Y+4)))))
  return(W*sumpart)
}


##########################################################################################################
#### Stacked Estimating Equations
##########################################################################################################

###for MSM 
ScoreMSMHeapNB<- function(A,Y,W,lam,pi1,odisp){
  c(ScoreHeapNB1covWtd(Y,W,lam,1,pi1,odisp),
    ScoreHeapNB1covWtd(Y,W,lam,A,pi1,odisp),
    ScoreHeapNB1cov.pi1Wtd(Y,W,lam,pi1,odisp),
    ScoreHeapNBodisp(Y,W,lam,pi1,odisp))}

###for PG
ScorePGHeapPois<- function(A,L1,Y,lam,pi1){
  c(ScoreHeapPois1cov(Y,lam,1,pi1),
    ScoreHeapPois1cov(Y,lam,L1,pi1),
    ScoreHeapPois1cov(Y,lam,A,pi1),
    ScoreHeapPois1cov.pi1(Y,lam,pi1))}


###for DR
ScoreDRHeapPois<- function(A,L1,W,Y,lam,pi1){
  c(ScoreHeapPois1cov(Y,lam,1,pi1),
    ScoreHeapPois1cov(Y,lam,A,pi1),
    ScoreHeapPois1cov(Y,lam,L1,pi1),
    ScoreHeapPois1cov(Y,lam,W,pi1),
    ScoreHeapPois1cov.pi1(Y,lam,pi1))}


##########################################################################################################
#### Find initial values, Heaping Models
##########################################################################################################

#MSM models 
MSMmaxLLNB<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  W<-data$IPTW 
  logitpi1=par[1]
  theta=par[2]
  lam=exp(par[3]+par[4]*A)
  negLL=(LLHeapNBWtd(Y,lam,logitpi1,theta,W))
  return(negLL)
}

ML.HeapNBMSM<-function(data,init.NB){
  MLEs<-optim(init.NB, MSMmaxLLNB, data=data, method = "L-BFGS-B", lower = c(-Inf, 0.0001, -Inf, -Inf), upper = c(Inf, Inf, Inf, Inf))
  MLEs2<-c(inv.logit(MLEs$par[1]),MLEs$par[2:4])
    return(MLEs2)
}


#PG
PGmaxLL<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  L1<-data$income

  logitpi1=par[1]
  lam=exp(par[2]+par[3]*L1+par[4]*A)
  negLL=(LLHeapPois(Y,lam,logitpi1))
  return(negLL)
}

ML.HeapPoisPG<-function(data,init.Pois){
  MLEs<-optim(init.Pois, PGmaxLL, data=data, method="BFGS")
  MLEs2<-c(inv.logit(MLEs$par[1]),MLEs$par[2:4])
  return(MLEs2)
}

PGmaxLL.mo<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  L1<-data$edu
  logitpi1=par[1]
  lam=exp(par[2]+par[3]*L1+par[4]*A)
  negLL=(LLHeapPois(Y,lam,logitpi1))
  return(negLL)
}


ML.HeapPoisPG.mo<-function(data,init.Pois){
  MLEs<-optim(init.Pois, PGmaxLL.mo, data=data, method="BFGS")
  MLEs2<-c(inv.logit(MLEs$par[1]),MLEs$par[2:4])
  return(MLEs2)
}

#DR
DRmaxLL<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  L1<-data$income
  W<-data$IPTW2
  
  logitpi1=par[1]
  lam=exp(par[2]+par[3]*A+par[4]*L1+par[5]*W)
  negLL=(LLHeapPois(Y,lam,logitpi1))
  return(negLL)
}

ML.HeapPoisDR<-function(data,init.Pois){
  MLEs<-optim(init.Pois, DRmaxLL, data=data, method="BFGS")
  MLEs2<-c(inv.logit(MLEs$par[1]),MLEs$par[2:5])
  return(MLEs2)
}

DRmaxLL.mo<-function(par,data){
  A<-data$inc
  Y<-data$CIG.heap
  L1<-data$edu
  W<-data$IPTW2
  
  logitpi1=par[1]
  lam=exp(par[2]+par[3]*A+par[4]*L1+par[5]*W)
  negLL=(LLHeapPois(Y,lam,logitpi1))
  return(negLL)
}

ML.HeapPoisDR.mo<-function(data,init.Pois){
  MLEs<-optim(init.Pois, DRmaxLL.mo, data=data, method="BFGS")
  MLEs2<-c(inv.logit(MLEs$par[1]),MLEs$par[2:5])
  return(MLEs2)
}



##########################################################################################################
#### MSM estimation functions  ####
##########################################################################################################

##wts known
geex_msm_fixed_heap_NB <- function(data, HEAP.IV, RR.IV){
  m_estimate(
    estFUN = msm_estFUN_fixed_heap_NB, 
    data   = data, 
       roots = c(HEAP.IV, RR.IV), 
       compute_roots = FALSE)
    #root_control = setup_root_control(start = c(HEAP.IV, RR.IV)))
}

msm_estFUN_fixed_heap_NB <- function(data){
  A <- data$exp
  Y <- data$out
  W <- data$wt

  function(theta){
    p<-length(theta)
    pi1=theta[1]
    odisp<-theta[2]
    lam=exp(theta[3]+theta[4]*A)
    c(ScoreMSMHeapNB(A,Y,W,lam,pi1,odisp),
      exp(theta[p-1])-theta[p])
  }}
  
  
  MSM_est_fixed_heap_NB <- function(data,exposure,outcome,weight,HEAP.IV,RR.IV) {
    data$exp<-exposure
    data$out<-outcome
    data$wt<-weight
    geex_results_fixed <- geex_msm_fixed_heap_NB(data,HEAP.IV,RR.IV)
    RR.MSM.fixed<-geex_results_fixed@estimates[length(geex_results_fixed@estimates)]
    seRR.MSM.fixed <- sqrt(geex_results_fixed@vcov[length(geex_results_fixed@estimates),length(geex_results_fixed@estimates)])
    MSM_ests_fixed<-cbind(RR.MSM.fixed,seRR.MSM.fixed)
  #  return(geex_results_fixed)
    return(MSM_ests_fixed)
  }
  

  ##wts estimated
   geex_msm_heap <- function(data, propensity_formula, HEAP.IV, RR.IV){
     e_model  <- glm(propensity_formula, data = data, family =binomial)
     models <- list(e = e_model)
    # nparms <- sum(unlist(lapply(models, function(x) length(coef(x))))) + 4
     
     m_estimate(
       estFUN = msm_estFUN_heap_NB, 
       data   = data, 
       #root_control = setup_root_control(start = c(coef(e_model), HEAP.IV, RR.IV)), 
       roots = c(coef(e_model), HEAP.IV, RR.IV), 
       compute_roots = FALSE,
       outer_args = list(models = models))
   }
   
   msm_estFUN_heap_NB <- function(data,models){
     A <- data$exp
     Y <- data$out
     Xe <- grab_design_matrix(data = data,
                              rhs_formula = grab_fixed_formula(models$e))
     e_pos  <- 1:ncol(Xe)
     e_scores  <- grab_psiFUN(models$e, data)
     lam_start<-ncol(Xe)+1
     
     function(theta){
       e <- plogis(Xe %*% theta[e_pos])
       W <-A*e^(-1)+(1-A)*(1-e)^(-1)
       pi1=theta[lam_start]
       odisp=theta[lam_start+1]
       lam=exp(theta[lam_start+2]+theta[lam_start+3]*A)
       p<-length(theta)
       
       c(e_scores(theta[e_pos]), 
         ScoreMSMHeapNB(A,Y,W,lam,pi1,odisp), 
         exp(theta[p-1])-theta[p])
     }}
  
   MSM_est_heap_NB <- function(data,exposure,outcome,propensity_formula,HEAP.IV,RR.IV) {
     data$exp<-exposure
     data$out<-outcome
     geex_results <- geex_msm_heap(data,propensity_formula,HEAP.IV,RR.IV)
     
     RR.MSM<-geex_results@estimates[length(geex_results@estimates)]
     seRR.MSM <- sqrt(geex_results@vcov[length(geex_results@estimates),length(geex_results@estimates)])
     MSM_ests<-cbind(RR.MSM,seRR.MSM)
     
     #return(geex_results)
     return(MSM_ests)
   }
  
  
  
  ##########################################################################################################
  #### Parametric g-formula ####
  ##########################################################################################################
  
   geex_PG_heap <- function(data, estfun, MOD.IV, CM0.IV, CM1.IV, RR.IV){
     
     m_estimate(
       estFUN = estfun, 
       data   = data, 
       roots = c(MOD.IV,CM1.IV,CM0.IV,RR.IV),
       compute_roots = FALSE)
    # root_control = setup_root_control(start = c(MOD.IV,CM1.IV,CM0.IV,RR.IV)))
   }
   
   estfun_gf_Pois_heap <- function(data){
     
     A<-data$inc
     L1<-data$income
     Y<-data$CIG.heap
     
     
     function(theta){
       p<-length(theta)
       pi1=theta[1]
       lam=exp(theta[2]+theta[3]*L1+theta[4]*A)
       CM0<-exp(theta[2]+theta[3]*L1)
       CM1<-exp(theta[2]+theta[3]*L1+theta[4])
       c(ScorePGHeapPois(A,L1,Y,lam,pi1),
         CM1-theta[p-2], 
         CM0-theta[p-1],
         theta[p-2]-theta[p]*theta[p-1])
     }
   }
   
   estfun_gf_Pois_heap_mis <- function(data){
     
     A<-data$inc
     L1<-data$edu
     Y<-data$CIG.heap
     
     
     function(theta){
       p<-length(theta)
       pi1=theta[1]
       lam=exp(theta[2]+theta[3]*L1+theta[4]*A)
       CM0<-exp(theta[2]+theta[3]*L1)
       CM1<-exp(theta[2]+theta[3]*L1+theta[4])
       c(ScorePGHeapPois(A,L1,Y,lam,pi1),
         CM1-theta[p-2], 
         CM0-theta[p-1],
         theta[p-2]-theta[p]*theta[p-1])
     }
   }
   
   
   PG_est_heap <- function(data, estfun, MOD.IV, CM0.IV, CM1.IV, RR.IV) {
     geex_resultspg <- geex_PG_heap(data, estfun, MOD.IV, CM0.IV, CM1.IV, RR.IV)
     
     RR.PG<-geex_resultspg@estimates[length(geex_resultspg@estimates)]
     seRR.PG <- sqrt(geex_resultspg@vcov[length(geex_resultspg@estimates),length(geex_resultspg@estimates)])
     PG_ests<-cbind(RR.PG,seRR.PG)
     
     return(PG_ests)
     #return(geex_resultspg)
   }
  
  
  
  ##########################################################################################################
  #### DR Estimation ####
  ##########################################################################################################
  
   geex_DR_heap <- function(data, estfun, propensity_formula, MOD.IV, CM0.IV, CM1.IV, RR.IV){
     
     e_model  <- glm(propensity_formula, data = data, family =binomial)
     models <- list(e = e_model)

     m_estimate(
       estFUN = estfun, 
       data   = data, 
      # root_control = setup_root_control(start = c(MOD.IV,coef(e_model), CM1.IV, CM0.IV, RR.IV)),
       roots = c(MOD.IV,coef(e_model), CM1.IV, CM0.IV, RR.IV),
       compute_roots = FALSE,
       outer_args = list(models = models))
   }
   
   estfun_dr_Pois_heap <- function(data,models){
     
     A<-data$inc
     L1<-data$income
     Y<-data$CIG.heap
     
     Xe <- grab_design_matrix(data = data,
                              rhs_formula = grab_fixed_formula(models$e))
     e_pos  <- 6:(5+ncol(Xe))
     e_scores  <- grab_psiFUN(models$e, data)
     
     function(theta){
       p<-length(theta)
       e <- plogis(Xe %*% theta[e_pos]) 
       W <-A*e^(-1)+(-1)*(1-A)*(1-e)^(-1)
       W0<-(-1)*(1-e)^(-1)
       W1<-e^(-1)
       pi1<-theta[1]
       lam=exp(theta[2]+theta[3]*A+theta[4]*L1+theta[5]*W)
       CM0<-exp(theta[2]+theta[4]*L1+theta[5]*W0)
       CM1<-exp(theta[2]+theta[3]+theta[4]*L1+theta[5]*W1)
       c(ScoreDRHeapPois(A,L1,W,Y,lam,pi1),  
         e_scores(theta[e_pos]),
         CM1-theta[p-2], 
         CM0-theta[p-1],
         theta[p]*theta[p-1]- theta[p-2])
     }
   }
   
   estfun_dr_Pois_heap_mis <- function(data,models){
     
     A<-data$inc
     L1<-data$edu
     Y<-data$CIG.heap
     
     Xe <- grab_design_matrix(data = data,
                              rhs_formula = grab_fixed_formula(models$e))
     e_pos  <- 6:(5+ncol(Xe))
     e_scores  <- grab_psiFUN(models$e, data)
     
     function(theta){
       p<-length(theta)
       e <- plogis(Xe %*% theta[e_pos]) 
       W <-A*e^(-1)+(-1)*(1-A)*(1-e)^(-1)
       W0<-(-1)*(1-e)^(-1)
       W1<-e^(-1)
       pi1<-theta[1]
       lam=exp(theta[2]+theta[3]*A+theta[4]*L1+theta[5]*W)
       CM0<-exp(theta[2]+theta[4]*L1+theta[5]*W0)
       CM1<-exp(theta[2]+theta[3]+theta[4]*L1+theta[5]*W1)
       c(ScoreDRHeapPois(A,L1,W,Y,lam,pi1),  
         e_scores(theta[e_pos]),
         CM1-theta[p-2], 
         CM0-theta[p-1],
         theta[p]*theta[p-1]- theta[p-2])
     }
   }
   
   
   DR_est_heap_Pois <- function(data,estfun,propensity_formula, MOD.IV, CM0.IV, CM1.IV, RR.IV) {
     
     geex_resultsdr <- geex_DR_heap(data,estfun,propensity_formula, MOD.IV, CM0.IV, CM1.IV, RR.IV)
     
     RR.DR<-geex_resultsdr@estimates[length(geex_resultsdr@estimates)]
     seRR.DR <- sqrt(geex_resultsdr@vcov[length(geex_resultsdr@estimates),length(geex_resultsdr@estimates)])
     DR_ests<-cbind(RR.DR,seRR.DR)
     
     return(DR_ests)
    # return(geex_resultsdr)
   }
   
  
  
  
  
  
  