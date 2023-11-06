#Program: 01_Analysis_no_heaping.R
#Developed by: BES
#Purpose: Analyze example dataset ('noheap.csv') using IPTW, PG, and DR estimators
#Created: last updated 10.23.23

#load required libraries
library(numDeriv)
library(boot)
library(geex)
library(MASS)

#bring in example dataset: 
### outcome=SP, exposure=inc, covariates=(agevis10, sxex, druguse)
### Estimate the effect of inc on SP, assuming covariates agevis10, sxex, druguse provide conditional exchangeability
### For outcome regression approaches, the Poisson distribution is assumed. Places that are adapted for other possible distributions are noted. 
dat<-read.csv('noheap.csv')


# load functions from source files - note: places that need updating with included covariates are indicated in this program. 
source('00_Estimators_no_Heaping_10.23.23.R')


  ###################################
  ####### IPTW Estimators############
  ###################################
  
  #estimate weights
  denom <- glm(inc ~ druguse + agevis10 + sxex, data = dat, family = binomial("logit")) #update with covariates included in the weight model
  denom.probs<-denom$fitted.values
  dat2<-as.data.frame(cbind(dat,denom.probs))
  dat2$incweight<-1/dat2$denom.probs
  dat2$noincweight<-1/(1-dat2$denom.probs)
  dat2$IPTW<-ifelse(dat2$inc==1,dat2$incweight, dat2$noincweight)

  #calculate the IPTW estimator
  CM0.iv<-sum(dat2$IPTW*(1-dat2$inc)*dat2$SP)/sum(dat2$IPTW*(1-dat2$inc))
  CM1.iv<-sum(dat2$IPTW*(dat2$inc)*dat2$SP)/sum(dat2$IPTW*(dat2$inc))
  log.RR.iv<-log(CM1.iv/CM0.iv)
  RR.iv<-CM1.iv/CM0.iv

  #IPTW estimator and corresponding SE, treating weights as known (this approach is not recommended - see manuscript)
  IPTW.ests.fixed<-MSM_est_fixed(dat2,dat2$inc,dat2$SP,dat2$IPTW, CM0.iv, CM1.iv,log.RR.iv,RR.iv)
  #this function takes in the dataset name, exposure name, outcome name, weight variable name, IPTW estimated causal mean under no exposure, IPTW estimated causal mean under exposure, log of the IPTW estimated CMR, and IPTW estimated CMR
  
  #IPTW estimator and corresponding SE, treating weights estimated (this is the recommended IPTW estimator)
  IPTW.ests<-MSM_est(dat2,dat2$inc,dat2$SP, inc ~ druguse + agevis10 + sxex, CM0.iv, CM1.iv,log.RR.iv,RR.iv)
  #this function takes in the dataset name, exposure name, outcome name, assumed propensity model, IPTW estimated causal mean under no exposure, IPTW estimated causal mean under exposure, log of the IPTW estimated CMR, and IPTW estimated CMR
  
  
  #############################################
  ####### PARAMETRIC G-FORMULA ESTIMATORS######
  #############################################
  
  #fit the outcome regression model (assuming Poisson distribution)
  out_model_Poisson  <- glm(SP ~ 1 + inc + druguse + sxex +  agevis10 , data = dat2, family = poisson) #update with covariates included in the outcome model, and switch distribution (as appropriate)

  #calculate the parametric g-formula estimator
  alltrt<-dat2
  alluntrt<-dat2
  alltrt$inc<-1
  alluntrt$inc<-0
  CM0.iv.pg<-mean(predict(out_model_Poisson, alluntrt, type = "response"))
  CM1.iv.pg<-mean(predict(out_model_Poisson, alltrt, type = "response"))
  RR.iv.pg<-CM1.iv.pg/CM0.iv.pg
  
  #Parametric g-formula estimator and corresponding SE, assuming a Poisson distribution
  PG.ests<-PG_est(dat2, estfun_gf_Pois, coef(out_model_Poisson), CM0.iv.pg, CM1.iv.pg, RR.iv.pg)
  #this function takes in the dataset name, name of the estimating function from the 00_Estimators_no_Heaping_10.23.23 program (estfun_gf_Pois, estfun_gf_NBin, estfun_gf_ZIP, or estfun_gf_ZINB), the estimated coefficients of the outcome regression model based on the assumed distribution, PG estimated causal mean under no exposure, PG estimated causal mean under exposure, and Parametric g-formula estimated CMR
  
  
   
  ###################################
  ####### DR Estimation #############
  ###################################
  
  #calculate the DR estimates of the causal means based on predicted values and estimated propensity scores, and the estimated CMR
  CM1.iv.dr<-mean(((dat2$inc*dat2$SP-(dat2$inc-denom$fitted.values)*predict(out_model_Poisson, alltrt, type = "response")))/denom$fitted.values)
  CM0.iv.dr<-mean((((1-dat2$inc)*dat2$SP+(dat2$inc-denom$fitted.values)*predict(out_model_Poisson, alluntrt, type = "response")))/(1-denom$fitted.values))
  RR.iv.dr<-CM1.iv.dr/CM0.iv.dr
  
  #DR estimator and corresponding SE, assuming a Poisson distribution
  DR.ests<-DR_est(dat2, estfun_dr_Pois, inc ~ druguse + agevis10 + sxex, coef(out_model_Poisson), CM0.iv.dr, CM1.iv.dr, RR.iv.dr)
  #this function takes in the dataset name, name of the estimating function from the 00_Estimators_no_Heaping_10.23.23 program (estfun_dr_Pois, estfun_dr_NBin, estfun_dr_ZIP, or estfun_dr_ZINB), the assumed propensity model, the estimated coefficients of the outcome regression model based on the assumed distribution, DR estimated causal mean under no exposure, DR estimated causal mean under exposure, and DR estimated CMR
  
  
  