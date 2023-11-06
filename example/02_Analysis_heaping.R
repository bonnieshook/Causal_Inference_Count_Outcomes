#Program: 02_Analysis_heaping.R
#Developed by: BES
#Purpose: Analyze example dataset ('heaping.csv') using IPTW, PG, and DR estimators
#Created: last updated 10.23.23

#load required libraries
library(numDeriv)
library(boot)
library(geex)
library(MASS)

#bring in example dataset: 
### outcome(heaped)=CIG.heap, exposure=inc, covariates=(income)
### Estimate the effect of inc on cigarettes, assuming covariate income provides conditional exchangeability
### For outcome regression approaches, the Poisson distribution is assumed. The NB is assumed for marginal heaping models.
dat<-read.csv('heaping.csv')

#assign heaping interval
dat$interval<-ceiling((dat$CIG.heap-4.9)/10)+1
dat$interval.col<-ifelse(dat$interval>2,2,dat$interval) #collapse anyone with a count > 14 into interval 2



# load functions from source files - note: places that need updating with included covariates are indicated in this program. 
source('00_Estimators_Heaping_10.23.23.R')


  #########################################
  ####### IPTW Estimators, HCAR############
  #########################################
  
  #estimate weights
  denom <- glm(inc ~ income, data = dat, family = binomial("logit"))
  denom.probs<-denom$fitted.values
  dat2<-as.data.frame(cbind(dat,denom.probs))
  dat2$incweight<-1/dat2$denom.probs
  dat2$noincweight<-1/(1-dat2$denom.probs)
  dat2$IPTW<-ifelse(dat2$inc==1,dat2$incweight, dat2$noincweight)
  
  ###Estimate E(Y_h^a)
  CM0.iv.naive<-sum(dat2$IPTW*(1-dat2$inc)*dat2$CIG.heap)/sum(dat2$IPTW*(1-dat2$inc))
  CM1.iv.naive<-sum(dat2$IPTW*(dat2$inc)*dat2$CIG.heap)/sum(dat2$IPTW*(dat2$inc))

  ###Estimate E(h(Y^a)) 
  dat2$CIG.heapall<-round(dat2$CIG.heap+0.00001,-1) #round all reported counts
  CM0.hY<-sum(dat2$IPTW*(1-dat2$inc)*dat2$CIG.heapall)/sum(dat2$IPTW*(1-dat2$inc))
  CM1.hY<-sum(dat2$IPTW*(dat2$inc)*dat2$CIG.heapall)/sum(dat2$IPTW*(dat2$inc))
  
  #estimate pi under HCAR
  naive.nbin<-glm.nb(CIG.heap ~ 1, data = dat2) #get smart starting values to find MLEs
  init.coefs.MSM<-c(0,1/naive.nbin$theta,coef(naive.nbin))
  HeapNB.IV<-ML.Heap.marginal.HCAR(dat2,init.coefs.MSM)
  
  #Estimate the causal means and the CMR using IPTW HCAR estimator
  CM0.IPTW.heap.HCAR<-((HeapNB.IV[1])^(-1))*(CM0.iv.naive-(1-(HeapNB.IV[1]))*CM0.hY)
  CM1.IPTW.heap.HCAR<-((HeapNB.IV[1])^(-1))*(CM1.iv.naive-(1-(HeapNB.IV[1]))*CM1.hY)
  CMR.IPTW.HCAR<-CM1.IPTW.heap.HCAR/CM0.IPTW.heap.HCAR
  
  #store roots to feed into geex when estimating the SE
  Heap.IPTW.IV.HCAR<-c(HeapNB.IV,CM1.iv.naive,CM0.iv.naive,CM1.hY,CM0.hY,CM1.IPTW.heap.HCAR,CM0.IPTW.heap.HCAR,CMR.IPTW.HCAR)
  
  #IPTW estimator and corresponding SE, treating weights as known (this approach is not recommended - see manuscript)
  IPTW.ests.fixed.heap.HCAR<-IPTW_est_fixed_heap_NB.HCAR(dat2,dat2$inc,dat2$CIG.heap,dat2$IPTW,Heap.IPTW.IV.HCAR)
  #this function takes in the dataset name, exposure name, (heaped) outcome name, weight variable name, and IPTW estimated CMR
  
  #IPTW estimator and corresponding SE, treating weights estimated (this is the recommended IPTW estimator)
  IPTW.ests.heap.HCAR<-IPTW_est_heap_NB.HCAR(dat2,dat2$inc,dat2$CIG.heap, inc ~ income, Heap.IPTW.IV.HCAR)
  #this function takes in the dataset name, exposure name, (heaped) outcome name, assumed propensity model, and IPTW estimated CMR
    

  ###########################################################
  ####### PARAMETRIC G-FORMULA ESTIMATORS, HCAR and IH ######
  ###########################################################
  
  #fit naive outcome regression model (assuming Poisson distribution) to use as initial values when finding MLEs in heaping model
  out_model.naive  <- glm(CIG.heap ~ income + inc, data = dat2, family = poisson)
  
  #### Fit HCAR heaping model
  HeapPoisPG.IV<-ML.HeapPoisPG.HCAR(dat2,c(0,coef(out_model.naive)))
  
  # Estimate causal means and CMR using PG HCAR estimator
  CM0.iv.pg.heap<-mean(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*dat2$income))
  CM1.iv.pg.heap<-mean(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*dat2$income+HeapPoisPG.IV[4]))
  RR.iv.pg.heap<-CM1.iv.pg.heap/CM0.iv.pg.heap
  
  #store roots to feed into geex when estimating the SE
  Heap.PG.IV.HCAR<-c(HeapPoisPG.IV,CM1.iv.pg.heap,CM0.iv.pg.heap,RR.iv.pg.heap)
  
  #Parametric g-formula HCAR estimator and corresponding SE, assuming a Poisson distribution
  PG.ests.heap.HCAR<-PG_est_heap.HCAR(dat2,estfun_gf_Pois_heap.HCAR, Heap.PG.IV.HCAR)
  #this function takes in the dataset name, name of the estimating function from the 00_Estimators_Heaping_10.23.23.R program, and the roots
  
  
  #### Fit informative heaping model
  HeapPoisPG.IV.IH<-ML.HeapPoisPG.IH(dat2,c(0,0,coef(out_model.naive)))
  
  # Estimate causal means and CMR using PG IH estimator
  CM0.iv.pg.heap.IH<-mean(exp(HeapPoisPG.IV.IH[3]+ HeapPoisPG.IV.IH[4]*dat2$income))
  CM1.iv.pg.heap.IH<-mean(exp(HeapPoisPG.IV.IH[3]+ HeapPoisPG.IV.IH[4]*dat2$income+ HeapPoisPG.IV.IH[5]))
  RR.iv.pg.heap.IH<-CM1.iv.pg.heap.IH/CM0.iv.pg.heap.IH
  
  #store roots to feed into geex when estimating the SE
  Heap.PG.IV.IH<-c(HeapPoisPG.IV.IH,CM1.iv.pg.heap.IH,CM0.iv.pg.heap.IH,RR.iv.pg.heap.IH)
  
  #Parametric g-formula IH estimator and corresponding SE, assuming a Poisson distribution
  PG.ests.heap.IH<-PG_est_heap.IH(dat2,estfun_gf_Pois_heap.IH, Heap.PG.IV.IH)
  #this function takes in the dataset name, name of the estimating function from the 00_Estimators_Heaping_10.23.23.R program, and the roots
  
  
  
  #########################################
  ####### DR Estimation, HCAR #############
  #########################################
  
  #compute hat(Y_h^a) and hat(h(Y^a)) for each participant
  dat2$hat.Y0<-(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*dat2$income))
  dat2$hat.Y1<-(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*dat2$income+HeapPoisPG.IV[4]))
  dat2$hat.h.Y0<-round(dat2$hat.Y0+0.00001,-1)
  dat2$hat.h.Y1<-round(dat2$hat.Y1+0.00001,-1)
  dat2$hat.Yh0<- Heap.PG.IV.HCAR[1]*dat2$hat.Y0+(1- Heap.PG.IV.HCAR[1])*dat2$hat.h.Y0
  dat2$hat.Yh1<- Heap.PG.IV.HCAR[1]*dat2$hat.Y1+(1- Heap.PG.IV.HCAR[1])*dat2$hat.h.Y1 
  
  #compute hat(E(Y_h^a))
  DR.EYh1<-mean((denom$fitted.values^(-1))*(dat2$inc*dat2$CIG.heap-(dat2$inc-denom$fitted.values)*dat2$hat.Yh1))
  DR.EYh0<-mean(((1-denom$fitted.values)^(-1))*((1-dat2$inc)*dat2$CIG.heap+(dat2$inc-denom$fitted.values)*dat2$hat.Yh0))
  
  #compute hat(E(h(Y^a)))
  DR.Eh.Y1<-mean((denom$fitted.values^(-1))*(dat2$inc*dat2$CIG.heapall-(dat2$inc-denom$fitted.values)*dat2$hat.h.Y1))
  DR.Eh.Y0<-mean(((1-denom$fitted.values)^(-1))*((1-dat2$inc)*dat2$CIG.heapall+(dat2$inc-denom$fitted.values)*dat2$hat.h.Y0))
  
  # Estimate causal means and CMR using DR HCAR estimator
  CM1.iv.dr.heap<-(Heap.IPTW.IV.HCAR[1]^(-1))*(DR.EYh1-(1-Heap.IPTW.IV.HCAR[1])*DR.Eh.Y1)
  CM0.iv.dr.heap<-(Heap.IPTW.IV.HCAR[1]^(-1))*(DR.EYh0-(1-Heap.IPTW.IV.HCAR[1])*DR.Eh.Y0)
  RR.iv.dr.heap<-CM1.iv.dr.heap/CM0.iv.dr.heap
  
  #store roots to feed into geex when estimating the SE
  DR.IV.HCAR<-c(HeapNB.IV,HeapPoisPG.IV,DR.EYh1,DR.EYh0,DR.Eh.Y1,DR.Eh.Y0,CM1.iv.dr.heap,CM0.iv.dr.heap,RR.iv.dr.heap)
  
  #call estimating function
  DR.ests.heap.HCAR<-DR_est_heap.HCAR(dat2, estfun_dr_heap.HCAR, inc ~ income, DR.IV.HCAR)
  #this function takes in the dataset name, name of the estimating function from the 00_Estimators_Heaping_10.23.23.R program,the assumed propensity model, and the roots
  
  