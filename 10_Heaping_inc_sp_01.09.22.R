#Program: 10_Heaping_inc_sp_09.02.21.R
#Developed by: BES
#Purpose: Simulate data similar to WIHS with data heaping, analyze using MSM, PG, and DR in geex

library(numDeriv)
library(boot)
library(geex)
library(rootSolve)
library(MASS)
library(numDeriv)

num.part<-800
num.sims<-2000
log.RR <- 0.25 
pi1<-.4 

sim.results<-rep(NA, num.sims)

for (k in 1:num.sims){
  #bring in estimating functions
   source('/nas/longleaf/home/bshooksa/R/libs/00_Estimators_11.03.21.R')
   source('/nas/longleaf/home/bshooksa/R/libs/00_Estimators_Heaping_01.09.22.R')
 
  skip_to_next=FALSE
  skip_to_next2=FALSE
  
  set.seed(020720+k)
  
  #generage confounder X 
  income<-as.vector(log(rgamma(num.part, shape=5, rate=0.5)))

  #generate participant id
  part.id<-seq(1:num.part)
  
  simdat1<-as.data.frame(cbind(part.id,income))
  
  #incarceration
  simdat1$incrate<-1-inv.logit(-.8+0.65*simdat1$income)
  simdat1$inc<-rbinom(num.part,1,simdat1$incrate)

  #create potential outcomes 
  simdat1$lambda0<-exp(simdat1$income-0.9)
  simdat1$lambda1<-exp(simdat1$income-0.9+log.RR)
  simdat1$CIG0<-rpois(num.part,simdat1$lambda0)
  simdat1$CIG1<-rpois(num.part,simdat1$lambda1)
  
  #assign rounding prob and heaping POs (rounding to nearest 10)
  simdat1$rpt.exact<-rbinom(num.part,1,pi1)
  
  simdat1$CIG0.heap<-ifelse(simdat1$rpt.exact==1,simdat1$CIG0,round(simdat1$CIG0+0.00001,-1))
  simdat1$CIG1.heap<-ifelse(simdat1$rpt.exact==1,simdat1$CIG1,round(simdat1$CIG1+0.00001,-1))
  
  
  #assign observed outcome
  simdat1$CIG.true<-ifelse(simdat1$inc==0,simdat1$CIG0,simdat1$CIG1)
  simdat1$CIG.heap<-ifelse(simdat1$inc==0,simdat1$CIG0.heap,simdat1$CIG1.heap)
  
  #create misspecified confounder
  simdat1$edu<-inv.logit(-3+1*simdat1$income+2*rnorm(num.part,0,1))

  
  ###################################
  ####### MSM with IPTW ############
  ###################################
  
  ############NAIVE#############
  #Weights
   denom <- glm(inc ~ income, data = simdat1, family = binomial("logit"))
   denom.probs<-denom$fitted.values
   simdat2<-as.data.frame(cbind(simdat1,denom.probs))
   simdat2$incweight<-1/simdat2$denom.probs
   simdat2$noincweight<-1/(1-simdat2$denom.probs)
   simdat2$IPTW<-ifelse(simdat2$inc==1,simdat2$incweight, simdat2$noincweight)
   simdat2$IPTW2<-ifelse(simdat2$inc==1,simdat2$incweight, (-1*simdat2$noincweight)) 
   
  ############HEAPING#############

   #Incorrectly specified weight models
   #Weights
   denom.m <- glm(inc ~ edu, data = simdat1, family = binomial("logit"))
   denom.probs.m<-denom.m$fitted.values
   simdat2.m<-as.data.frame(cbind(simdat1,denom.probs.m))
   simdat2.m$incweight<-1/simdat2.m$denom.probs
   simdat2.m$noincweight<-1/(1-simdat2.m$denom.probs)
   simdat2.m$IPTW<-ifelse(simdat2.m$inc==1,simdat2.m$incweight, simdat2.m$noincweight)
   simdat2.m$IPTW2<-ifelse(simdat2.m$inc==1,simdat2.m$incweight, (-1*simdat2.m$noincweight))
   
   naive.nbin.mw<-glm.nb(CIG.heap ~ 1 + inc, data = simdat2.m, weights=IPTW)
   init.coefs.MSM.mw<-c(0,1/naive.nbin.mw$theta,coef(naive.nbin.mw))
   HeapNB.IV.MSM.mw<-ML.HeapNBMSM(simdat2.m,init.coefs.MSM.mw)
   RR.iv.heapNB.mw<-exp(HeapNB.IV.MSM.mw[4])
   MSM.ests.fixed.heap.mw<-MSM_est_fixed_heap_NB(simdat2.m,simdat2.m$inc,simdat2.m$CIG.heap,simdat2.m$IPTW,HeapNB.IV.MSM.mw, RR.iv.heapNB.mw)
   colnames(MSM.ests.fixed.heap.mw) <- c('RR.MSM.fixed.heap.mw', 'seRR.MSM.fixed.heap.mw') 
   
   ### weights estimated
   MSM.ests.heap.mw<-MSM_est_heap_NB(simdat2.m,simdat2.m$inc,simdat2.m$CIG.heap, inc ~ edu, HeapNB.IV.MSM.mw, RR.iv.heapNB.mw)
   colnames(MSM.ests.heap.mw) <- c('RR.MSM.heap.mw', 'seRR.MSM.heap.mw') 
   
   ###################################
   ####### PARAMETRIC G-FORMULA ######
   ###################################
 
   
   ############HEAPING#############

   ###Incorrectly specified outcome model
   out_model.naive.mo  <- glm(CIG.heap ~ edu + inc, data = simdat2.m, family = poisson)
   HeapPoisPG.IV.mo<-ML.HeapPoisPG.mo(simdat2.m,c(0,coef(out_model.naive.mo)))
   CM0.iv.pg.heap.mo<-mean(exp(HeapPoisPG.IV.mo[2]+HeapPoisPG.IV.mo[3]*simdat2.m$edu))
   CM1.iv.pg.heap.mo<-mean(exp(HeapPoisPG.IV.mo[2]+HeapPoisPG.IV.mo[3]*simdat2.m$edu+HeapPoisPG.IV.mo[4]))
   RR.iv.pg.heap.mo<-CM1.iv.pg.heap.mo/CM0.iv.pg.heap.mo
   PG.ests.heap.mo<-PG_est_heap(simdat2.m,estfun_gf_Pois_heap_mis, HeapPoisPG.IV.mo,CM0.iv.pg.heap.mo,CM1.iv.pg.heap.mo,RR.iv.pg.heap.mo)
   colnames(PG.ests.heap.mo) <- c('RR.PG.heap.mo', 'seRR.PG.heap.mo') 
   
   
   ###################################
   ####### DR Estimation #############
   ###################################
   

   ############HEAPING#############
 
   ###Incorrectly specified models
   #correct specification of outcome model, incorect specification of weight model
   out_model.DR.naive.mw<-glm(CIG.heap ~ 1 + inc + income + IPTW2, data = simdat2.m, family = poisson)
   HeapPoisDR.IV.mw<-ML.HeapPoisDR(simdat2.m,c(0,coef(out_model.DR.naive.mw)))
   
   CM0.iv.dr.heap.mw<-mean(exp(HeapPoisDR.IV.mw[2]+HeapPoisDR.IV.mw[4]*simdat2.m$income+HeapPoisDR.IV.mw[5]*simdat2.m$noincweight*(-1)))
   CM1.iv.dr.heap.mw<-mean(exp(HeapPoisDR.IV.mw[2]+HeapPoisDR.IV.mw[3]+HeapPoisDR.IV.mw[4]*simdat2.m$income+HeapPoisDR.IV.mw[5]*simdat2.m$incweight))
   RR.iv.dr.heap.mw<-CM1.iv.dr.heap.mw/CM0.iv.dr.heap.mw
   
   tryCatch(DR.ests.heap.mw<-DR_est_heap_Pois(simdat2.m, estfun_dr_Pois_heap, inc ~ edu, HeapPoisDR.IV.mw, CM0.iv.dr.heap.mw, CM1.iv.dr.heap.mw, RR.iv.dr.heap.mw), error = function(e) { skip_to_next <<- TRUE})
   if(skip_to_next==TRUE) {
     DR.ests.heap.mw<-matrix(c(NA,NA),nrow=1)
     skip_to_next<-FALSE}
   
   colnames(DR.ests.heap.mw) <- c('RR.DR.heap.mw', 'seRR.DR.heap.mw')
   
   
   #correct specification of weight model, incorect specification of outcome model
   out_model.DR.naive.mo<-glm(CIG.heap ~ 1 + inc + edu + IPTW2, data = simdat2, family = poisson)
   HeapPoisDR.IV.mo<-ML.HeapPoisDR.mo(simdat2,c(0,coef(out_model.DR.naive.mo)))
   
   CM0.iv.dr.heap.mo<-mean(exp(HeapPoisDR.IV.mo[2]+HeapPoisDR.IV.mo[4]*simdat2$edu+HeapPoisDR.IV.mo[5]*simdat2$noincweight*(-1)))
   CM1.iv.dr.heap.mo<-mean(exp(HeapPoisDR.IV.mo[2]+HeapPoisDR.IV.mo[3]+HeapPoisDR.IV.mo[4]*simdat2$edu+HeapPoisDR.IV.mo[5]*simdat2$incweight))
   RR.iv.dr.heap.mo<-CM1.iv.dr.heap.mo/CM0.iv.dr.heap.mo
   
   tryCatch(DR.ests.heap.mo<-DR_est_heap_Pois(simdat2, estfun_dr_Pois_heap_mis, inc ~ income, HeapPoisDR.IV.mo, CM0.iv.dr.heap.mo, CM1.iv.dr.heap.mo, RR.iv.dr.heap.mo), error = function(e) { skip_to_next <<- TRUE})
   if(skip_to_next==TRUE) {
     DR.ests.heap.mo<-matrix(c(NA,NA),nrow=1)
     skip_to_next<-FALSE}
   
   colnames(DR.ests.heap.mo) <- c('RR.DR.heap.mo', 'seRR.DR.heap.mo') 
   
   #misspecify both models
   out_model.DR.naive.mb<-glm(CIG.heap ~ 1 + inc + edu + IPTW2, data = simdat2.m, family = poisson)
   tryCatch(HeapPoisDR.IV.mb<-ML.HeapPoisDR.mo(simdat2.m,c(0,coef(out_model.DR.naive.mb))), error = function(e) { skip_to_next <<- TRUE})
   if(skip_to_next==FALSE){ 
      CM0.iv.dr.heap.mb<-mean(exp(HeapPoisDR.IV.mb[2]+HeapPoisDR.IV.mb[4]*simdat2.m$edu+HeapPoisDR.IV.mb[5]*simdat2.m$noincweight*(-1)))
      CM1.iv.dr.heap.mb<-mean(exp(HeapPoisDR.IV.mb[2]+HeapPoisDR.IV.mb[3]+HeapPoisDR.IV.mb[4]*simdat2.m$edu+HeapPoisDR.IV.mb[5]*simdat2.m$incweight))
      RR.iv.dr.heap.mb<-CM1.iv.dr.heap.mb/CM0.iv.dr.heap.mb
   tryCatch(DR.ests.heap.mb<-DR_est_heap_Pois(simdat2.m, estfun_dr_Pois_heap_mis, inc ~ edu, HeapPoisDR.IV.mb, CM0.iv.dr.heap.mb, CM1.iv.dr.heap.mb, RR.iv.dr.heap.mb), error = function(e) { skip_to_next2 <<- TRUE})}
   if(skip_to_next==TRUE | skip_to_next2==TRUE) {
      DR.ests.heap.mb<-matrix(c(NA,NA),nrow=1)
      }
   colnames(DR.ests.heap.mb) <- c('RR.DR.heap.mb', 'seRR.DR.heap.mb') 
   
   
  #combine results for three estimators
  simres<-cbind(MSM.ests.fixed.heap.mw, MSM.ests.heap.mw, PG.ests.heap.mo, DR.ests.heap.mw, DR.ests.heap.mo, DR.ests.heap.mb)
  sim.results<-rbind(sim.results,simres)
  
  
rm(list= ls()[!(ls() %in% c('num.sims','sim.results','num.part','log.RR','pi1'))])
  
}

sim.results2<-as.data.frame(sim.results[-1,])
write.csv(sim.results2, "HeapSims_IS_010922_n800.csv")