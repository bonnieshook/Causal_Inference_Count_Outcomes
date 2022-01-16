#Program: 10_Poisson_Heaping_04.08.21.R
#Developed by: BES
#Purpose: Simulate data similar to WIHS, analyze using MSM, PG, and DR in geex

library(numDeriv)
library(boot)
library(geex)
library(rootSolve)
library(MASS)
library(numDeriv)

num.part<-800
num.sims<-2000
log.RR <- 0.25 
pi1<-.4 #prob of exact count

sim.results<-rep(NA, num.sims)

for (k in 1:num.sims){
  #bring in estimating functions
   source('/nas/longleaf/home/bshooksa/R/libs/00_Estimators_11.03.21.R')
   source('/nas/longleaf/home/bshooksa/R/libs/00_Estimators_Heaping_01.09.22.R')
  
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
   
   #calculate initial values for geex
   CM0.iv.naive<-sum(simdat2$IPTW*(1-simdat2$inc)*simdat2$CIG.heap)/sum(simdat2$IPTW*(1-simdat2$inc))
   CM1.iv.naive<-sum(simdat2$IPTW*(simdat2$inc)*simdat2$CIG.heap)/sum(simdat2$IPTW*(simdat2$inc))
   log.RR.iv.naive<-log(CM1.iv.naive/CM0.iv.naive)
   RR.iv.naive<-CM1.iv.naive/CM0.iv.naive

   #call estimating function, weights known
   MSM.ests.fixed.naive<-MSM_est_fixed(simdat2,simdat2$inc,simdat2$CIG.heap,simdat2$IPTW,CM0.iv.naive, CM1.iv.naive,log.RR.iv.naive,RR.iv.naive)
   colnames(MSM.ests.fixed.naive) <- c('RR.MSM.fixed.naive', 'seRR.MSM.fixed.naive') 
   
   #call estimating function, weights estimated
    MSM.ests.naive<-MSM_est(simdat2,simdat2$inc,simdat2$CIG.heap, inc ~ income, CM0.iv.naive, CM1.iv.naive,log.RR.iv.naive,RR.iv.naive)
    colnames(MSM.ests.naive) <- c('RR.MSM.naive', 'seRR.MSM.naive') 
    
    ############HEAPING#############
   ### Correct specification
   naive.nbin<-glm.nb(CIG.heap ~ 1 + inc, data = simdat2, weights=IPTW)
   init.coefs.MSM<-c(0,1/naive.nbin$theta,coef(naive.nbin))
   HeapNB.IV.MSM<-ML.HeapNBMSM(simdat2,init.coefs.MSM)
   
   RR.iv.heapNB<-exp(HeapNB.IV.MSM[4])
   MSM.ests.fixed.heap<-MSM_est_fixed_heap_NB(simdat2,simdat2$inc,simdat2$CIG.heap,simdat2$IPTW,HeapNB.IV.MSM, RR.iv.heapNB)
   colnames(MSM.ests.fixed.heap) <- c('RR.MSM.fixed.heap', 'seRR.MSM.fixed.heap') 
 
   ### weights estimated
   MSM.ests.heap<-MSM_est_heap_NB(simdat2,simdat2$inc,simdat2$CIG.heap, inc ~ income, HeapNB.IV.MSM, RR.iv.heapNB)
   colnames(MSM.ests.heap) <- c('RR.MSM.heap', 'seRR.MSM.heap') 
   
   
   ###################################
   ####### PARAMETRIC G-FORMULA ######
   ###################################
   
   ############NAIVE#############
   ###Correctly specified outcome model
   out_model.naive  <- glm(CIG.heap ~ income + inc, data = simdat2, family = poisson)
   alltrt<-simdat2
   alluntrt<-simdat2
   alltrt$inc<-1
   alluntrt$inc<-0
   CM0.iv.pg.naive<-mean(predict(out_model.naive, alluntrt, type = "response"))
   CM1.iv.pg.naive<-mean(predict(out_model.naive, alltrt, type = "response"))
   RR.iv.pg.naive<-CM1.iv.pg.naive/CM0.iv.pg.naive
   
   #call estimating function
   PG.ests.naive<-PG_est(simdat2, estfun_gf_Pois_cig, coef(out_model.naive), CM0.iv.pg.naive, CM1.iv.pg.naive, RR.iv.pg.naive)
   colnames(PG.ests.naive) <- c('RR.PG.naive', 'seRR.PG.naive') 
   
   
   ############HEAPING#############
   ####Correct specificatoin
   HeapPoisPG.IV<-ML.HeapPoisPG(simdat2,c(0,coef(out_model.naive)))

   CM0.iv.pg.heap<-mean(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income))
   CM1.iv.pg.heap<-mean(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income+HeapPoisPG.IV[4]))
   RR.iv.pg.heap<-CM1.iv.pg.heap/CM0.iv.pg.heap
   PG.ests.heap<-PG_est_heap(simdat2,estfun_gf_Pois_heap, HeapPoisPG.IV,CM0.iv.pg.heap,CM1.iv.pg.heap,RR.iv.pg.heap)
   colnames(PG.ests.heap) <- c('RR.PG.heap', 'seRR.PG.heap') 
   
   
   ###################################
   ####### DR Estimation #############
   ###################################
   
   ############NAIVE#############
   #correct specification of both models
   #calculate the initial values for geex
   CM1.iv.dr.naive<-mean(((simdat2$inc*simdat2$CIG.heap-(simdat2$inc-denom$fitted.values)*predict(out_model.naive, alltrt, type = "response")))/denom$fitted.values)
   CM0.iv.dr.naive<-mean((((1-simdat2$inc)*simdat2$CIG.heap+(simdat2$inc-denom$fitted.values)*predict(out_model.naive, alluntrt, type = "response")))/(1-denom$fitted.values))
   RR.iv.dr.naive<-CM1.iv.dr.naive/CM0.iv.dr.naive
   #call estimating function
   DR.ests.naive<-DR_est(simdat2, estfun_dr_Pois_Cig, inc ~ income, coef(out_model.naive), CM0.iv.dr.naive, CM1.iv.dr.naive, RR.iv.dr.naive)
   colnames(DR.ests.naive) <- c('RR.DR.naive', 'seRR.DR.naive') 
   
   ############HEAPING#############
   ###Correct specification
   out_model.DR.naive<-glm(CIG.heap ~ 1 + inc + income + IPTW2, data = simdat2, family = poisson)
   HeapPoisDR.IV<-ML.HeapPoisDR(simdat2,c(0,coef(out_model.DR.naive)))
   
   CM0.iv.dr.heap<-mean(exp(HeapPoisDR.IV[2]+HeapPoisDR.IV[4]*simdat2$income+HeapPoisDR.IV[5]*simdat2$noincweight*(-1)))
   CM1.iv.dr.heap<-mean(exp(HeapPoisDR.IV[2]+HeapPoisDR.IV[3]+HeapPoisDR.IV[4]*simdat2$income+HeapPoisDR.IV[5]*simdat2$incweight))
   RR.iv.dr.heap<-CM1.iv.dr.heap/CM0.iv.dr.heap
   
   #call estimating function
   DR.ests.heap<-DR_est_heap_Pois(simdat2, estfun_dr_Pois_heap, inc ~ income, HeapPoisDR.IV, CM0.iv.dr.heap, CM1.iv.dr.heap, RR.iv.dr.heap)
   colnames(DR.ests.heap) <- c('RR.DR.heap', 'seRR.DR.heap') 
   
 
  #combine results for three estimators
  simres<-cbind(MSM.ests.fixed.naive, MSM.ests.naive, MSM.ests.fixed.heap, MSM.ests.heap, PG.ests.naive, PG.ests.heap, DR.ests.naive, DR.ests.heap)
  
    
 sim.results<-rbind(sim.results,simres)
  
  
rm(list= ls()[!(ls() %in% c('num.sims','sim.results','num.part','log.RR','pi1'))])
  
}

sim.results2<-as.data.frame(sim.results[-1,])
write.csv(sim.results2, "HeapSims_cs_010922_n800.csv")