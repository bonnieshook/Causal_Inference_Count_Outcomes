#Program: 05_Heaping_inc_sp_HCAR.R
#Developed by: BES
#Purpose: Simulate data similar to WIHS with data heaping, analyze using MSM, PG, and DR in geex

library(numDeriv)
library(boot)
library(geex)
library(rootSolve)
library(MASS)
library(numDeriv)

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
subdir_name <- as.character(args[1])
num.part <- as.numeric(args[2])
num.sims<-as.numeric(args[3])

user_home_dir <- " " # enter name for top-level directory for the project
dir_path <- paste(user_home_dir,subdir_name,"/",sep="")

# load functions from source files, 
# all of which are in the top-level user_home_dir
source('00_Estimators_Heaping_07.17.23.R')
source('00_Estimators_11.03.21.R')

sim <- Sys.getenv("SLURM_ARRAY_TASK_ID") 

setwd(dir_path)

log.RR <- 0.25 
 
  skip_to_next=FALSE
  skip_to_next2=FALSE
  
#generate the sample
  
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
  
  #assign observed outcome
  simdat1$CIG.true<-ifelse(simdat1$inc==0,simdat1$CIG0,simdat1$CIG1)
  
  #set probability of reporting an exact count. Under HCAR, pi1=pi2
  pi1<-0.4
  pi2<-0.4
  
  
  #find cut points in data (5, 15, 25,...) and assign probability of reporting exactly within each interval
  max.cp<-floor((max(simdat1$CIG.true)-5)/10)*10+5
  cps<-seq(from=5,to=max.cp,by=10)
  num.cp.ints<-(max.cp-5)/10+2
  prob.exact<-c(pi1,pi2)
  
  #assign rounding prob and heaping POs and observed heaped outcome (rounding to nearest 10)
  simdat1$interval<-ceiling((simdat1$CIG.true-4.9)/10)+1
  simdat1$interval.col<-ifelse(simdat1$interval>2,2,simdat1$interval) #collapse anyone with a count > 14 into interval 2
  simdat1$prob.exact<-prob.exact[simdat1$interval.col]
  simdat1$rpt.exact<-rbinom(num.part,1,simdat1$prob.exact)
  
  simdat1$CIG0.heap<-ifelse(simdat1$rpt.exact==1,simdat1$CIG0,round(simdat1$CIG0+0.00001,-1))
  simdat1$CIG1.heap<-ifelse(simdat1$rpt.exact==1,simdat1$CIG1,round(simdat1$CIG1+0.00001,-1))
  simdat1$CIG.heap<-ifelse(simdat1$inc==0,simdat1$CIG0.heap,simdat1$CIG1.heap)
  simdat1$CIG.heapall<-round(simdat1$CIG.heap+0.00001,-1)
  
  #create misspecified confounder
  simdat1$edu<-inv.logit(-3+1*simdat1$income+2*rnorm(num.part,0,1))
  
  #calculate the proportion of the sample in each interval
  prop.int1<-nrow(simdat1[simdat1$interval.col==1,])/num.part
  prop.int2<-nrow(simdat1[simdat1$interval.col==2,])/num.part
  
  
  ###################################
  ####### IPTW ######################
  ###################################
  
  ############NAIVE#############
  #Weights
   denom <- glm(inc ~ income, data = simdat1, family = binomial("logit"))
   denom.probs<-denom$fitted.values
   simdat2<-as.data.frame(cbind(simdat1,denom.probs))
   simdat2$incweight<-1/simdat2$denom.probs
   simdat2$noincweight<-1/(1-simdat2$denom.probs)
   simdat2$IPTW<-ifelse(simdat2$inc==1,simdat2$incweight, simdat2$noincweight)

   
  ############HEAPING#############

   #Incorrectly specified weight models
   #Weights
   denom.m <- glm(inc ~ edu, data = simdat1, family = binomial("logit"))
   denom.probs.m<-denom.m$fitted.values
   simdat2.m<-as.data.frame(cbind(simdat1,denom.probs.m))
   simdat2.m$incweight<-1/simdat2.m$denom.probs
   simdat2.m$noincweight<-1/(1-simdat2.m$denom.probs)
   simdat2.m$IPTW<-ifelse(simdat2.m$inc==1,simdat2.m$incweight, simdat2.m$noincweight)

   ###Estimate E(h(Y^a)) (Note: E(Y_h^a) are CM0.iv.naive and CM1.iv.naive)
   CM0.iv.naive<-sum(simdat2.m$IPTW*(1-simdat2.m$inc)*simdat2.m$CIG.heap)/sum(simdat2.m$IPTW*(1-simdat2.m$inc))
   CM1.iv.naive<-sum(simdat2.m$IPTW*(simdat2.m$inc)*simdat2.m$CIG.heap)/sum(simdat2.m$IPTW*(simdat2.m$inc))   
   CM0.hY<-sum(simdat2.m$IPTW*(1-simdat2.m$inc)*simdat2.m$CIG.heapall)/sum(simdat2.m$IPTW*(1-simdat2.m$inc))
   CM1.hY<-sum(simdat2.m$IPTW*(simdat2.m$inc)*simdat2.m$CIG.heapall)/sum(simdat2.m$IPTW*(simdat2.m$inc))
   
   #estimate pi under HCAR
   naive.nbin<-glm.nb(CIG.heap ~ 1, data = simdat2.m)
   init.coefs.MSM<-c(0,1/naive.nbin$theta,coef(naive.nbin))
   HeapNB.IV<-ML.Heap.marginal.HCAR(simdat2.m,init.coefs.MSM)
   
   #calculate IPTW estimate of CMR
   CM0.IPTW.heap<-((HeapNB.IV[1])^(-1))*(CM0.iv.naive-(1-(HeapNB.IV[1]))*CM0.hY)
   CM1.IPTW.heap<-((HeapNB.IV[1])^(-1))*(CM1.iv.naive-(1-(HeapNB.IV[1]))*CM1.hY)
   CMR.IPTW<-CM1.IPTW.heap/CM0.IPTW.heap
   
   Heap.IPTW.IV<-c(HeapNB.IV,CM1.iv.naive,CM0.iv.naive,CM1.hY,CM0.hY,CM1.IPTW.heap,CM0.IPTW.heap,CMR.IPTW)
   
   ### weights fixed
   IPTW.ests.fixed.heap.HCAR.mw<-IPTW_est_fixed_heap_NB.HCAR(simdat2.m,simdat2.m$inc,simdat2.m$CIG.heap,simdat2.m$IPTW,Heap.IPTW.IV)
   colnames(IPTW.ests.fixed.heap.HCAR.mw) <- c('RR.IPTW.fixed.heap.HCAR.mw', 'seRR.IPTW.fixed.heap.HCAR.mw') 
   
   ### weights estimated
   IPTW.ests.heap.HCAR.mw<-IPTW_est_heap_NB.HCAR(simdat2.m,simdat2.m$inc,simdat2.m$CIG.heap, inc ~ edu, Heap.IPTW.IV)
   colnames(IPTW.ests.heap.HCAR.mw) <- c('RR.IPTW.heap.HCAR.mw', 'seRR.IPTW.heap.HCAR.mw') 
   
   
   ###################################
   ####### PARAMETRIC G-FORMULA ######
   ###################################
 
   ############HEAPING#############

   ###HCAR
   ###Incorrectly specified outcome model
   out_model.naive.mo  <- glm(CIG.heap ~ edu + inc, data = simdat2.m, family = poisson)
   HeapPoisPG.IV.mo<-ML.HeapPoisPG.mo.HCAR(simdat2.m,c(0,coef(out_model.naive.mo)))
   CM0.iv.pg.heap.mo<-mean(exp(HeapPoisPG.IV.mo[2]+HeapPoisPG.IV.mo[3]*simdat2.m$edu))
   CM1.iv.pg.heap.mo<-mean(exp(HeapPoisPG.IV.mo[2]+HeapPoisPG.IV.mo[3]*simdat2.m$edu+HeapPoisPG.IV.mo[4]))
   RR.iv.pg.heap.mo<-CM1.iv.pg.heap.mo/CM0.iv.pg.heap.mo
   
   Heap.PG.IV.mo<-c(HeapPoisPG.IV.mo,CM1.iv.pg.heap.mo,CM0.iv.pg.heap.mo,RR.iv.pg.heap.mo)
   
   PG.ests.heap.HCAR.mo<-PG_est_heap.HCAR(simdat2,estfun_gf_Pois_heap_mis.HCAR, Heap.PG.IV.mo)
   colnames(PG.ests.heap.HCAR.mo) <- c('RR.PG.heap.HCAR.mo', 'seRR.PG.heap.HCAR.mo') 
   
   ###IH
   ###Incorrectly specified outcome model
   HeapPoisPG.IV.IH.mo<-ML.HeapPoisPG.mo.IH(simdat2.m,c(0,0,coef(out_model.naive.mo)))
   CM0.iv.pg.heap.IH.mo<-mean(exp(HeapPoisPG.IV.IH.mo[3]+ HeapPoisPG.IV.IH.mo[4]*simdat2.m$edu))
   CM1.iv.pg.heap.IH.mo<-mean(exp(HeapPoisPG.IV.IH.mo[3]+ HeapPoisPG.IV.IH.mo[4]*simdat2.m$edu+ HeapPoisPG.IV.IH.mo[5]))
   RR.iv.pg.heap.IH.mo<-CM1.iv.pg.heap.IH.mo/CM0.iv.pg.heap.IH.mo
   
   Heap.PG.IV.IH.mo<-c(HeapPoisPG.IV.IH.mo,CM1.iv.pg.heap.IH.mo,CM0.iv.pg.heap.IH.mo,RR.iv.pg.heap.IH.mo)
   
   PG.ests.heap.IH.mo<-PG_est_heap.IH(simdat2,estfun_gf_Pois_heap_mis.IH, Heap.PG.IV.IH.mo)
   colnames(PG.ests.heap.IH.mo) <- c('RR.PG.heap.IH.mo', 'seRR.PG.heap.IH.mo') 
   
   
   ###################################
   ####### DR Estimation #############
   ###################################
   

   ############HEAPING#############
   
   ###Incorrectly specified models

   ##correct specification of outcome model, incorrect specification of weight model

   #compute hat(Y_h^a) and hat(h(Y^a)) for each participant
   out_model.naive  <- glm(CIG.heap ~ income + inc, data = simdat2, family = poisson)
   HeapPoisPG.IV<-ML.HeapPoisPG.HCAR(simdat2,c(0,coef(out_model.naive)))
   
   simdat2$hat.Y0<-(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income))
   simdat2$hat.Y1<-(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income+HeapPoisPG.IV[4]))
   simdat2$hat.h.Y0<-round(simdat2$hat.Y0+0.00001,-1)
   simdat2$hat.h.Y1<-round(simdat2$hat.Y1+0.00001,-1)
   simdat2$hat.Yh0<-HeapPoisPG.IV[1]*simdat2$hat.Y0+(1-HeapPoisPG.IV[1])*simdat2$hat.h.Y0
   simdat2$hat.Yh1<-HeapPoisPG.IV[1]*simdat2$hat.Y1+(1-HeapPoisPG.IV[1])*simdat2$hat.h.Y1 
   
   #compute hat(E(Y_h^a))
   DR.EYh1.mw<-mean((denom.m$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heap-(simdat2$inc-denom.m$fitted.values)*simdat2$hat.Yh1))
   DR.EYh0.mw<-mean(((1-denom.m$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heap+(simdat2$inc-denom.m$fitted.values)*simdat2$hat.Yh0))
   
   #compute hat(E(h(Y^a)))
   DR.Eh.Y1.mw<-mean((denom.m$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heapall-(simdat2$inc-denom.m$fitted.values)*simdat2$hat.h.Y1))
   DR.Eh.Y0.mw<-mean(((1-denom.m$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heapall+(simdat2$inc-denom.m$fitted.values)*simdat2$hat.h.Y0))
   
   #estimate pi under HCAR
   naive.nbin<-glm.nb(CIG.heap ~ 1, data = simdat2)
   init.coefs.MSM<-c(0,1/naive.nbin$theta,coef(naive.nbin))
   HeapNB.IV<-ML.Heap.marginal.HCAR(simdat2,init.coefs.MSM)
   
   #compile into estimator
   CM1.iv.dr.heap.mw<-(HeapNB.IV[1]^(-1))*(DR.EYh1.mw-(1-HeapNB.IV[1])*DR.Eh.Y1.mw)
   CM0.iv.dr.heap.mw<-(HeapNB.IV[1]^(-1))*(DR.EYh0.mw-(1-HeapNB.IV[1])*DR.Eh.Y0.mw)
   
   RR.iv.dr.heap.mw<-CM1.iv.dr.heap.mw/CM0.iv.dr.heap.mw
   
   DR.IV.mw<-c(HeapNB.IV,HeapPoisPG.IV,DR.EYh1.mw,DR.EYh0.mw,DR.Eh.Y1.mw,DR.Eh.Y0.mw,CM1.iv.dr.heap.mw,CM0.iv.dr.heap.mw,RR.iv.dr.heap.mw)
   
   #call estimating function
   DR.ests.heap.HCAR.mw<-DR_est_heap.HCAR(simdat2, estfun_dr_heap.HCAR, inc ~ edu, DR.IV.mw)
   colnames(DR.ests.heap.HCAR.mw) <- c('RR.DR.heap.HCAR.mw', 'seRR.DR.heap.HCAR.mw') 

   ##correct specification of weight model, incorrect specification of outcome model
   
   #compute hat(Y_h^a) and hat(h(Y^a)) for each participant
   out_model.naive.mo  <- glm(CIG.heap ~ edu + inc, data = simdat2, family = poisson)
   HeapPoisPG.IV.mo<-ML.HeapPoisPG.mo.HCAR(simdat2,c(0,coef(out_model.naive.mo)))
   
   simdat2.m$hat.Y0<-(exp(HeapPoisPG.IV.mo[2]+HeapPoisPG.IV.mo[3]*simdat2.m$edu))
   simdat2.m$hat.Y1<-(exp(HeapPoisPG.IV.mo[2]+HeapPoisPG.IV.mo[3]*simdat2.m$edu+HeapPoisPG.IV.mo[4]))
   simdat2.m$hat.h.Y0<-round(simdat2.m$hat.Y0+0.00001,-1)
   simdat2.m$hat.h.Y1<-round(simdat2.m$hat.Y1+0.00001,-1)
   simdat2.m$hat.Yh0<-HeapPoisPG.IV.mo[1]*simdat2.m$hat.Y0+(1-HeapPoisPG.IV.mo[1])*simdat2.m$hat.h.Y0
   simdat2.m$hat.Yh1<-HeapPoisPG.IV.mo[1]*simdat2.m$hat.Y1+(1-HeapPoisPG.IV.mo[1])*simdat2.m$hat.h.Y1 
   
   #compute hat(E(Y_h^a))
   DR.EYh1.mo<-mean((denom$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heap-(simdat2$inc-denom$fitted.values)*simdat2.m$hat.Yh1))
   DR.EYh0.mo<-mean(((1-denom$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heap+(simdat2$inc-denom$fitted.values)*simdat2.m$hat.Yh0))
   
   #compute hat(E(h(Y^a)))
   DR.Eh.Y1.mo<-mean((denom$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heapall-(simdat2$inc-denom$fitted.values)*simdat2.m$hat.h.Y1))
   DR.Eh.Y0.mo<-mean(((1-denom$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heapall+(simdat2$inc-denom$fitted.values)*simdat2.m$hat.h.Y0))
   
   #compile into estimator
   CM1.iv.dr.heap.mo<-(HeapNB.IV[1]^(-1))*(DR.EYh1.mo-(1-HeapNB.IV[1])*DR.Eh.Y1.mo)
   CM0.iv.dr.heap.mo<-(HeapNB.IV[1]^(-1))*(DR.EYh0.mo-(1-HeapNB.IV[1])*DR.Eh.Y0.mo)
   
   RR.iv.dr.heap.mo<-CM1.iv.dr.heap.mo/CM0.iv.dr.heap.mo
   
   DR.IV.mo<-c(HeapNB.IV,HeapPoisPG.IV.mo,DR.EYh1.mo,DR.EYh0.mo,DR.Eh.Y1.mo,DR.Eh.Y0.mo,CM1.iv.dr.heap.mo,CM0.iv.dr.heap.mo,RR.iv.dr.heap.mo)
   
   #call estimating function
   DR.ests.heap.HCAR.mo<-DR_est_heap.HCAR(simdat2, estfun_dr_heap_mis.HCAR, inc ~ income, DR.IV.mo)
   colnames(DR.ests.heap.HCAR.mo) <- c('RR.DR.heap.HCAR.mo', 'seRR.DR.heap.HCAR.mo') 
 
   ##misspecify both models
   
   #compute hat(E(Y_h^a))
   DR.EYh1.mb<-mean((denom.m$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heap-(simdat2$inc-denom.m$fitted.values)*simdat2.m$hat.Yh1))
   DR.EYh0.mb<-mean(((1-denom.m$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heap+(simdat2$inc-denom.m$fitted.values)*simdat2.m$hat.Yh0))
   
   #compute hat(E(h(Y^a)))
   DR.Eh.Y1.mb<-mean((denom.m$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heapall-(simdat2$inc-denom.m$fitted.values)*simdat2.m$hat.h.Y1))
   DR.Eh.Y0.mb<-mean(((1-denom.m$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heapall+(simdat2$inc-denom.m$fitted.values)*simdat2.m$hat.h.Y0))
   
   #compile into estimator
   CM1.iv.dr.heap.mb<-(HeapNB.IV[1]^(-1))*(DR.EYh1.mb-(1-HeapNB.IV[1])*DR.Eh.Y1.mb)
   CM0.iv.dr.heap.mb<-(HeapNB.IV[1]^(-1))*(DR.EYh0.mb-(1-HeapNB.IV[1])*DR.Eh.Y0.mb)
   
   RR.iv.dr.heap.mb<-CM1.iv.dr.heap.mb/CM0.iv.dr.heap.mb
   
   DR.IV.mb<-c(HeapNB.IV,HeapPoisPG.IV.mo,DR.EYh1.mb,DR.EYh0.mb,DR.Eh.Y1.mb,DR.Eh.Y0.mb,CM1.iv.dr.heap.mb,CM0.iv.dr.heap.mb,RR.iv.dr.heap.mb)
   
   #call estimating function
   DR.ests.heap.HCAR.mb<-DR_est_heap.HCAR(simdat2, estfun_dr_heap_mis.HCAR, inc ~ edu, DR.IV.mb)
   colnames(DR.ests.heap.HCAR.mb) <- c('RR.DR.heap.HCAR.mb', 'seRR.DR.heap.HCAR.mb') 
   
  #combine results for three estimators  
  simres<-cbind(IPTW.ests.fixed.heap.HCAR.mw, IPTW.ests.heap.HCAR.mw, PG.ests.heap.HCAR.mo, PG.ests.heap.IH.mo, DR.ests.heap.HCAR.mw, DR.ests.heap.HCAR.mo, DR.ests.heap.HCAR.mb)

  output_filename <- paste(dir_path,"/results_", sim, ".csv", sep="")
  write.csv(simres, output_filename)