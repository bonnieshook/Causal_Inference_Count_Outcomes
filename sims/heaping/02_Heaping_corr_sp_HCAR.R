#Program: 02_Heaping_corr_sp_HCAR.R
#Developed by: BES
#Purpose: Simulate data similar to WIHS, analyze using MSM, PG, and DR in geex

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

   #calculate initial values for geex
   CM0.iv.naive<-sum(simdat2$IPTW*(1-simdat2$inc)*simdat2$CIG.heap)/sum(simdat2$IPTW*(1-simdat2$inc))
   CM1.iv.naive<-sum(simdat2$IPTW*(simdat2$inc)*simdat2$CIG.heap)/sum(simdat2$IPTW*(simdat2$inc))
   log.RR.iv.naive<-log(CM1.iv.naive/CM0.iv.naive)
   RR.iv.naive<-CM1.iv.naive/CM0.iv.naive

   #call estimating function, weights known
   IPTW.ests.fixed.naive<-MSM_est_fixed(simdat2,simdat2$inc,simdat2$CIG.heap,simdat2$IPTW,CM0.iv.naive, CM1.iv.naive,log.RR.iv.naive,RR.iv.naive)
   colnames(IPTW.ests.fixed.naive) <- c('RR.IPTW.fixed.naive', 'seRR.IPTW.fixed.naive') 
   
   #call estimating function, weights estimated
   IPTW.ests.naive<-MSM_est(simdat2,simdat2$inc,simdat2$CIG.heap, inc ~ income, CM0.iv.naive, CM1.iv.naive,log.RR.iv.naive,RR.iv.naive)
   colnames(IPTW.ests.naive) <- c('RR.IPTW.naive', 'seRR.IPTW.naive') 
    
   ###########HEAPING#############

   ### Correct specification, HCAR
   
   ###Estimate E(h(Y^a)) (Note: E(Y_h^a) are CM0.iv.naive and CM1.iv.naive)
   CM0.hY<-sum(simdat2$IPTW*(1-simdat2$inc)*simdat2$CIG.heapall)/sum(simdat2$IPTW*(1-simdat2$inc))
   CM1.hY<-sum(simdat2$IPTW*(simdat2$inc)*simdat2$CIG.heapall)/sum(simdat2$IPTW*(simdat2$inc))
   
   #estimate pi under HCAR
   naive.nbin<-glm.nb(CIG.heap ~ 1, data = simdat2)
   init.coefs.MSM<-c(0,1/naive.nbin$theta,coef(naive.nbin))
   HeapNB.IV<-ML.Heap.marginal.HCAR(simdat2,init.coefs.MSM)
   
   #calculate IPTW estimate of CMR
   CM0.IPTW.heap.HCAR<-((HeapNB.IV[1])^(-1))*(CM0.iv.naive-(1-(HeapNB.IV[1]))*CM0.hY)
   CM1.IPTW.heap.HCAR<-((HeapNB.IV[1])^(-1))*(CM1.iv.naive-(1-(HeapNB.IV[1]))*CM1.hY)
   CMR.IPTW.HCAR<-CM1.IPTW.heap.HCAR/CM0.IPTW.heap.HCAR
   
   Heap.IPTW.IV.HCAR<-c(HeapNB.IV,CM1.iv.naive,CM0.iv.naive,CM1.hY,CM0.hY,CM1.IPTW.heap.HCAR,CM0.IPTW.heap.HCAR,CMR.IPTW.HCAR)
   
   IPTW.ests.fixed.heap.HCAR<-IPTW_est_fixed_heap_NB.HCAR(simdat2,simdat2$inc,simdat2$CIG.heap,simdat2$IPTW,Heap.IPTW.IV.HCAR)
   colnames(IPTW.ests.fixed.heap.HCAR) <- c('RR.IPTW.fixed.heap.HCAR', 'seRR.IPTW.fixed.heap.HCAR') 
 
   ### weights estimated
   IPTW.ests.heap.HCAR<-IPTW_est_heap_NB.HCAR(simdat2,simdat2$inc,simdat2$CIG.heap, inc ~ income, Heap.IPTW.IV.HCAR)
   colnames(IPTW.ests.heap.HCAR) <- c('RR.IPTW.heap.HCAR', 'seRR.IPTW.heap.HCAR') 
   

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
   ####Correct specification - HCAR
   HeapPoisPG.IV<-ML.HeapPoisPG.HCAR(simdat2,c(0,coef(out_model.naive)))

   CM0.iv.pg.heap<-mean(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income))
   CM1.iv.pg.heap<-mean(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income+HeapPoisPG.IV[4]))
   RR.iv.pg.heap<-CM1.iv.pg.heap/CM0.iv.pg.heap
   
   Heap.PG.IV.HCAR<-c(HeapPoisPG.IV,CM1.iv.pg.heap,CM0.iv.pg.heap,RR.iv.pg.heap)
  
   PG.ests.heap.HCAR<-PG_est_heap.HCAR(simdat2,estfun_gf_Pois_heap.HCAR, Heap.PG.IV.HCAR)
   colnames(PG.ests.heap.HCAR) <- c('RR.PG.heap.HCAR', 'seRR.PG.heap.HCAR') 
   
   ####Correct specification - Informative Heaping
   HeapPoisPG.IV.IH<-ML.HeapPoisPG.IH(simdat2,c(0,0,coef(out_model.naive)))
   
   CM0.iv.pg.heap.IH<-mean(exp(HeapPoisPG.IV.IH[3]+ HeapPoisPG.IV.IH[4]*simdat2$income))
   CM1.iv.pg.heap.IH<-mean(exp(HeapPoisPG.IV.IH[3]+ HeapPoisPG.IV.IH[4]*simdat2$income+ HeapPoisPG.IV.IH[5]))
   RR.iv.pg.heap.IH<-CM1.iv.pg.heap.IH/CM0.iv.pg.heap.IH
   
   Heap.PG.IV.IH<-c(HeapPoisPG.IV.IH,CM1.iv.pg.heap.IH,CM0.iv.pg.heap.IH,RR.iv.pg.heap.IH)

   PG.ests.heap.IH<-PG_est_heap.IH(simdat2,estfun_gf_Pois_heap.IH, Heap.PG.IV.IH)
   colnames(PG.ests.heap.IH) <- c('RR.PG.heap.IH', 'seRR.PG.heap.IH') 
      
   
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
   
   #compute hat(Y_h^a) and hat(h(Y^a)) for each participant
   simdat2$hat.Y0<-(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income))
   simdat2$hat.Y1<-(exp(HeapPoisPG.IV[2]+HeapPoisPG.IV[3]*simdat2$income+HeapPoisPG.IV[4]))
   simdat2$hat.h.Y0<-round(simdat2$hat.Y0+0.00001,-1)
   simdat2$hat.h.Y1<-round(simdat2$hat.Y1+0.00001,-1)
   simdat2$hat.Yh0<- Heap.PG.IV.HCAR[1]*simdat2$hat.Y0+(1- Heap.PG.IV.HCAR[1])*simdat2$hat.h.Y0
   simdat2$hat.Yh1<- Heap.PG.IV.HCAR[1]*simdat2$hat.Y1+(1- Heap.PG.IV.HCAR[1])*simdat2$hat.h.Y1 
   
   #compute hat(E(Y_h^a))
   DR.EYh1<-mean((denom$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heap-(simdat2$inc-denom$fitted.values)*simdat2$hat.Yh1))
   DR.EYh0<-mean(((1-denom$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heap+(simdat2$inc-denom$fitted.values)*simdat2$hat.Yh0))

   #compute hat(E(h(Y^a)))
   DR.Eh.Y1<-mean((denom$fitted.values^(-1))*(simdat2$inc*simdat2$CIG.heapall-(simdat2$inc-denom$fitted.values)*simdat2$hat.h.Y1))
   DR.Eh.Y0<-mean(((1-denom$fitted.values)^(-1))*((1-simdat2$inc)*simdat2$CIG.heapall+(simdat2$inc-denom$fitted.values)*simdat2$hat.h.Y0))
   
   #compile into estimator
   CM1.iv.dr.heap<-(Heap.IPTW.IV.HCAR[1]^(-1))*(DR.EYh1-(1-Heap.IPTW.IV.HCAR[1])*DR.Eh.Y1)
   CM0.iv.dr.heap<-(Heap.IPTW.IV.HCAR[1]^(-1))*(DR.EYh0-(1-Heap.IPTW.IV.HCAR[1])*DR.Eh.Y0)
   
   RR.iv.dr.heap<-CM1.iv.dr.heap/CM0.iv.dr.heap
   
   DR.IV.HCAR<-c(HeapNB.IV,HeapPoisPG.IV,DR.EYh1,DR.EYh0,DR.Eh.Y1,DR.Eh.Y0,CM1.iv.dr.heap,CM0.iv.dr.heap,RR.iv.dr.heap)
      
   #call estimating function
   DR.ests.heap.HCAR<-DR_est_heap.HCAR(simdat2, estfun_dr_heap.HCAR, inc ~ income, DR.IV.HCAR)
   colnames(DR.ests.heap.HCAR) <- c('RR.DR.heap.HCAR', 'seRR.DR.heap.HCAR') 

 
  #combine results for three estimators
  simres<-cbind(IPTW.ests.fixed.naive, IPTW.ests.naive, IPTW.ests.fixed.heap.HCAR, IPTW.ests.heap.HCAR, PG.ests.naive, PG.ests.heap.HCAR,PG.ests.heap.IH, DR.ests.naive,DR.ests.heap.HCAR)


  output_filename <- paste(dir_path,"/results_", sim, ".csv", sep="")
  write.csv(simres, output_filename)