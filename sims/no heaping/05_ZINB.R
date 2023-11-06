#Program: 05_ZINB.R
#Developed by: BES
#Purpose: Simulate data similar to WIHS as ZINB distribution, analyze using MSM, PG, and DR in geex
#Created: 02.08.20, last updated 02.05.23

library(numDeriv)
library(boot)
library(geex)
library(MASS)
library(pscl)


## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
subdir_name <- as.character(args[1])
num.part <- as.numeric(args[2])
num.sims<-as.numeric(args[3])

user_home_dir <- " " # enter path for top-level directory for the project
dir_path <- paste(user_home_dir,subdir_name,"/",sep="")

# load functions from source files, 
# all of which are in the top-level user_home_dir
source('00_Estimators_11.03.21.R')

sim <- Sys.getenv("SLURM_ARRAY_TASK_ID") 

setwd(dir_path)

log.RR <- 0.5


  skip_to_next <- FALSE
  skip_to_next2 <- FALSE

  #generage age of each woman at start of study
  agevis<-as.vector(runif(num.part, min = 20, max = 40))
  
  #generate participant id
  part.id<-seq(1:num.part)
  
  simdat1<-as.data.frame(cbind(part.id,agevis))
  
  #create a drug use and SE rate for each woman 
  simdat1$drugrate<-inv.logit(0.08-(simdat1$agevis/100)) 
  simdat1$druguse<-rbinom(num.part, 1, simdat1$drugrate) 
  simdat1$sxexrate<-inv.logit(-2.9-simdat1$agevis/100+1.2*simdat1$druguse)
  simdat1$sxex<-rbinom(num.part,1,simdat1$sxexrate) 
  
  #incarceration
  simdat1$incrate<-inv.logit(-1.5+(1-simdat1$agevis/100)+0.5*simdat1$druguse+0.5*simdat1$sxex)
  simdat1$inc<-rbinom(num.part,1,simdat1$incrate)
  
  #create prob of excess zero
  simdat1$exczero0p<-inv.logit(-2.5+(simdat1$agevis/100)-0.3*simdat1$druguse-2*simdat1$sxex)
  simdat1$exczero1p<-inv.logit(-2.5+(simdat1$agevis/100)-0.3*simdat1$druguse-2*simdat1$sxex)
  simdat1$exczero0<-rbinom(num.part,1,simdat1$exczero0p)
  simdat1$exczero1<-rbinom(num.part,1,simdat1$exczero1p)
  
  #create means for susceptible pop
  simdat1$lambda0<-exp(-1-0.005*simdat1$agevis+0.7*simdat1$druguse+3.5*simdat1$sxex)
  simdat1$lambda1<-exp(-1-0.005*simdat1$agevis+0.7*simdat1$druguse+3.5*simdat1$sxex+log.RR)
  
  simdat1$SP0<-ifelse(simdat1$exczero0==1,0,rnegbin(num.part,simdat1$lambda0, theta=2))
  simdat1$SP1<-ifelse(simdat1$exczero1==1,0,rnegbin(num.part,simdat1$lambda1, theta=2))
  
  #assign observed outcome
  simdat1$SP<-ifelse(simdat1$inc==0,simdat1$SP0,simdat1$SP1)
  
  
  #age age/10
  simdat1$agevis10<-simdat1$agevis/10

  ###################################
  ####### MSM with IPTW ############
  ###################################
  
  #Weights
  denom <- glm(inc ~ druguse + agevis10 + sxex, data = simdat1, family = binomial("logit"))
  denom.probs<-denom$fitted.values
  simdat2<-as.data.frame(cbind(simdat1,denom.probs))
  simdat2$incweight<-1/simdat2$denom.probs
  simdat2$noincweight<-1/(1-simdat2$denom.probs)
  simdat2$IPTW<-ifelse(simdat2$inc==1,simdat2$incweight, simdat2$noincweight)

  #Correctly specified weight models
  #calculate initial values for geex
  CM0.iv<-sum(simdat2$IPTW*(1-simdat2$inc)*simdat2$SP)/sum(simdat2$IPTW*(1-simdat2$inc))
  CM1.iv<-sum(simdat2$IPTW*(simdat2$inc)*simdat2$SP)/sum(simdat2$IPTW*(simdat2$inc))
  log.RR.iv<-log(CM1.iv/CM0.iv)
  RR.iv<-CM1.iv/CM0.iv

  #call estimating function, weights known
  MSM.ests.fixed<-MSM_est_fixed(simdat2,simdat2$inc,simdat2$SP,simdat2$IPTW,CM0.iv, CM1.iv,log.RR.iv,RR.iv)

  #call estimating function, weights estimated
  MSM.ests<-MSM_est(simdat2,simdat2$inc,simdat2$SP, inc ~ druguse + agevis10 + sxex, CM0.iv, CM1.iv,log.RR.iv,RR.iv)

    
  #Incorrectly specified weight models
  #calculate initial values for geex
  #Weights
  denom.m <- glm(inc ~ agevis10 + sxex, data = simdat1, family = binomial("logit"))
  denom.probs.m<-denom.m$fitted.values
  simdat2.m<-as.data.frame(cbind(simdat1,denom.probs.m))
  simdat2.m$incweight<-1/simdat2.m$denom.probs.m
  simdat2.m$noincweight<-1/(1-simdat2.m$denom.probs.m)
  simdat2.m$IPTW<-ifelse(simdat2.m$inc==1,simdat2.m$incweight, simdat2.m$noincweight)
  
  CM0.iv.mw<-sum(simdat2.m$IPTW*(1-simdat2.m$inc)*simdat2.m$SP)/sum(simdat2.m$IPTW*(1-simdat2.m$inc))
  CM1.iv.mw<-sum(simdat2.m$IPTW*(simdat2.m$inc)*simdat2.m$SP)/sum(simdat2.m$IPTW*(simdat2.m$inc))
  log.RR.iv.mw<-log(CM1.iv.mw/CM0.iv.mw)
  RR.iv.mw<-CM1.iv.mw/CM0.iv.mw
  
  #call estimating function, weights known
  MSM.ests.fixed.mw<-MSM_est_fixed(simdat2.m,simdat2.m$inc,simdat2.m$SP,simdat2.m$IPTW,CM0.iv.mw, CM1.iv.mw,log.RR.iv.mw,RR.iv.mw)
  colnames(MSM.ests.fixed.mw) <- c('RR.MSM.fixed.mw', 'seRR.MSM.fixed.mw') 
    #data,exposure,outcome,weight,CM0,CM1,logRR,RR  
  
  #call estimating function, weights estimated
  MSM.ests.mw<-MSM_est(simdat2.m,simdat2.m$inc,simdat2.m$SP, inc ~ agevis10 + sxex, CM0.iv.mw, CM1.iv.mw,log.RR.iv.mw,RR.iv.mw)
  colnames(MSM.ests.mw) <- c('RR.MSM.mw', 'seRR.MSM.mw') 
  
  ###################################
  ####### PARAMETRIC G-FORMULA ######
  ###################################
  
  ###Correctly specified outcome model
  #calculate the initial values for geex
  tryCatch(thetamod  <- zeroinfl(SP ~ 1 + inc + druguse + sxex +  agevis10 | 1 + druguse + sxex +  agevis10, data = simdat2, dist = "negbin", link="logit"), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next==FALSE) { 
  init.ZINB<-c(coef(thetamod),log(1/thetamod$theta))
  mod.iv<-ML.ZINB(simdat2,init.ZINB)
  
  CM0.iv.pg<-mean(exp(mod.iv[1]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10)))
  CM1.iv.pg<-mean(exp(mod.iv[1]+mod.iv[2]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10)))
  RR.iv.pg<-CM1.iv.pg/CM0.iv.pg
  
  #call estimating function
    tryCatch(PG.ests<-PG_est(simdat2, estfun_gf_ZINB, mod.iv, CM0.iv.pg, CM1.iv.pg, RR.iv.pg), error = function(e) { skip_to_next2 <<- TRUE})
      
      if(skip_to_next2==TRUE) {
      PG.ests<-matrix(c(NA,NA),nrow=1)
      skip_to_next2<-FALSE}
  
      colnames(PG.ests)<- c('RR.PG', 'seRR.PG') 
  
    
  
  ###################################
  ####### DR Estimation #############
  ###################################
  
  #correct specification of both models
  #calculate the initial values for geex
  pred0.dr<-exp(mod.iv[1]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10))
  pred1.dr<-exp(mod.iv[1]+mod.iv[2]+mod.iv[3]*simdat2$druguse+mod.iv[4]*simdat2$sxex+mod.iv[5]*simdat2$agevis10)*(1-inv.logit(mod.iv[6]+mod.iv[7]*simdat2$druguse+mod.iv[8]*simdat2$sxex+mod.iv[9]*simdat2$agevis10))
  CM1.iv.dr<-mean(((simdat2$inc*simdat2$SP-(simdat2$inc-denom$fitted.values)*pred1.dr))/denom$fitted.values)
  CM0.iv.dr<-mean((((1-simdat2$inc)*simdat2$SP+(simdat2$inc-denom$fitted.values)*pred0.dr))/(1-denom$fitted.values))
  RR.iv.dr<-CM1.iv.dr/CM0.iv.dr
  #call estimating function
  tryCatch(DR.ests<-DR_est(simdat2, estfun_dr_ZINB, inc ~ druguse + agevis10 + sxex, mod.iv, CM0.iv.dr, CM1.iv.dr, RR.iv.dr), error = function(e) { skip_to_next2 <<- TRUE})
  
  if(skip_to_next2==TRUE) {
    DR.ests<-matrix(c(NA,NA),nrow=1)
    skip_to_next2<-FALSE}
  
  colnames(DR.ests)<- c('RR.DR', 'seRR.DR') 
  
  #correct specification of outcome model, incorect specification of weight model
  #calculate the initial values for geex
  CM1.iv.dr.mw<-mean(((simdat2.m$inc*simdat2.m$SP-(simdat2.m$inc-denom.m$fitted.values)*pred1.dr))/denom.m$fitted.values)
  CM0.iv.dr.mw<-mean((((1-simdat2.m$inc)*simdat2.m$SP+(simdat2.m$inc-denom.m$fitted.values)*pred0.dr))/(1-denom.m$fitted.values))
  RR.iv.dr.mw<-CM1.iv.dr.mw/CM0.iv.dr.mw
  
  #call estimating function
  tryCatch(DR.ests.mw<-DR_est(simdat2, estfun_dr_ZINB, inc ~ agevis10 + sxex, mod.iv, CM0.iv.dr.mw, CM1.iv.dr.mw, RR.iv.dr.mw), error = function(e) { skip_to_next2 <<- TRUE})

    if(skip_to_next2==TRUE) {
    DR.ests.mw<-matrix(c(NA,NA),nrow=1)
    skip_to_next2<-FALSE}
  
  colnames(DR.ests.mw)<- c('RR.DR.mw', 'seRR.DR.mw') 
           
  }
  
  if(skip_to_next==TRUE) {
    PG.ests<-matrix(c(NA,NA),nrow=1)
    colnames(PG.ests)<- c('RR.PG', 'seRR.PG') 
    DR.ests<-matrix(c(NA,NA),nrow=1)
    colnames(DR.ests) <- c('RR.DR', 'seRR.DR') 
    DR.ests.mw<-matrix(c(NA,NA),nrow=1)
    colnames(DR.ests.mw) <- c('RR.DR.mw', 'seRR.DR.mw') 
    skip_to_next<-FALSE}
  
  
  ###################################
  ####### PARAMETRIC G-FORMULA ######
  ###################################
  
  ###Incorrectly specified outcome model
  #calculate the initial values for geex
  tryCatch(thetamod.m  <- zeroinfl(SP ~ 1 + inc + sxex +  agevis10 | 1 + sxex +  agevis10, data = simdat2, dist = "negbin", link="logit"), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next==FALSE) { 
  init.ZINB.ms<-c(coef(thetamod.m),log(1/thetamod.m$theta))
  mod.iv.mo<-ML.ZINB.ms(simdat2)
  
  CM0.iv.pg.mo<-mean(exp(mod.iv.mo[1]+mod.iv.mo[3]*simdat2$sxex+mod.iv.mo[4]*simdat2$agevis10)*(1-inv.logit(mod.iv.mo[5]+mod.iv.mo[6]*simdat2$sxex+mod.iv.mo[7]*simdat2$agevis10)))
  CM1.iv.pg.mo<-mean(exp(mod.iv.mo[1]+mod.iv.mo[2]+mod.iv.mo[3]*simdat2$sxex+mod.iv.mo[4]*simdat2$agevis10)*(1-inv.logit(mod.iv.mo[5]+mod.iv.mo[6]*simdat2$sxex+mod.iv.mo[7]*simdat2$agevis10)))
  RR.iv.pg.mo<-CM1.iv.pg.mo/CM0.iv.pg.mo
  
  #call estimating function
  tryCatch(PG.ests.mo<-PG_est(simdat2, estfun_gf_inc_ZINB, mod.iv.mo, CM0.iv.pg.mo, CM1.iv.pg.mo, RR.iv.pg.mo), error = function(e) { skip_to_next2 <<- TRUE})
  
  if(skip_to_next2==TRUE) {
    PG.ests.mo<-matrix(c(NA,NA),nrow=1)
    skip_to_next2<-FALSE}
  
  colnames(PG.ests.mo)<- c('RR.PG.mo', 'seRR.PG.mo') 

   
  ###################################
  ####### DR Estimation #############
  ###################################
  
  #correct specification of weight model, incorect specification of outcome model
  #calculate the initial values for geex
  pred0.dr.mo<-exp(mod.iv.mo[1]+mod.iv.mo[3]*simdat2$sxex+mod.iv.mo[4]*simdat2$agevis10)*(1-inv.logit(mod.iv.mo[5]+mod.iv.mo[6]*simdat2$sxex+mod.iv.mo[7]*simdat2$agevis10))
  pred1.dr.mo<-exp(mod.iv.mo[1]+mod.iv.mo[2]+mod.iv.mo[3]*simdat2$sxex+mod.iv.mo[4]*simdat2$agevis10)*(1-inv.logit(mod.iv.mo[5]+mod.iv.mo[6]*simdat2$sxex+mod.iv.mo[7]*simdat2$agevis10))
  CM1.iv.dr.mo<-mean(((simdat2$inc*simdat2$SP-(simdat2$inc-denom$fitted.values)*pred1.dr.mo))/denom$fitted.values)
  CM0.iv.dr.mo<-mean((((1-simdat2$inc)*simdat2$SP+(simdat2$inc-denom$fitted.values)*pred0.dr.mo))/(1-denom$fitted.values))
  RR.iv.dr.mo<-CM1.iv.dr.mo/CM0.iv.dr.mo
  #call estimating function
  tryCatch(DR.ests.mo<-DR_est(simdat2, estfun_dr_inc_ZINB, inc ~ druguse + agevis10 + sxex, mod.iv.mo, CM0.iv.dr.mo, CM1.iv.dr.mo, RR.iv.dr.mo), error = function(e) { skip_to_next2 <<- TRUE})

  if(skip_to_next2==TRUE) {
    DR.ests.mo<-matrix(c(NA,NA),nrow=1)
    skip_to_next2<-FALSE}
  
  colnames(DR.ests.mo)<- c('RR.DR.mo', 'seRR.DR.mo')   
  
  #incorrect specification of both models
  CM1.iv.dr.mb<-mean(((simdat2.m$inc*simdat2.m$SP-(simdat2.m$inc-denom.m$fitted.values)*pred1.dr.mo))/denom.m$fitted.values)
  CM0.iv.dr.mb<-mean((((1-simdat2.m$inc)*simdat2.m$SP+(simdat2.m$inc-denom.m$fitted.values)*pred0.dr.mo))/(1-denom.m$fitted.values))
  RR.iv.dr.mb<-CM1.iv.dr.mb/CM0.iv.dr.mb
  #call estimating function
  tryCatch(DR.ests.mb<-DR_est(simdat2, estfun_dr_inc_ZINB, inc ~ agevis10 + sxex, mod.iv.mo, CM0.iv.dr.mb, CM1.iv.dr.mb, RR.iv.dr.mb), error = function(e) { skip_to_next2 <<- TRUE})

  if(skip_to_next2==TRUE) {
    DR.ests.mb<-matrix(c(NA,NA),nrow=1)
    skip_to_next2<-FALSE}
  
  colnames(DR.ests.mb)<- c('RR.DR.mb', 'seRR.DR.mb') 
           
    }
  
  if(skip_to_next==TRUE) {
    PG.ests.mo<-matrix(c(NA,NA),nrow=1)
    colnames(PG.ests.mo) <- c('RR.PG.mo', 'seRR.PG.mo') 
    DR.ests.mo<-matrix(c(NA,NA),nrow=1)
    colnames(DR.ests.mo) <- c('RR.DR.mo', 'seRR.DR.mo') 
    DR.ests.mb<-matrix(c(NA,NA),nrow=1)
    colnames(DR.ests.mb) <- c('RR.DR.mb', 'seRR.DR.mb') 
    skip_to_next<-FALSE}
  
 #combine results for three estimators
  simres<-cbind(MSM.ests.fixed, MSM.ests,MSM.ests.fixed.mw, MSM.ests.mw, PG.ests, PG.ests.mo, DR.ests, DR.ests.mw, DR.ests.mo, DR.ests.mb)

  output_filename <- paste(dir_path,"/results_", sim, ".csv", sep="")
  write.csv(simres, output_filename)