#Program: Combine Sim results and output table
#Developed by: BES
#Created: 12.03.19, last updated 02.07.23

log.RR<-0.5
nsim<-5000

#bring in results
Pois<-as.data.frame(read.csv("results_Poisson800.csv"))
Nbin<-as.data.frame(read.csv("results_Nbin800.csv"))
ZIP<-as.data.frame(read.csv("results_ZIP800.csv"))
ZINB<-as.data.frame(read.csv("results_ZINB800.csv"))

#calculate bias
calc_table <- function(RR,seRR){
  RRbias<-mean(RR, na.rm = TRUE)-exp(log.RR)
  RRpctbias<-(RRbias/exp(log.RR))*100
  ll <- na.omit(RR)-1.96*na.omit(seRR)
  ul <- na.omit(RR)+1.96*na.omit(seRR)
  cov <-100*mean((ll <= exp(log.RR) & ul >= exp(log.RR)))
  MSE<-median(na.omit(seRR))
  ESE<-sd(na.omit(RR))
  SER<-MSE/ESE
  results<-c(round(RRpctbias,digits=1),"&", round(MSE, digits=2),"&", round(ESE, digits=2),"&", round(SER, digits=2),"&", round(cov, digits=0),'\\')
}

calc_table_bcor <- function(Dist,label){
Dist.MSM.fixed.cor<-c(label, "&", "IPTW, fixed", "&", calc_table(Dist$RR.MSM.fixed,Dist$seRR.MSM.fixed))
Dist.MSM.est.cor<-c("", "&", "IPTW, estimated", "&", calc_table(Dist$RR.MSM,Dist$seRR.MSM))
Dist.PG.cor<-c("", "&", "PG", "&", calc_table(Dist$RR.PG,Dist$seRR.PG))
Dist.DR.bcor<-c("", "&", "DR", "&", calc_table(Dist$RR.DR,Dist$seRR.DR))
all.cor<-rbind(Dist.MSM.fixed.cor,Dist.MSM.est.cor,Dist.PG.cor,Dist.DR.bcor)
}

#tables
Pois.cor<-calc_table_bcor(Pois, "Poisson")
NBin.cor<-calc_table_bcor(Nbin,"NB")
ZIP.cor<-calc_table_bcor(ZIP, "ZIP")
ZINB.cor<-calc_table_bcor(ZINB, "ZINB")

bcorrect<-rbind(Pois.cor,NBin.cor,ZIP.cor,ZINB.cor)

###misspecified weights
calc_table_mw <- function(Dist,label){
  Dist.MSM.fixed.mw<-c(label, "&", "IPTW, fixed, MW", "&", calc_table(Dist$RR.MSM.fixed.mw,Dist$seRR.MSM.fixed.mw))
  Dist.MSM.est.mw<-c("", "&", "IPTW, estimated, MW", "&", calc_table(Dist$RR.MSM.mw,Dist$seRR.MSM.mw))
  Dist.DR.mw<-c("", "&", "DR, MW", "&", calc_table(Dist$RR.DR.mw,Dist$seRR.DR.mw))
  mw<-rbind(Dist.MSM.fixed.mw,Dist.MSM.est.mw,Dist.DR.mw)
}

Pois.mw<-calc_table_mw(Pois, "Poisson")
NBin.mw<-calc_table_mw(Nbin,"NB")
ZIP.mw<-calc_table_mw(ZIP, "ZIP")
ZINB.mw<-calc_table_mw(ZINB, "ZINB")

mw<-rbind(Pois.mw,NBin.mw,ZIP.mw,ZINB.mw)

###misspecified outcome
calc_table_mo <- function(Dist,label){
  Dist.PG.mo<-c(label, "&", "PG, MO", "&", calc_table(Dist$RR.PG.mo,Dist$seRR.PG.mo))
  Dist.DR.mo<-c("", "&", "DR, MO", "&", calc_table(Dist$RR.DR.mo,Dist$seRR.DR.mo))
  mo<-rbind(Dist.PG.mo,Dist.DR.mo)
}

Pois.mo<-calc_table_mo(Pois, "Poisson")
NBin.mo<-calc_table_mo(Nbin,"NB")
ZIP.mo<-calc_table_mo(ZIP, "ZIP")
ZINB.mo<-calc_table_mo(ZINB, "ZINB")

mo<-rbind(Pois.mo,NBin.mo,ZIP.mo,ZINB.mo)

###misspecified both
calc_table_mb <- function(Dist,label){
  Dist.DR.mb<-c(label, "&", "DR, MB", "&", calc_table(Dist$RR.DR.mb,Dist$seRR.DR.mb))
}

Pois.mb<-calc_table_mb(Pois, "Poisson")
NBin.mb<-calc_table_mb(Nbin,"NB")
ZIP.mb<-calc_table_mb(ZIP, "ZIP")
ZINB.mb<-calc_table_mb(ZINB, "ZINB")

mb<-rbind(Pois.mb,NBin.mb,ZIP.mb,ZINB.mb)

misspec.all<-rbind(mw,mo,mb)

write.csv(bcorrect, "SimSummaryBCorrect.csv")
write.csv(misspec.all, "SimSummaryMisspec.csv")




