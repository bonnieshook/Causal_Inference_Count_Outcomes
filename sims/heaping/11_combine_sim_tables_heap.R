#Program: Combine heaping sim results and output tables
#Developed by: BES
#Created: 12.03.19, last updated 07.18.23

log.RR<-0.25
nsim<-5000

#bring in results
CS.S1<-as.data.frame(read.csv("results_HeapCorr800_HCAR.csv"))
IS.S1<-as.data.frame(read.csv("results_HeapInc800_HCAR.csv"))


#calculate emp bias and coverage
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

##################################
# SCENARIO 1 (HCAR)
##################################
#correct specification 
IPTW.naive<-calc_table(CS.S1$RR.IPTW.naive,CS.S1$seRR.IPTW.naive)
IPTW.naive.fixed<-calc_table(CS.S1$RR.IPTW.fixed.naive,CS.S1$seRR.IPTW.fixed.naive)
IPTW.HCAR<-calc_table(CS.S1$RR.IPTW.heap.HCAR,CS.S1$seRR.IPTW.heap.HCAR)
IPTW.HCAR.fixed<-calc_table(CS.S1$RR.IPTW.fixed.heap.HCAR,CS.S1$seRR.IPTW.fixed.heap.HCAR)

PG.naive<-calc_table(CS.S1$RR.PG.naive,CS.S1$seRR.PG.naive)
PG.HCAR<-calc_table(CS.S1$RR.PG.heap.HCAR,CS.S1$seRR.PG.heap.HCAR)
PG.IH<-calc_table(CS.S1$RR.PG.heap.IH,CS.S1$seRR.PG.heap.IH)

DR.naive<-calc_table(CS.S1$RR.DR.naive,CS.S1$seRR.DR.naive)
DR.HCAR<-calc_table(CS.S1$RR.DR.heap.HCAR,CS.S1$seRR.DR.heap.HCAR)

### incorrect specification
IPTW.fixed.heap.HCAR.mw<-calc_table(IS.S1$RR.IPTW.fixed.heap.HCAR.mw,IS.S1$seRR.IPTW.fixed.heap.HCAR.mw)
IPTW.heap.HCAR.mw<-calc_table(IS.S1$RR.IPTW.heap.HCAR.mw,IS.S1$seRR.IPTW.heap.HCAR.mw)
PG.Heap.HCAR.mo<-calc_table(IS.S1$RR.PG.heap.HCAR.mo,IS.S1$seRR.PG.heap.HCAR.mo)
PG.Heap.IH.mo<-calc_table(IS.S1$RR.PG.heap.IH.mo,IS.S1$seRR.PG.heap.IH.mo)
DR.Heap.HCAR.mw<-calc_table(IS.S1$RR.DR.heap.HCAR.mw,IS.S1$seRR.DR.heap.HCAR.mw)
DR.Heap.HCAR.mo<-calc_table(IS.S1$RR.DR.heap.HCAR.mo,IS.S1$seRR.DR.heap.HCAR.mo)
DR.Heap.HCAR.mb<-calc_table(IS.S1$RR.DR.heap.HCAR.mb,IS.S1$seRR.DR.heap.HCAR.mb)

Scenario1<-rbind(IPTW.naive.fixed, IPTW.HCAR.fixed,IPTW.fixed.heap.HCAR.mw, IPTW.naive,IPTW.HCAR,IPTW.heap.HCAR.mw,PG.naive,PG.HCAR,PG.Heap.HCAR.mo,PG.IH,PG.Heap.IH.mo,DR.naive,DR.HCAR,DR.Heap.HCAR.mw,DR.Heap.HCAR.mo,DR.Heap.HCAR.mb)

write.csv(Scenario1, "Heapsims_summary_S1_n800.csv") 

rm(list= ls()[!(ls() %in% c('nsim','log.RR','calc_table'))])


##################################
# SCENARIO 2 (IH)
##################################

CS.S2<-as.data.frame(read.csv("results_HeapCorr800_IH.csv"))
IS.S2<-as.data.frame(read.csv("results_HeapInc800_IH.csv"))

### correct specification 
IPTW.naive<-calc_table(CS.S2$RR.IPTW.naive,CS.S2$seRR.IPTW.naive)
IPTW.naive.fixed<-calc_table(CS.S2$RR.IPTW.fixed.naive,CS.S2$seRR.IPTW.fixed.naive)
IPTW.HCAR<-calc_table(CS.S2$RR.IPTW.heap.HCAR,CS.S2$seRR.IPTW.heap.HCAR)
IPTW.HCAR.fixed<-calc_table(CS.S2$RR.IPTW.fixed.heap.HCAR,CS.S2$seRR.IPTW.fixed.heap.HCAR)

PG.naive<-calc_table(CS.S2$RR.PG.naive,CS.S2$seRR.PG.naive)
PG.HCAR<-calc_table(CS.S2$RR.PG.heap.HCAR,CS.S2$seRR.PG.heap.HCAR)
PG.IH<-calc_table(CS.S2$RR.PG.heap.IH,CS.S2$seRR.PG.heap.IH)

DR.naive<-calc_table(CS.S2$RR.DR.naive,CS.S2$seRR.DR.naive)
DR.HCAR<-calc_table(CS.S2$RR.DR.heap.HCAR,CS.S2$seRR.DR.heap.HCAR)

### incorrect specification
IPTW.fixed.heap.HCAR.mw<-calc_table(IS.S2$RR.IPTW.fixed.heap.HCAR.mw,IS.S2$seRR.IPTW.fixed.heap.HCAR.mw)
IPTW.heap.HCAR.mw<-calc_table(IS.S2$RR.IPTW.heap.HCAR.mw,IS.S2$seRR.IPTW.heap.HCAR.mw)
PG.Heap.HCAR.mo<-calc_table(IS.S2$RR.PG.heap.HCAR.mo,IS.S2$seRR.PG.heap.HCAR.mo)
PG.Heap.IH.mo<-calc_table(IS.S2$RR.PG.heap.IH.mo,IS.S2$seRR.PG.heap.IH.mo)
DR.Heap.HCAR.mw<-calc_table(IS.S2$RR.DR.heap.HCAR.mw,IS.S2$seRR.DR.heap.HCAR.mw)
DR.Heap.HCAR.mo<-calc_table(IS.S2$RR.DR.heap.HCAR.mo,IS.S2$seRR.DR.heap.HCAR.mo)
DR.Heap.HCAR.mb<-calc_table(IS.S2$RR.DR.heap.HCAR.mb,IS.S2$seRR.DR.heap.HCAR.mb)

Scenario2<-rbind(IPTW.naive.fixed, IPTW.HCAR.fixed,IPTW.fixed.heap.HCAR.mw, IPTW.naive,IPTW.HCAR,IPTW.heap.HCAR.mw,PG.naive,PG.HCAR,PG.Heap.HCAR.mo,PG.IH,PG.Heap.IH.mo,DR.naive,DR.HCAR,DR.Heap.HCAR.mw,DR.Heap.HCAR.mo,DR.Heap.HCAR.mb)
write.csv(Scenario2, "Heapsims_summary_S2_n800.csv")

rm(list= ls()[!(ls() %in% c('nsim','log.RR'))])



