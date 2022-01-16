#Program: Combine Sim results and output table
#Developed by: BES
#Created: 12.03.19, last updated 01.09.22

log.RR<-0.25

#bring in results
CS<-as.data.frame(read.csv("HeapSims_cs_010922_n800.csv"))
IS<-as.data.frame(read.csv("HeapSims_IS_010922_n800.csv"))

#for IS, make estimate missing if SE is missing
IS$RR.DR.heap.mw2<-ifelse(IS$seRR.DR.heap.mw=="NA","NA",IS$RR.DR.heap.mw)
IS$RR.DR.heap.mo2<-ifelse(IS$seRR.DR.heap.mo=="NA","NA",IS$RR.DR.heap.mo)
IS$RR.DR.heap.mb2<-ifelse(IS$seRR.DR.heap.mb=="NA","NA",IS$RR.DR.heap.mb)

#count missings
pct.missing.DR.MW<-100*(sum(is.na(IS$RR.DR.heap.mw2))/2000)
pct.missing.DR.MO<-100*(sum(is.na(IS$RR.DR.heap.mo2))/2000)
pct.missing.DR.MB<-100*(sum(is.na(IS$RR.DR.heap.mb2))/2000)

pct.missing.DR.MW
pct.missing.DR.MO
pct.missing.DR.MB

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

#correct specification
MSM.naive<-calc_table(CS$RR.MSM.naive,CS$seRR.MSM.naive)
MSM.naive.fixed<-calc_table(CS$RR.MSM.fixed.naive,CS$seRR.MSM.fixed.naive)
MSM.heap<-calc_table(CS$RR.MSM.heap,CS$seRR.MSM.heap)
MSM.heap.fixed<-calc_table(CS$RR.MSM.fixed.heap,CS$seRR.MSM.fixed.heap)

PG.naive<-calc_table(CS$RR.PG.naive,CS$seRR.PG.naive)
PG.Heap<-calc_table(CS$RR.PG.heap,CS$seRR.PG.heap)

DR.naive<-calc_table(CS$RR.DR.naive,CS$seRR.DR.naive)
DR.Heap<-calc_table(CS$RR.DR.heap,CS$seRR.DR.heap)


combine.cs<-rbind(MSM.naive.fixed, MSM.heap.fixed, MSM.naive,MSM.heap,PG.naive,PG.Heap,DR.naive,DR.Heap)

write.csv(combine.cs, "Heapsims_summary_CS_010922_n800.csv") #Web Table 5


### incorrect specification
MSM.fixed.heap.mw<-calc_table(IS$RR.MSM.fixed.heap.mw,IS$seRR.MSM.fixed.heap.mw)
MSM.heap.mw<-calc_table(IS$RR.MSM.heap.mw,IS$seRR.MSM.heap.mw)
PG.Heap.mo<-calc_table(IS$RR.PG.heap.mo,IS$seRR.PG.heap.mo)
DR.Heap.mw<-calc_table(IS$RR.DR.heap.mw2,IS$seRR.DR.heap.mw)
DR.Heap.mo<-calc_table(IS$RR.DR.heap.mo2,IS$seRR.DR.heap.mo)
DR.Heap.mb<-calc_table(IS$RR.DR.heap.mb2,IS$seRR.DR.heap.mb)


combine.is<-rbind(MSM.fixed.heap.mw,MSM.heap.mw,DR.Heap.mw,PG.Heap.mo,DR.Heap.mo,DR.Heap.mb)

write.csv(combine.is, "Heapsims_summary_IS_010922_n800.csv") #Web Table 6
