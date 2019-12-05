setwd("C:\\Users\\Sara\\Documents\\HHMI_projects\\Molecular_clock\\Analysis")
library(ggplot2)
library(magrittr)
library(tidyverse)
qubitdf=read.csv("qubitcsv.csv")
qubitdf%<>%mutate(FragRep=paste(Fragment, Replicate, sep="_"))
colnames(qubitdf)[1]="Sample"

qubit1300=qubitdf[which(qubitdf$Sample=="VA 1300"),]

qubit_all=qubitdf[which(qubitdf$Sample!="VA 1300"),]



# Grouped
ggplot(qubit1300[1:6,], aes(fill=(FragRep), y=Undiluted.concentration, x=Sample)) + 
  geom_bar(position="dodge", stat="identity") +  
  theme_bw() + ylab("Concentration, ng/ul")

ggplot(qubit1300[7:10,], aes(fill=(FragRep), y=Undiluted.concentration, x=Sample)) + 
  geom_bar(position="dodge", stat="identity") +  
  theme_bw() + ylab("Concentration, ng/ul")



ggplot(qubit_all, aes(fill=(FragRep), y=Undiluted.concentration, x=Sample)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + scale_fill_manual(values=c("#B2EBF2","#00BCD4", "#C8E6C9", "#66BB6A", "#FFF9C4", "#FFEB3B")) + theme_bw() + ylab("Concentration, ng/ul")



ggplot(qubit_all, aes(fill=(FragRep), y=Undiluted.concentration, x=Sample)) + 
geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c("#FFFF00","#FFA500", "#ADFF2F", "#32CD32", "#FF00FF", "#8B008B"))