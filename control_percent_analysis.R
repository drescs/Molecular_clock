library(magrittr)

setwd("O:/JOLabShared/SaraDrescher/Molecular_clock/nano")
control_names=read.csv("O:/JOLabShared/SaraDrescher/Molecular_clock/nano/control_sample_names.csv", header=FALSE)

#we'll blast all these reads against LAI and Q23, here in our "Reference database":
#setwd O:/JOLabShared/SaraDrescher/Molecular_clock/LAI_Q23_DB
#makeblastdb -in K03455_and_Q23.fasta -dbtype nucl -parse_seqids -out LAI_Q23_DB

#datapath="/fh/fast/bloom_j/SR/ngs/illumina/sdresche/191014_M04866_0296_000000000-CM883/Data/Intensities/BaseCalls"
datapath="O:/JOLabShared/SaraDrescher/Molecular_clock/nano"

#write the commands for use in terminal
commands=vector()
outputfiles=vector()
for (i in 1:length(control_names[,1])){
  sample_num=control_names[i, 1]
  sample_name=control_names[i,2]
  file1=paste0(datapath, "/", sample_num, "_S", sample_num, "_L001_R1_001.fastq.gz")
  file2=paste0(datapath, "/", sample_num, "_S", sample_num, "_L001_R2_001.fastq.gz")
  command1=paste0("pandaseq -f ", file1," -r ", file2," -w ", sample_name, "_merged2.fasta")
  command2=paste0("blastn -db /home/sdresche/JOLabShared/SaraDrescher/Molecular_clock/LAI_Q23_DB/LAI_Q23_DB -query ", sample_name, "_merged.fasta -out ",sample_name,"_blastout2.txt -outfmt 6 -max_target_seqs 1")
  output=paste0(sample_name,"_blastout2.txt")
  outputfiles=c(outputfiles,output)
  commands=c(commands, command1, command2)
}

writeLines(c(noquote("#!/bin/sh"),commands), con="O:/JOLabShared/SaraDrescher/Molecular_clock/second_run/control_blast_commands.txt")

#write_lines(commands, "blastcommands.txt")
#run this in the terminal

control_names=read.csv("O:/JOLabShared/SaraDrescher/Molecular_clock/nano/control_sample_names.csv", header=FALSE)

#count how many sequences aligned to LAI (HXB2) or Q23
HXB2counts=vector()
Q23counts=vector()
for (i in 1:length(outputfiles)){
  blast=read.delim(outputfiles[i], sep="\t", header=FALSE)
  blast_unique=unique(blast[,1:2])
  HXB2=length(which(blast_unique[,2]=="K03455|HIVHXB2CG"))
  Q23=length(which(blast_unique[,2]=="AF004885.1"))
  HXB2counts=c(HXB2counts, HXB2)
  Q23counts=c(Q23counts, Q23)
}

percent_Q23=Q23counts/(HXB2counts+Q23counts)*100
round_percent_Q23=round(percent_Q23, digits=3)
control_percent=cbind(control_names, HXB2counts, Q23counts, round_percent_Q23)
control_percent %<>% separate(V2, c("sample", "frag"), sep="_", remove=FALSE, fill="left")


control_percent=control_percent[order(control_percent$V3),]
#corr_percent=select(control_percent, rep, frag_sample, round_percent_Q23)
#corr_percent=spread(corr_percent, rep, round_percent_Q23)
#ggplot(corr_percent, aes(x=R1, y=R2)) + geom_point() + geom_smooth(method=lm, se=FALSE)
control_percent%<>%rename("V3"="Expected_Percent")
control_percent$Expected_Percent=as.numeric(control_percent$Expected_Percent)
control_percent=control_percent[order(control_percent$Expected_Percent),]
control_percent$Percent_Expected=factor(control_percent$Expected_Percent, levels = c('0','1','5','10','100'))
ggplot(control_percent, aes(fill=Percent_Expected, x=(reorder(V2, Expected_Percent)), y=round_percent_Q23)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("First Run") +
 theme_classic() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
  
  labs(y="Percent Q23")



write.csv(control_percent, "Control percent 2nd.csv")
control_perect=read.csv("O:/JOLabShared/SaraDrescher/Molecular_clock/nano/Control percent.csv")

