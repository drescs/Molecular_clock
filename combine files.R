sample_names=read.csv("O:/Documents/molecular_clock/sample_names.csv", stringsAsFactors = FALSE)

sample_names%<>%mutate(fragment=unlist(strsplit(Sample_Name, split="_"))[grepl("F.", unlist(strsplit(Sample_Name, split="_")))])
sample_names%<>%mutate(replicate=unlist(strsplit(Sample_Name, split="_"))[grepl("R.", unlist(strsplit(Sample_Name, split="_")))])
sample_names%<>% mutate(Name=unlist(strsplit(sample_names$Sample_Name, split="_F._R.")))

pt_names=sample_names[grepl("M", sample_names$Name),]
pt_names%<>%mutate(Month=unlist(strsplit(Sample_Name, split="_"))[grepl("M.", unlist(strsplit(Sample_Name, split="_")))])



pt_names_R1=pt_names[which(pt_names$replicate=="R1"),]
pt_names_R2=pt_names[which(pt_names$replicate=="R2"),]

pt_names_R1_F1=pt_names_R1[which(pt_names_R1$fragment=="F1"),]
pt_names_R2_F1=pt_names_R2[which(pt_names_R2$fragment=="F1"),]
pt_names_R1_F2=pt_names_R1[which(pt_names_R1$fragment=="F2"),]
pt_names_R2_F2=pt_names_R2[which(pt_names_R2$fragment=="F2"),]

F1andF2_R1=full_join(pt_names_R1_F1, pt_names_R1_F2, by="Name")
F1andF2_R2=full_join(pt_names_R2_F1, pt_names_R2_F2, by="Name")

F1andF2_both=rbind(F1andF2_R1, F1andF2_R2)

combine_files=function(F1andF2_both){
  line0="#!/bin/sh"
  mytext=noquote(line0)
  myfiles=vector()
  for (i in 1:length(F1andF2_both$Name)){
    samp1=F1andF2_both$Sample_ID.x[i]
    samp2=F1andF2_both$Sample_ID.y[i]
    expname="nano"
    file1=paste0(expname, "/", samp1, "_S", samp1, "_L001_R1_001.fastq.gz")
    file2=paste0(expname, "/", samp1, "_S", samp1, "_L001_R2_001.fastq.gz")
    file3=paste0(expname, "/", samp2, "_S", samp2, "_L001_R1_001.fastq.gz")
    file4=paste0(expname, "/", samp2, "_S", samp2, "_L001_R2_001.fastq.gz")
    file5=paste0(expname, "/", F1andF2_both$Name[i],"_R", F1andF2_both$replicate.x[i],   "R1.fastq.gz")
    file6=paste0(expname, "/", F1andF2_both$Name[i],"_R", F1andF2_both$replicate.x[i],  "R2.fastq.gz")
    line1=paste0("cat ", file1, " ",  file3," > ", file5)
    line2=paste0("cat ", file2, " ",  file4," > ", file6)
    mytext=c(mytext, line1, line2)
    myfiles=c(myfiles, file5, file6)
  }
  noquote(mytext)
}

myfiles=combine_files(F1andF2_both)

fileConn=file("O:/JOLabShared/SaraDrescher/Molecular_clock/combine_files.txt")
writeLines(combine_files(F1andF2_both), fileConn)
close(fileConn)

#Then go to linux terminal, home/sdresche/JOLabShared/SaraDrescher/Molecular_clock/
# awk '{ sub ("\r$", ""); print }' combine_files.txt > combine_files.sh
# sbatch combine_files.sh
# or: sh combine_files.sh



combine_files2=function(F1andF2_both){
  line0="#!/bin/sh"
  mytext=noquote(line0)
  myfiles=vector()
  mytext2=vector()
  for (i in 1:length(F1andF2_both$Name)){
    samp1=F1andF2_both$Sample_ID.x[i]
    samp2=F1andF2_both$Sample_ID.y[i]
    expname="nano"
    file1=paste0(expname, "/", samp1, "_S", samp1, "_L001_R1_001.fastq.gz")
    file2=paste0(expname, "/", samp1, "_S", samp1, "_L001_R2_001.fastq.gz")
    file3=paste0(expname, "/", samp2, "_S", samp2, "_L001_R1_001.fastq.gz")
    file4=paste0(expname, "/", samp2, "_S", samp2, "_L001_R2_001.fastq.gz")
    file5=paste0(expname, "/", F1andF2_both$Name[i],"_R", F1andF2_both$replicate.x[i],   "R1.fastq.gz")
    file6=paste0(expname, "/", F1andF2_both$Name[i],"_R", F1andF2_both$replicate.x[i],  "R2.fastq.gz")
    line1=paste0("cat ", file1, " ",  file3," > ", file5)
    line2=paste0("cat ", file2, " ",  file4," > ", file6)
    line3=paste0("hivmmer --id ", paste0(F1andF2_both$Name[i], "_",F1andF2_both$replicate.x[i]) , "_F1F2 --fq1 ", file5, " --fq2 ", file6, " --ref PolACD.trimmed.aa --region int")
    mytext=c(mytext, line1, line2)
    myfiles=c(myfiles, line3)
  }
  noquote(mytext)
  return(myfiles)
}


text=combine_files2(F1andF2_both)
