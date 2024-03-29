# pull out fragments corresponding to fragment numbers from sample spreadsheet.
library(stringr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(gridExtra)
sample_names=read.csv("O:/JOLabShared/SaraDrescher/Molecular_clock/sample_names_nano.csv", stringsAsFactors = FALSE)

sample_names%<>%mutate(fragment=unlist(strsplit(Sample_Name, split="_"))[grepl("F.", unlist(strsplit(Sample_Name, split="_")))])
#this pulls the particular fragment out of the sample name

#nano data path: X:\\fast\\bloom_j\\SR\\ngs\\illumina\\sdresche\\190503_M03100_0427_000000000-CG8MD\\Data\\Intensities\\BaseCalls"
#second run data path: X:\\fast\\bloom_j\\SR\\ngs\\illumina\\sdresche\\191014_M04866_0296_000000000-CM883\\Data\\Intensities\\BaseCalls"
#unix path: /fh/fast/bloom_j/SR/ngs/illumina/sdresche/191014_M04866_0296_000000000-CM883/Data/Intensities/BaseCalls

#this function writes a file containing all the code to run HIVMMER (a list of commands, one set per sample)
#datapath is the location where the illumina reads are 
#sample_namesdb is the (already loaded) object containing the columns "fragment", "Sample_ID" "Sample_Names"
#this assumes that the samples are in order by sample number in the Sample_Names columns, as the reads are labeled by number, not name!
#if they aren't, you'll want to order by Sample_ID before beginning
#Otherwise, the column Sample_ID isn't used
#probably would be more bug proof if I checked this first--maybe for a later iteration


#mypath="/fh/fast/bloom_j/SR/ngs/illumina/sdresche/191014_M04866_0296_000000000-CM883/Data/Intensities/BaseCalls"
mypath="/home/sdresche/JOLabShared/SaraDrescher/Molecular_clock/nano"
write_HIVMMER_Command=function(sample_num, output_name=default_name, expname="", sample_namesdb, datapath) {
  outvector=vector()
  frag=sample_namesdb$fragment[sample_num]
  default_name=paste(sample_namesdb$Sample_Name[sample_num], sample_num, sep=".")
  outfile=paste(expname,  "results", default_name, sep="/")
  file1=paste0(datapath,"/", sample_num, "_S", sample_num, "_L001_R1_001.fastq.gz")
  file2=paste0(datapath,"/", sample_num, "_S", sample_num, "_L001_R2_001.fastq.gz")
  ref=ifelse(frag=="F1", "GagACD.trimmed.aa.hmm", "PolACD.trimmed.aa.hmm")
  out_name=paste("/home/sdresche/JOLabShared/SaraDrescher/Molecular_clock", outfile, sep="/")
  out_aavf=paste0( outfile, ".hmmsearch2.aavf")
  mycommand=paste0("hivmmer --id ", out_name," --fq1 ", file1," --fq2 ", file2," --ref ", ref, " --region int")
  noquote(mycommand)
}


write_HIVMMER_outfiles=function(sample_num, output_name=default_name, expname="", sample_namesdb, datapath) {
  outvector=vector()
  frag=sample_namesdb$fragment[sample_num]
  default_name=paste(sample_namesdb$Sample_Name[sample_num], sample_num, sep=".")
  outfile=paste(expname,  "results", output_name, sep="/")
  out_name=paste0("/home/sdresche/JOLabShared/SaraDrescher/Molecular_clock", outfile)
  out_aavf=paste0( outfile, ".hmmsearch2.aavf")
  outvector=c(outvector, out_aavf)
  return(outvector)
}

out_files_list=write_HIVMMER_outfiles(1:length(sample_names$Sample_Name),datapath=mypath,expname="nano", sample_namesdb=sample_names)
out_files_list_name=paste0("out_files_list_", expname, ".csv")
write.csv(out_files_list, out_files_list_name)

#hivmmer --id /home/sdresche/JOLabShared/SaraDrescher/p1_4_f1_Neher --fq1 /home/sdresche/JOLabShared/SaraDrescher/testout1.fastq --fq2  /home/sdresche/JOLabShared/SaraDrescher/testout2.fastq --ref GagACD.trimmed.aa.hmm --region int

#write_HIVMMER_Command(c(1,2,10,20,30,40,50,60,70,80,90), datapath=mypath, expname="nano", sample_namesdb=sample_names)

fileConn<-file("O:/JOLabShared/SaraDrescher/Code/HIVMMERCommands_nano.txt")
writeLines(c(noquote("#!/bin/sh"), noquote("cd /home/sdresche/JOLabShared/SaraDrescher/Molecular_clock"), write_HIVMMER_Command(18:19, expname="nano", sample_namesdb=sample_names, datapath=mypath)), fileConn)
close(fileConn)



#awk '{ sub("\r$", ""); print }'  ~/JOLabShared/SaraDrescher/Code/HIVMMERCommands_nano.txt> ~/JOLabShared/SaraDrescher/Code/HIVMMER_nano.sh
#sbatch HIVMMER_nano.sh
#awk '{ sub("\r$", ""); print }' pt_cons_blastcommands2.txt> pt_cons_blastcommands2.sh
#sbatch ~/JOLabShared/SaraDrescher/Code/HIVMMER_nano_filtered.sh
#c-36720921.txt
#also made a list of all the output AAVF files to use to do the analysis
#out_files_list=write_HIVMMER_Command(1:length(sample_names$Sample_Name), expname="nano/filtered", sample_namesdb=sample_names)
#out_files_list=list.files("O://JOLabShared//SaraDrescher//Molecular_clock//nano//filtered//results", pattern="hmmsearch2.aavf")


source("C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/Code_for_github/Molecular_clock/AAVF_3rdcodon_function.R")


sample_names=read.csv("O:/JOLabShared/SaraDrescher/Molecular_clock/sample_names_nano.csv", stringsAsFactors = FALSE)

sample_names%<>%mutate(fragment=unlist(strsplit(Sample_Name, split="_"))[grepl("F.", unlist(strsplit(Sample_Name, split="_")))])

#to make plots
pdf(file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/coverage_plots_nano.pdf")
coverage_min=5000
third_output=vector()
plots_list=list()
for (i in 1:length(sample_names$Sample_Name)){
  fragment=sample_names$fragment[i]
  filepath=paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/", out_files_list[i])
  x=AAVF.3rdcodon(filepath, fragment=fragment)
  atitle=paste0(sample_names$Sample_Name[i], " APD= ", round(as.numeric(x[[2]][2]), 4))
  print(x[[1]] + ggtitle(atitle) + geom_hline(yintercept=5000, color="red", linetype="dashed"))
  third_output=rbind(third_output, x[[2]])
  plots_list=c(plots_list, x[1])
}
dev.off()

#Same thing without the plots:
third_output=vector()
for (i in 1:length(sample_names$Sample_Name)){
  fragment=sample_names$fragment[i]
  filepath=paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/", out_files_list[i])
  x=AAVF.3rdcodon(filepath, fragment=fragment)
  myvector=cbind(x[[2]], i)
  third_output=rbind(third_output,myvector)
}


third_output=as.data.frame(third_output)

colnames(third_output)=c("Name", "APD", "AAstart", "AAend", "Good AAs", "Seq continuous", "HXB2nt_start", "HXB2nt_end", "consensus", "Fragment", "sample_num")
write.csv(third_output, file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/third_output_nano_unfiltered.csv")


