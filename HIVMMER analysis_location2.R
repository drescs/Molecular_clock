# pull out fragments corresponding to fragment numbers from sample spreadsheet.
library(stringr)
library(tidyverse)
library(magrittr)
sample_names=read.csv("O:/Documents/molecular_clock/sample_names.csv", stringsAsFactors = FALSE)

sample_names%<>%mutate(fragment=unlist(strsplit(Sample_Name, split="_"))[grepl("F.", unlist(strsplit(Sample_Name, split="_")))])
#this is pretty badly written as you have to load sample names before you call this function or it won't work...
#tried to fix it by requiring the sample_names to be passed as sample_namesdb
#nano data path: X:\\fast\\bloom_j\\SR\\ngs\\illumina\\sdresche\\190503_M03100_0427_000000000-CG8MD\\Data\\Intensities\\BaseCalls"
#second run data path: X:\\fast\\bloom_j\\SR\\ngs\\illumina\\sdresche\\191014_M04866_0296_000000000-CM883\\Data\\Intensities\\BaseCalls"
#unix path: /fh/fast/bloom_j/SR/ngs/illumina/sdresche/191014_M04866_0296_000000000-CM883/Data/Intensities/BaseCalls

#this function writes a file containing all the code to run HIVMMER (a list of commands, one set per sample)
write_HIVMMER_Command=function(sample_num, output_name=default_name, expname="", sample_namesdb, datapath) {
  outvector=vector()
  frag=sample_namesdb$fragment[sample_num]
  default_name=paste(sample_namesdb$Sample_Name[sample_num], sample_num, sep=".")
  outfile=paste(expname,  "results", output_name, sep="/")
  file1=paste0(expname, sample_num, "_S", sample_num, "_L001_R1_001.fastq.gz")
  file2=paste0(expname, sample_num, "_S", sample_num, "_L001_R2_001.fastq.gz")
  ref=ifelse(frag=="F1", "GagACD.trimmed.aa.hmm", "PolACD.trimmed.aa.hmm")
  out_name=paste0("/home/sdresche/JOLabShared/SaraDrescher/Molecular_clock", outfile)
  out_aavf=paste0( outfile, ".hmmsearch2.aavf")
  outvector=c(outvector, out_aavf)
  mycommand=paste0("hivmmer --id ", out_name," --fq1 ", file1," --fq2 ", file2," --ref ", ref, " --region int")
  noquote(mycommand)
  return(outvector)
}

#hivmmer --id /home/sdresche/JOLabShared/SaraDrescher/p1_4_f1_Neher --fq1 /home/sdresche/JOLabShared/SaraDrescher/testout1.fastq --fq2  /home/sdresche/JOLabShared/SaraDrescher/testout2.fastq --ref GagACD.trimmed.aa.hmm --region int

#write_HIVMMER_Command(1, expname="nano/filtered/filtered_" )

#fileConn<-file("O:/Documents/molecular_clock/HIVMMERCommands2.txt")
#writeLines(c(noquote("#!/bin/sh"), write_HIVMMER_Command(1:length(sample_names$Sample_Name), expname="filtered/filtered_", datapath=currentpath)), fileConn)
#close(fileConn)

#also made a list of all the output AAVF files to use to do the analysis
out_files_list=write_HIVMMER_Command(1:length(sample_names$Sample_Name), expname="nano", sample_namesdb=sample_names)


source("O:/Documents/Code/Sequence Analysis/AAVF_parser_function.R")
source("O:/Documents/Code/Sequence Analysis/AAVF_3rdcodon_function.R")


#filepath="O:/Documents/molecular_clock/nano/results/pt_93_M6_F1_R1.1.hmmsearch2.aavf"
parsed_output=vector()

for (i in 1:94){
  fragment=sample_names$fragment[i]
  filepath=paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/", out_files_list[i])
  print(filepath)
  x=c(AAVF.parser(filepath, fragment=fragment, codon_variant_frac_min=0.05)[[1]], fragment, i)
  parsed_output=rbind(parsed_output, x)
}

parsed_output=as.data.frame(parsed_output)
colnames(parsed_output)=c("Name", "APD", "Sites_covered", "Is_continuous", "HXB2nt_start", "HXB2nt_end", "Fragment", "sample_num")
write.csv(parsed_output, file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/synonymous_5000_01_new.csv")

parsed_output%<>%separate(Name, c("Name", "Replicate"), sep="R")




third_output=vector()

for (i in 1:94){
  fragment=sample_names$fragment[i]
  filepath=paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/", out_files_list[i])
  x=c(AAVF.3rdcodon(filepath, fragment=fragment), i )
  third_output=rbind(third_output, x)
}

third_output=as.data.frame(third_output)
colnames(third_output)=c("Name", "APD", "HXB2nt_start", "HXB2nt_end", "consensus", "Fragment", "sample_num")
write.csv(third_output, file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/third_output_10_14.csv")




samp_num=7
coverage_min=5000
codon_variant_cov_min=50
codon_variant_frac_min=0.01

#AAVF.parser=function(filepath, codon_variant_frac_min=0.01, coverage_min=5000, fragment)
#AAVF.3rdcodon=function(filepath, fragment, codon_variant_frac_min=0.01, coverage_min=5000) 
filepath=paste0("O:/Documents/molecular_clock/", out_files_list[samp_num])
output=AAVF.parser(filepath, fragment=fragment)
output3rd=AAVF.3rdcodon(filepath, fragment=fragment)

x=AAVF.parser(paste0("O:/Documents/molecular_clock/", out_files_list[samp_num]), sample_names$fragment[samp_num])
x

