#make pt specific consensus fasta for use in BLAST
library(magrittr)
library(seqinr)
library(ShortRead)
library(tidyverse)

#load the consensus sequences that are in your ouput from the AAVF parser
  my_cons=read.csv("C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/third_output_10_14.csv")
  my_cons%<>%mutate(fasta_names=paste0(Name, ".", sample_num))

#make a fasta file from these that will serve as your local blast database 
  write.fasta(as.list(my_cons$consensus), my_cons$Name, open = "w", file.out="O:/JOLabShared/SaraDrescher/Molecular_clock/second_run/pt_ref_sequences2.fasta", as.string=TRUE)

#make the blast database in linux:
# makeblastdb -in pt_ref_sequences2.fasta -dbtype nucl -parse_seqids -out pt_ref_DB2
  
#write the commands you will use to merge your files and  then blast them against your local database
#note that the blast command requires merged files
    
  cons_names=my_cons[,c("Name","sample_num")]
  index=1:94
  commands=vector()
  outputfiles=vector()
  forward_list=vector()
  reverse_list=vector()
  for (i in (1:length(cons_names[,1]))){
    sample_num=cons_names[i, 2]
    sample_name=cons_names[i,1]
    path="X:/fast/bloom_j/SR/ngs/illumina/sdresche/191014_M04866_0296_000000000-CM883/Data/Intensities/BaseCalls/"
    file1=paste0( sample_num, "_S", sample_num, "_L001_R1_001.fastq.gz")
    file2=paste0( sample_num, "_S", sample_num, "_L001_R2_001.fastq.gz")
    command1=paste0("pandaseq -f ", path, file1," -r ",path, file2," -w ", sample_name, "_merged.fasta")
    command2=paste0("blastn -db pt_ref_DB2 -query ", sample_name, "_merged.fasta -out ",sample_name,"_blastout.txt -outfmt 6 -max_target_seqs 1")
    output=paste0(sample_name,"_blastout.txt")
    outputfiles=c(outputfiles,output)
    commands=c(commands, command1, command2)
    forward_list=c(forward_list, file1)
    reverse_list=c(reverse_list, file2)
  }
  

  
  
  files_list=as.data.frame(cbind(cons_names[,2], as.character(cons_names[,1]), forward_list, reverse_list))
  
  write.csv(files_list, "O:/JOLabShared/SaraDrescher/Molecular_clock/sample_names_and_file_names.csv")
  #write these commands to a text file
  write_lines(c(noquote("#!/bin/sh"), noquote("cd ~/JOLabShared/SaraDrescher/Molecular_clock/second_run"), commands), "O:/JOLabShared/SaraDrescher/Code/pt_cons_blastcommands.txt")
  
#bash this file in the linux console on Nomachine--don't forget to use slurm or it takes forever
#need to make sure sequences align to the right fragment and patient; but timepoint doesn't matter.
#unfortunately I named the files poorly so have to do a lot of text parsing to get that info from the
#samples and controls into the same column
  
#this is something that will likely change (hopefully become simpler) as I name future files more carefully.
mycols=as.data.frame(files_list)

#second_run
# mycols %<>% separate(V2, c("Pt_name", "rest"), sep="_F", remove=FALSE)
# mycols %<>% separate(rest, c("Frag", "rest"), sep="_R")
# mycols %<>% separate(rest, c("Rep", "Time"), sep="_[T,M]", fill="right")
# mycols %<>% separate(Pt_name, c("Pt_name", "Gel_extracted"), sep="gel", fill="right")
# mycols %<>% separate(Pt_name, c("percent", "Pt_name"), sep="%", fill="left")

#nano
mycols %<>% separate(V2, c("Mix", "therest"), sep="ix_", remove=FALSE, fill="left")
mycols %<>% separate(therest, c("therest"), sep="_blastout.txt")
mycols %<>% separate(therest, c("therest", "Rep"), sep="_R", fill="right")
mycols %<>% separate(therest, c( "therest", "Frag"), sep=("_F"), fill="right")
mycols %<>% separate(therest, c("P", "therest"), sep="t_", remove=FALSE, fill="left")
mycols %<>% separate(therest, c( "Pt_name", "Month"), sep=("_M"), remove=FALSE, fill="right")
mycols %<>% separate(Frag, c( "Frag", "more"), sep=("_"), remove=FALSE, fill="right")
mycols %<>% separate(Frag, c( "Frag", "evenmore"), sep=("-"), remove=FALSE, fill="right")
mycols%<>%select(V1, V2, Frag, Month, Pt_name, Rep, forward_list, reverse_list)

totals=vector()
not_contam=vector()
for (i in 1:94){
blast=read.delim(paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/nano/",outputfiles[i]), sep="\t", header=FALSE)
x=(unique(blast[,1:2]))

#nano
#mycols %<>% separate(V2, c("Mix", "therest"), sep="ix_", remove=FALSE, fill="left")
#mycols %<>% separate(therest, c("therest"), sep="_blastout.txt")
x %<>% separate(V2, c("Mix", "therest"), sep="ix_", remove=FALSE, fill="left")
x %<>% separate(therest, c("therest"), sep="_blastout.txt")
x %<>% separate(therest, c("therest", "Rep"), sep="_R", fill="right")
x %<>% separate(therest, c( "therest", "Frag"), sep=("_F"), fill="right")
x %<>% separate(therest, c("P", "therest"), sep="t_", remove=FALSE, fill="left")
x %<>% separate(therest, c( "Pt_name", "Month"), sep=("_M"), remove=FALSE, fill="right")
x %<>% separate(Frag, c( "Frag", "more"), sep=("_"), remove=FALSE, fill="right")
x %<>% separate(Frag, c( "Frag", "evenmore"), sep=("-"), remove=FALSE, fill="right")


#second_run
# x %<>% separate(V2, c("Pt_name", "rest"), sep="_F", remove=FALSE)
# x %<>% separate(rest, c("Frag", "rest"), sep="_R")
# x %<>% separate(rest, c("Rep", "Time"), sep="_[T,M]", fill="right")
# x %<>% separate(Pt_name, c("Pt_name", "Gel_extracted"), sep="gel", fill="right")
# x %<>% separate(Pt_name, c("percent", "Pt_name"), sep="%", fill="left")
good_seqs=x[which(x$Pt_name==mycols$Pt_name[i] & x$Frag==mycols$Frag[i]),]
reverses=readFastq((paste0(path, as.character(files_list$reverse_list[i]))))
forwards=readFastq((paste0(path, as.character(files_list$forward_list[i]))))
mynames <-  as.character(ShortRead::id(forwards))
mynamesr <-  as.character(ShortRead::id(reverses))
mynames_f=(gsub(" [1-2]:N:0:[0-9]*", "", mynames))
mynames_r=(gsub(" [1-2]:N:0:[0-9]*", "", mynamesr))
z=nchar(i)+1
if(z==2){
good_seq_ids=gsub('.{2}$', '', as.character(good_seqs$V1))
} else if (z==3) {
  good_seq_ids=gsub('.{3}$', '', as.character(good_seqs$V1))
}

#these actually should be identical unless something is really wrong!
good_for_indices=which(mynames_f%in%good_seq_ids)
good_rev_indices=which(mynames_r%in%good_seq_ids)
totals=c(totals, length(forwards))
not_contam=c(not_contam, length(good_for_indices))
#forwards1 <- forwards[good_for_indices]
#reverses1 <- reverses[good_rev_indices]
#writeFastq(forwards1, paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/second_run/filtered/filtered_", as.character(files_list$forward_list[i])), full=TRUE)
#writeFastq(reverses1, paste0("O:/JOLabShared/SaraDrescher/Molecular_clock/second_run/filtered/filtered_", as.character(files_list$reverse_list[i])), full=TRUE)
}