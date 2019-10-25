#Let's play with HIVMMER output!
library(readr)
library(stringr)
library(magrittr)
library(tidyverse)
library(Biostrings)
library(seqinr)
#update these variables for each analysis
filepath="O:/JOLabShared/SaraDrescher/Molecular_clock/nano/results/pt_93_M6_F2_R1.3.hmmsearch2.aavf"
#consensusfilepath="O:/Documents/molecular clock/HMsample05.hmmsearch2.ref.consensus.fasta"
coverage_min=500
codon_variant_cov_min=50
codon_variant_frac_min=0.01
fragment="F2"

AAVF_metadata=read_lines(filepath, n_max=10)
AAVF_raw=read.table(filepath, sep="\t", header=FALSE)
AAVF_names=unlist(strsplit(AAVF_metadata[10], "\\t"))
AAVF_names[1]="CHROM"

#pull out the name of the sample from line 4
sample_name=unlist(strsplit((unlist(strsplit(AAVF_metadata[4], "[.]"))[1]),"="))[2]

colnames(AAVF_raw)=AAVF_names
AAVF_raw$POS=as.numeric(AAVF_raw$POS)
AAVF_raw$COVERAGE=as.numeric(AAVF_raw$COVERAGE)

#For each AA, this lists the number of codons found. There are up to 6 codons per AA (leucine, arginine) so we'll split this column into 6. Unused columns are filled with NA to the right.
AAVF_raw%<>%separate(INFO, c("AC", "ACC", "ACF"), sep=";")
AAVF_raw=separate(AAVF_raw, AC, c("AC1", "AC2", "AC3", "AC4", "AC5", "AC6"), sep=",")
Variants=separate(AAVF_raw, ACC, c("ACC1", "ACC2", "ACC3", "ACC4", "ACC5", "ACC6"), sep=",")
Variants=separate(Variants, ACF, c("ACF1", "ACF2", "ACF3", "ACF4", "ACF5", "ACF6"), sep=",")
#need to remove strings "AC=", "ACC=" "ACF=" from the first set of columns
Variants$AC1=str_replace(Variants$AC1, "AC=", "")
Variants$ACC1=str_replace(Variants$ACC1, "ACC=", "")
Variants$ACF1=str_replace(Variants$ACF1, "ACF=", "")

#right now covereage is total coverage at a site. Since I'm only analyzing synonymous mutations, I want to only look at synonmous codons and how ofthen they are sequenced
#At least, right now I think so
#So in addition to per site coverage, I want to calculate per codon coverage.
#that's not included in initial file so I'll make it


Variants$ACC1%<>%as.numeric()
Variants$ACC2%<>%as.numeric()
Variants$ACC3%<>%as.numeric()
Variants$ACC4%<>%as.numeric()
Variants$ACC5%<>%as.numeric()
Variants$ACC6%<>%as.numeric()

Variants%<>% mutate(Syn_Coverage=rowSums(Variants[,c("ACC1", "ACC2", "ACC3", "ACC4", "ACC5", "ACC6")],  na.rm=T))



#Need to restructure this data and I'm getting annyoed messing around with reshape and melt and whatnot. Stata would do an awesome job here but whatever.
#I can just make 6 copies of this data frame, each including only the three relevant end columns, then stack them
#How annoying.
#Ok maybe we can combine columns and then use gather, and then use separate
#annoyingly this converts all my numeric values back to string but I'll fix that later.

Variants%<>%unite(Codon_1, c(AC1, ACC1, ACF1), sep = "_", remove = TRUE)
Variants%<>%unite(Codon_2, c(AC2, ACC2, ACF2), sep = "_", remove = TRUE)
Variants%<>%unite(Codon_3, c(AC3, ACC3, ACF3), sep = "_", remove = TRUE)
Variants%<>%unite(Codon_4, c(AC4, ACC4, ACF4), sep = "_", remove = TRUE)
Variants%<>%unite(Codon_5, c(AC5, ACC5, ACF5), sep = "_", remove = TRUE)
Variants%<>%unite(Codon_6, c(AC6, ACC6, ACF6), sep = "_", remove = TRUE)

Variants %<>% gather(key=codon_number, value=codon_info, Codon_1, Codon_2, Codon_3, Codon_4, Codon_5, Codon_6)

Variants=separate(Variants, codon_info, c("Codon", "Codon_coverage", "Codon_frequency"), sep="_")
Variants=Variants[which(Variants$Codon!="NA"),]

Variants$Codon_coverage=as.numeric(Variants$Codon_coverage)
Variants$Codon_frequency=as.numeric(Variants$Codon_frequency)
Variants%<>% group_by(POS, ALT) %>% mutate(dominant_allele_coverage=max(Codon_coverage,na.rm=T))
Variants%<>%mutate(freq_in_syn=Codon_coverage/Syn_Coverage)
Variants%<>%mutate(freq_syn_dom=dominant_allele_coverage/Syn_Coverage)

#Select only those that pass the filter
Pass_vars=Variants[which(Variants$FILTER=="PASS"),]
Pass_vars=Variants[which(Variants$COVERAGE>=coverage_min),]
#use this to ignore sites where diversity doesn't meet our minimum variant calling threshhold
Pass_vars%<>%mutate(to_use=ifelse(1-freq_syn_dom-codon_variant_frac_min > 0, 1, 0 )) 
Pass_vars%<>%mutate(contribution=freq_in_syn*(1-freq_in_syn)*to_use)
#Now only those that are synonymous
Syn_pass_vars=Pass_vars[which(Pass_vars$REF==Pass_vars$ALT),]
#need to know number of sites to divde by
covered_sites=length(unique(Syn_pass_vars$POS))

APD=sum(Syn_pass_vars$contribution)/covered_sites



#the coverage column lists overall sequencing depth at each position
#We want to just check quickly that that is reasonable--but it's repeated if there are different AAs at the same pos, so we just want unique position/coverage values
POS_COVERAGE=unique(AAVF_raw[,c("POS", "COVERAGE", "REF")])
POS_COVERAGE%<>%arrange(POS)
POS_COVERAGE$REF=as.character(POS_COVERAGE$REF)

#plot coverage by positionls
covplot=ggplot(POS_COVERAGE, aes(y=COVERAGE, x=POS)) + 
  geom_bar(position="dodge", stat="identity") +  
  theme_bw()
covplot
consensus=paste(POS_COVERAGE$REF, collapse="")
good_consensus=paste(POS_COVERAGE$REF[(which(POS_COVERAGE$COVERAGE>coverage_min))], collapse="")

#make sure our sequence of amino acids with significant coverage is continuous, and extract 1st and last codons with good covereage, relative to 
#sample specific reference
check_continuous=function (POS_COVERAGE, coverage_min) {
  goodAAs=POS_COVERAGE$POS[(which(POS_COVERAGE$COVERAGE>coverage_min))]
  start_pos=min(goodAAs)
  print (paste0("starting position ", start_pos))
  end_pos=max(goodAAs)
  print (paste0("ending position ", end_pos))
  print(paste("There are",  length(goodAAs), "amino acids with >5000x coverage for this sample"))
  is_cont=ifelse(end_pos-start_pos+1==length(goodAAs), "Continuous sequence", "Not continuous sequence" )
  seq_ok=ifelse(end_pos-start_pos+1==length(goodAAs), 1, 0)
  print(is_cont)
  return_vector=c(start_pos, end_pos, length(goodAAs), seq_ok)
}

check_vector=check_continuous(POS_COVERAGE, coverage_min)

#translate those start and end codons to the corresponding HXB2 nucleotide start and stop positions
#this is important for plugging them into Neher's clock, found at:
#https://hiv.biozentrum.unibas.ch/ETI/

#now to figure out where our sample specific consensus lives relative to HXB2. This was tricky, I suspect I'm missing something easy.
#but, I aligned them
#found where the aligned portion of HXB2 sits in the HXB2 translated pol gene 
#(use str split, size of the first split is how many AAs skipped before alignment starts.)
#take that length +1 to get first aa, multiply by 3 to get nt position, then add nt position of beginning of gene, subtracting 1..



get_start_HXB2_AA_coord=function(consensus, frag){
  if (frag=="F1"){
    #this name might not be right!
    myRef=read.fasta(file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/gag_HXB2_aa.fasta")
  } else {
    myRef=read.fasta(file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/HXB2_pol_aa.fasta")
  }
  consens_string=toupper(c2s(consensus[[1]]))
  myref_string=toupper(c2s(myRef[[1]]))
  myalign=pairwiseAlignment(consens_string, myref_string, type="local")
  myalign
  search_string=as.character(subject(myalign))
  HXB2aastart=nchar(str_split(myref_string, search_string)[[1]][1])
  return(HXB2aastart)
}

HXB2aastart=get_start_HXB2_AA_coord(good_consensus, fragment)

#get hxb2 gene start site depending on fragment
gene_start=ifelse(fragment=="F1", 790, 2085)

HXB_nt_start=(HXB2aastart+1)*3+gene_start-1
HXB_nt_end=HXB_nt_start + check_vector[3]*3
print(HXB_nt_start)
print(HXB_nt_end)
print(APD)

#So need to determine when to include a variant in my calculation and how to calculate "diversity"--will we use fraction of a particular codon/total codons for that AA?
#seems best, will need to add this up as we get total coverage, which includes other AA variants 

myRef=read.fasta(file="C:/Users/Sara/Documents/HHMI_projects/Molecular_clock/Analysis/HXB2_pol_aa.fasta")









