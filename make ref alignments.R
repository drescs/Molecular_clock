#modify a lanl reference compendium to include A, C, D and HXB2(for indexing)
library (seqinr)

make_ref=function(filename, outfilename){
fullfasta=read.fasta (file=filename)
subtypes=(strsplit(names(fullfasta), "\\."))
subtype.vector= sapply (subtypes, "[[", 1)
#gotta keep HXB2 even though it's the wrong subtype
HXB2=grep("HXB2",(names(fullfasta)))

#the whole one is too big and kills the hmm build command, so i will subset it.

use.samples=c(HXB2, sample(which(subtype.vector %in% c("A1", "A2", "C", "D")), 400))
subtypesACD.fasta=fullfasta[use.samples]
write.fasta(sequences=subtypesACD.fasta, names=names(subtypesACD.fasta), file.out=outfilename)
}

gag_ref="O:/Documents/molecular_clock/HIV1_ALL_2010_790-2350_DNA.fasta"
gag_outfile="O:/Documents/molecular_clock/GagACD_w_HXB2_ref.fasta"

pol_ref="O:/Documents/molecular_clock/HIV1_ALL_2010_pol_DNA.fasta"
pol_outfile="O:/Documents/molecular_clock/PolACD_w_HXB2_ref.fasta"

make_ref(pol_ref, pol_outfile)
make_ref(gag_ref, gag_outfile)

#to build a HIVMMER index file, use the command:
#hivmmer-trim-reference PolACD_w_HXB2_ref.fasta > PolACD.trimmed.aa.fasta
#hmmbuild  PolACD.trimmed.aa.hmm PolACD.trimmed.aa.fasta


#Next step is to use the HIVMMER command:

#hivmmer --id ID --fq1 FASTQ1 --fq2 FASTQ2 --ref REFERENCE  --region int

#ID specifies a name for the analysis that will be used as the basename for all output.
#FASTQ1 and FASTQ2 are the forward and reverse Illumina reads.
#seems odd but always put the region as "int." Then it does nothing to the output.
#If you put "pol" it filters for known resistance mutations which isn't useful for looking at APD
