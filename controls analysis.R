#analyse controls

#first step will be to find the places where Q23 and HXB2 differ

#would need to align them for this

library(Biostrings)
library(seqinr)

Q23=read.fasta("O:/Documents/molecular_clock/Q23.fa")
HXB2=read.fasta("O:/Documents/molecular_clock/nano/K03455.fasta")

pairwiseAlignment(Q23, HXB2)