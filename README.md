# Molecular_clock
Code for analyzing miseq data used for the molecular clock project

Right now I've just thrown all the code up; it'll need a lot of reorganizing before someone can use it but gotta start somewhere.

Basic premise:

Using miseq data, perform some quality control (eliminating cross contamination, etc) calculate average pairwise distance for patient viral samples

This depends on using HIVMMER from the Kantor Lab (https://github.com/kantorlab). 

Right now it requires going back and forth between R and console. Probably there's a way to run it all via R or all via the console. Someone smarter than me could probably do that quickly. I'll try it if I have time.

Step 1: Build hivmmer index for HIV sequences/areas of interest: using make_ref_alignments.fasta 
Step 2: run hivmmer on the files from miseq (needs two files, one forward and one reverse, for each)--can use HIVMMER analysis location 2 to write this code in bulk if it's helpful, but probably it won't be unless all your sample names are stored just like mine. (CSV with sample number in one column, sample name in another).
Step 3: Calculate APD You'll get several files from each sample you put through HIVMMER. The file of interest is samplename.hmmsearch2.aavf. The aavf files is in a weird format and you'll want to parse that--use an AAVF_parser. The most current is AAVF_parser_3rd_codon. It calculates average pairwise distance (APD) (and a few other metrics) using the third codon position only--so it looks at quasi-synonymous changes, consistant with the method used by Puller et all to calculate APD. Another version, AAVF_synonymous, will give you something quite similar but using synonymous codons (compared to reference). The difference is subtle and likely doesn't matter but it's important to know which you used.
Step: 4 To check that all reads with a particular barcode are truly from the sample in question and not contamination from another sample in the batch, we put all of the consensus sequence that came out of hivmmer (they are hidden in that AAVF file) and make a local blast database with them. Then we reblast all the sequences against these and remove any that align to anythign but the correct consensus. Step 1 is to make a the blast database--that's in make consensus fasta.R. This makes the blast database and then writes the blast commands to a text file. A tricky bit is that the blast software (from BLAST+) requires you to first merge the paired-end reads, so that is in there too. Set the blast to allow a maximum of one (the best) match.

If like me, you are running R studio on windows, it will write th carriage returns differently than linux. Thus, before running the commands you need to run:
'''
awk '{ sub("\r$", ""); print }' Commands.txt> Commands.sh
'''
(replace Commands.txt with whatever file contains your blast commands).

Step 5: After you've blasted everything on the commands line, we need to analyze the output. That's in the second half of make consensus fasta.R. It will look at what sequence each read aligns to best, and pull out the identifiers of th reads that aligned to the correct sequences (and the correct fragment). This will get all mixed up for controls that use mixtures of pure strains--but we know that and can ignore that output. 
Then, it imports the original sequencing data, filters for only those that align to the correct consensus, and writes new forward and reverse fastq files (titled filtered_originialfilename.fastq.gz).

When this is compelete, the next step would be to run these new, decomtaminated files in HIVMMER and analyze the new hmmsearch2.aavf files using AAVF parser again.

