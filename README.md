# Molecular_clock
Code for analyzing miseq data used for the molecular clock project

Right now I've just thrown all the code up; it'll need a lot of reorganizing before someone can use it but gotta start somewhere.

Basic premise:

Using miseq data, perform some quality control (eliminating cross contamination, etc) calculate average pairwise distance for patient viral samples

This depends on using HIVMMER from the Kantor Lab (https://github.com/kantorlab). 

Right now it requires going back and forth between R and console. Probably there's a way to run it all via R or all via the console. Someone smarter than me could probably do that quickly. 

#Step 1: Build hivmmer index for HIV sequences/areas of interest: using make_ref_alignments.fasta 

#Step 2: run hivmmer on the files from miseq (needs two files, one forward and one reverse, for each)--can use HIVMMER analysis location 2 to write this code in bulk if it's helpful, but probably it won't be unless all your sample names are stored just like mine. (CSV with sample number in one column, sample name in another). The files are input using the original name (e.g., 1_S1_001_L001_R1_01.fastq.gz, but the results use the sample name instead, e.g., pt93_T1_F1_R1). That's the main helpful purpose of writing the code in R. 
Write HIVMMER command also makes an object called "outputfiles" which is the list (well, vector in R!) of output files you will want to use later, IN ORDER BY SAMPLE NUMBER. The reason that last is important is if you just grab all the files in the results directory with the right file ending, they won't be in order by sample number (if you sort by time changed they may be out of order if you ran it on multiple nodes at once). The sample number is used later in the code so you'll want to use the script to write the vector of output files.

#Step 3: Calculate APD. You'll get several files from each sample you put through HIVMMER. The file of interest is samplename.hmmsearch2.aavf. The aavf files is in a weird format and you'll want to parse that--use an AAVF_parser. The most current is AAVF_parser_3rd_codon. (There is a bit of code to source and call this in write_HIVMMER_commands). It calculates average pairwise distance (APD) (and a few other metrics) using the third codon position only--so it looks at quasi-synonymous changes, consistant with the method used by Puller et all to calculate APD. Another version, AAVF_synonymous, will give you something quite similar but using synonymous codons (compared to reference). The difference is subtle and likely doesn't matter but it's important to know which you used.

#Step: 4 To check that all reads with a particular barcode are truly from the sample in question and not contamination from another sample in the batch, we put all of the consensus sequence that came out of hivmmer (they are hidden in that AAVF file) and make a local blast database with them. Then we reblast all the sequences against these and remove any that align to anythign but the correct consensus. =First you need to make a the blast database--that's in make consensus fasta.R. This makes the blast database and then writes the blast commands to a text file. A tricky bit is that the blast software (from BLAST+) requires you to first merge the paired-end reads, so that is in there too. Set the blast to allow a maximum of one (the best) match.

If like me, you are running R studio on windows, it will write th carriage returns differently than linux. Thus, before running the commands you need to run:
```
awk '{ sub("\r$", ""); print }' Commands.txt> Commands.sh
```
(replace Commands.txt with whatever file contains your blast commands).

#Step 5: After you've blasted everything on the command line, we need to analyze the output. That's in the second half of make consensus fasta.R. It will look at what sequence each read aligns to best, and pull out the identifiers of th reads that aligned to the correct sequences (and the correct fragment). This will get all mixed up for controls that use mixtures of pure strains--but we know that and can ignore that output. 
Then, it imports the original sequencing data, filters for only those that align to the correct consensus, and writes new forward and reverse fastq files (titled filtered_originialfilename.fastq.gz).

Step 6. When this is compelete, the next step would be to run these new, decomtaminated files in HIVMMER and analyze the new hmmsearch2.aavf files using AAVF parser again.

In parallel: In addition to this we'll want to see if the controls we made (mixes of Q23 (subtype A) and LAI (subtype B, also called HXB2) worked properly. control_percent_analysis does this: first make the "local blast database" which is just a fasta file with a Q23 and HXB2 sequence in it, then run makeblastdb on that on the command line. Then you need to merge the illumina reads using pandaseq and then align them to the database. Then the analysis will tell you what percent of the aligned reads for that sample aligned to Q23.


