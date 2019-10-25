# Molecular_clock
Code for analyzing miseq data used for the molecular clock project

Right now I've just thrown all the code up; it'll need a lot of reorganizing before someone can use it but gotta start somewhere.

Basic premise:

Using miseq data, perform some quality control (eliminating cross contamination, etc) calculate average pairwise distance for patient viral samples

This depends heavily on the use of HIVMMER from the Kantor Lab (https://github.com/kantorlab). 

Right now it requires going back and forth between R and console. Probably there's a way to run it all via R or all via the console. Someone smarter than me could probably do that quickly. I'll try it if I have time.

Step 1: Build hivmmer index for HIV sequences/areas of interest: using make_ref_alignments.fasta 
run hivmmer on the files from miseq (needs two files, one forward and one reverse, for each)
