#!/bin/sh
pandaseq -f 68_S68_L001_R1_001.fastq.gz -r 68_S68_L001_R2_001.fastq.gz -w VA_1300_F1_30_cycles_R1_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F1_30_cycles_R1_merged.fasta -out VA_1300_F1_30_cycles_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 69_S69_L001_R1_001.fastq.gz -r 69_S69_L001_R2_001.fastq.gz -w VA_1300_F1_35-pre_gel_R1_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F1_35-pre_gel_R1_merged.fasta -out VA_1300_F1_35-pre_gel_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 70_S70_L001_R1_001.fastq.gz -r 70_S70_L001_R2_001.fastq.gz -w VA_1300_F1_35-pre_gel_R2_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F1_35-pre_gel_R2_merged.fasta -out VA_1300_F1_35-pre_gel_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 71_S71_L001_R1_001.fastq.gz -r 71_S71_L001_R2_001.fastq.gz -w pt_127_M9_F1_R1_merged.fasta
blastn -db pt_ref_DB -query pt_127_M9_F1_R1_merged.fasta -out pt_127_M9_F1_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 72_S72_L001_R1_001.fastq.gz -r 72_S72_L001_R2_001.fastq.gz -w pt_127_M9_F1_R2_merged.fasta
blastn -db pt_ref_DB -query pt_127_M9_F1_R2_merged.fasta -out pt_127_M9_F1_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 73_S73_L001_R1_001.fastq.gz -r 73_S73_L001_R2_001.fastq.gz -w pt_127_M9_F2_R1_merged.fasta
blastn -db pt_ref_DB -query pt_127_M9_F2_R1_merged.fasta -out pt_127_M9_F2_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 74_S74_L001_R1_001.fastq.gz -r 74_S74_L001_R2_001.fastq.gz -w pt_127_M9_F2_R2_merged.fasta
blastn -db pt_ref_DB -query pt_127_M9_F2_R2_merged.fasta -out pt_127_M9_F2_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 75_S75_L001_R1_001.fastq.gz -r 75_S75_L001_R2_001.fastq.gz -w pt_127_M9_F3_R1_merged.fasta
blastn -db pt_ref_DB -query pt_127_M9_F3_R1_merged.fasta -out pt_127_M9_F3_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 76_S76_L001_R1_001.fastq.gz -r 76_S76_L001_R2_001.fastq.gz -w pt_127_M9_F3_R2_merged.fasta
blastn -db pt_ref_DB -query pt_127_M9_F3_R2_merged.fasta -out pt_127_M9_F3_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 77_S77_L001_R1_001.fastq.gz -r 77_S77_L001_R2_001.fastq.gz -w VA_1300_F2_R1_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F2_R1_merged.fasta -out VA_1300_F2_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 78_S78_L001_R1_001.fastq.gz -r 78_S78_L001_R2_001.fastq.gz -w VA_1300_F2_R2_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F2_R2_merged.fasta -out VA_1300_F2_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 79_S79_L001_R1_001.fastq.gz -r 79_S79_L001_R2_001.fastq.gz -w VA_1300_F3-20_15_R1_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F3-20_15_R1_merged.fasta -out VA_1300_F3-20_15_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 80_S80_L001_R1_001.fastq.gz -r 80_S80_L001_R2_001.fastq.gz -w VA_1300_F3-25_10_R1_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F3-25_10_R1_merged.fasta -out VA_1300_F3-25_10_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 81_S81_L001_R1_001.fastq.gz -r 81_S81_L001_R2_001.fastq.gz -w VA_1300_F1_35extract_R1_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F1_35extract_R1_merged.fasta -out VA_1300_F1_35extract_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 82_S82_L001_R1_001.fastq.gz -r 82_S82_L001_R2_001.fastq.gz -w VA_1300_F1_35extract_R2_merged.fasta
blastn -db pt_ref_DB -query VA_1300_F1_35extract_R2_merged.fasta -out VA_1300_F1_35extract_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 83_S83_L001_R1_001.fastq.gz -r 83_S83_L001_R2_001.fastq.gz -w pt_116_M24_F1_R1_merged.fasta
blastn -db pt_ref_DB -query pt_116_M24_F1_R1_merged.fasta -out pt_116_M24_F1_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 84_S84_L001_R1_001.fastq.gz -r 84_S84_L001_R2_001.fastq.gz -w pt_116_M24_F1_R2_merged.fasta
blastn -db pt_ref_DB -query pt_116_M24_F1_R2_merged.fasta -out pt_116_M24_F1_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 85_S85_L001_R1_001.fastq.gz -r 85_S85_L001_R2_001.fastq.gz -w pt_116_M24_F2_R1_merged.fasta
blastn -db pt_ref_DB -query pt_116_M24_F2_R1_merged.fasta -out pt_116_M24_F2_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 86_S86_L001_R1_001.fastq.gz -r 86_S86_L001_R2_001.fastq.gz -w pt_116_M24_F2_R2_merged.fasta
blastn -db pt_ref_DB -query pt_116_M24_F2_R2_merged.fasta -out pt_116_M24_F2_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 87_S87_L001_R1_001.fastq.gz -r 87_S87_L001_R2_001.fastq.gz -w pt_116_M24_F3_R1_merged.fasta
blastn -db pt_ref_DB -query pt_116_M24_F3_R1_merged.fasta -out pt_116_M24_F3_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 88_S88_L001_R1_001.fastq.gz -r 88_S88_L001_R2_001.fastq.gz -w pt_116_M24_F3_R2_merged.fasta
blastn -db pt_ref_DB -query pt_116_M24_F3_R2_merged.fasta -out pt_116_M24_F3_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 89_S89_L001_R1_001.fastq.gz -r 89_S89_L001_R2_001.fastq.gz -w pt_127_M18_F1_R1_merged.fasta
blastn -db pt_ref_DB -query pt_127_M18_F1_R1_merged.fasta -out pt_127_M18_F1_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 90_S90_L001_R1_001.fastq.gz -r 90_S90_L001_R2_001.fastq.gz -w pt_127_M18_F1_R2_merged.fasta
blastn -db pt_ref_DB -query pt_127_M18_F1_R2_merged.fasta -out pt_127_M18_F1_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 91_S91_L001_R1_001.fastq.gz -r 91_S91_L001_R2_001.fastq.gz -w pt_127_M18_F2_R1_merged.fasta
blastn -db pt_ref_DB -query pt_127_M18_F2_R1_merged.fasta -out pt_127_M18_F2_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 92_S92_L001_R1_001.fastq.gz -r 92_S92_L001_R2_001.fastq.gz -w pt_127_M18_F2_R2_merged.fasta
blastn -db pt_ref_DB -query pt_127_M18_F2_R2_merged.fasta -out pt_127_M18_F2_R2_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 93_S93_L001_R1_001.fastq.gz -r 93_S93_L001_R2_001.fastq.gz -w pt_127_M18_F3_R1_merged.fasta
blastn -db pt_ref_DB -query pt_127_M18_F3_R1_merged.fasta -out pt_127_M18_F3_R1_blastout.txt -outfmt 6 -max_target_seqs 1
pandaseq -f 94_S94_L001_R1_001.fastq.gz -r 94_S94_L001_R2_001.fastq.gz -w pt_127_M18_F3_R2_merged.fasta
blastn -db pt_ref_DB -query pt_127_M18_F3_R2_merged.fasta -out pt_127_M18_F3_R2_blastout.txt -outfmt 6 -max_target_seqs 1
