
#see if I can pull the useful info out of the bowtie console output

# sample1  <------THIS

# 120183 reads; of these:    <------THIS

#   120183 (100.00%) were paired; of these:
#   55336 (46.04%) aligned concordantly 0 times
# 64739 (53.87%) aligned concordantly exactly 1 time
# 108 (0.09%) aligned concordantly >1 times
# ----
#   55336 pairs aligned concordantly 0 times; of these:
#   19706 (35.61%) aligned discordantly 1 time
# ----
#   35630 pairs aligned 0 times concordantly or discordantly; of these:
#   71260 mates make up the pairs; of these:
#   70332 (98.70%) aligned 0 times
# 848 (1.19%) aligned exactly 1 time
# 80 (0.11%) aligned >1 times
# 70.74% overall alignment rate <------THIS

#Where's the text output saved? a slurm.out or a file you piped it to.
bowtie_txt_path="O:/Documents/bowtie_output/slurm-36817750.out"
bowtietext=read_lines(bowtie_txt_path)

find_names=grepl("sample", bowtietext)
find_total_reads=grepl("reads; of these", bowtietext)
find_overall_aligned=grepl("overall alignment rate", bowtietext)

btsamplenames=bowtietext[which(find_names)]
total_reads=bowtietext[which(find_total_reads)]
overall_aligned=bowtietext[which(find_overall_aligned)]

sampleNumber=unlist(strsplit(btsamplenames, split="sample"))[1:length(btsamplenames)*2]
total_reads_num=unlist(strsplit(total_reads, split=" "))[(1:length(total_reads))*4-3]
overall_aligned_num=unlist(strsplit(overall_aligned, split=" "))[(1:length(total_reads))*4-3]


bowtie_text_info=data.frame(sampleNumber, total_reads_num, overall_aligned_num, stringsAsFactors = FALSE)
#just save this somehwere so you can use it, using write.csv(file, "path/filename")
write.csv(bowtie_text_info, "O:/Documents/bowtie_output/bowtie_info.csv")