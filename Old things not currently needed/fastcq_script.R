#Looking at sequencing data
#First step will be quality control
#I'll try to run fastqc from R

#this won't work on a windows machine because my life is too hard
#so let's get it ready to run on the lab Mac instead
#that is too annoying because I have to install a bunch of stuff
#I'll try to call the .bat file from R? Can I do that?

library(tidyverse)
library(fastqcr)
datapath="X:\\fast\\bloom_j\\SR\\ngs\\illumina\\sdresche\\190503_M03100_0427_000000000-CG8MD\\Data\\Intensities\\BaseCalls"
resultspath="O:\\JO_LabShard\\Sara Drescher\\Molecular_clock\\Analysis\\fastqc"

fastqc(fq.dir = datapath, # FASTQ files director
       qc.dir = resultspath, # Results direcory
       threads = 4                    # Number of threads
)


#files are named like this: 5_S5_L001_I1_001.fastq.gz 
#there are 94 samples so the last is 94_S94_L001_**_001.fastq.gz
#the ** is either I1, I2, R1, or R2
#I only care about R1 and R2

#so.....*_S*_L001_R*_001.fastq.gz
#this gets all the names
myseqfiles=dir(datapath, pattern= ".*_S.*_L001_R[1-2]_001.fastq.gz")

# make them one string
myseqfiles_path=paste(datapath, myseqfiles, sep="\\")
myseqstring=paste(myseqfiles_path, collapse=" ")


library (readr)
library (stringr)

#bat_code=str_remove(read_file("c:\\users\\sara\\downloads\\fastqc_v0.11.8\\FastQC\\fastqc.bat"), "%*")
#that's really not necessary

bat_code="java -Xmx250m -Dfastqc.output_dir=C:\\Users\\Sara\\Documents\\HHMI_projects\
\Molecular_clock\\Analysis\\fastqc -classpath .;./sam-1.103.jar;./jbzip2-0.9.jar 
uk.ac.babraham.FastQC.FastQCApplication %*"

new_bat_code=paste("java -Xmx250m -Dfastqc.output_dir=", resultspath," -classpath .;./sam-1.103.jar;./jbzip2-0.9.jar  uk.ac.babraham.FastQC.FastQCApplication", myseqstring, sep=" ")

#pasted "new_bat_code" into a text file named fastqnew.bat
#ran fastqnew at the cmd
###new_bat_code=str_remove(string, pattern)
###shell.exec ("c:\\users\\sara\\downloads\\fastqc_v0.11.8\\FastQC\\fastqcA.bat")


qc.dir="C:\\Users\\Sara\\Documents\\HHMI_projects\\Molecular_clock\\Analysis\\fastqc"
qc <- qc_aggregate(qc.dir)


