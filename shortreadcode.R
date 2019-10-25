#trying to play with shortread library to look at our files

library(ShortRead)


fls <- dir("C:\\Users\\Sara\\Documents\\HHMI_projects\\Molecular_clock\\Analysis", "*fastq.gz$", full=TRUE)
#collecting statistics over the files
qaSummary <- qa(fls, type="fastq")
#and creating and viewing a report
browseURL(report(qaSummary))
#By default, the report is based on a sample of 1M reads.
#These QA facilities are easily augmented by writing custom f

writeFasta(fls)
writeFasta(readFastq("C:\\Users\\Sara\\Documents\\HHMI_projects\\Molecular_clock\\Analysis", "1_S1_L001_R2_001.fastq.gz"), "1_S1_L001_R2_001.fastq.gz")