library(ggplot2)
library(ggpubr)
#make some plots for lab meeting.
filtered1=read.csv("C:\\Users\\Sara\\Documents\\HHMI_projects\\Molecular_clock\\Analysis\\third_output_nano_filtered_1.csv")
filtered1$Fragment=filtered1$V6

filtered1$Month=(strsplit(as.character(filtered1$Name), split="_"))

filtered1%<>%separate( Name, c("Pt", "Sample_name", "Month", "Fragement", "Replicate"), sep="_", remove=FALSE)


filtered2=read.csv("C:\\Users\\Sara\\Documents\\HHMI_projects\\Molecular_clock\\Analysis\\third_output_second_run_10.csv")


filtered_noVA1=filtered1[(grep("M", filtered1$Month)),]
filtered_noVA1$Month=as.numeric(sub("M", "", filtered_noVA1$Month))
#filtered_noVA1%<>%rename("APD"="V2")

filtered_noVA1$Sample_name[15]="pt_116_"
filtered_noVA1%<>%mutate(Experiment="1")
#filtered_noVA1=separate(filtered_noVA1, Name, into="Sample_name", sep = "_", remove = F,
         #convert = FALSE, extra = "drop")


filtered_noVA2=filtered2[which(!is.na(filtered2$Month)),]
#filtered_noVA$Sample_name[15]="pt_116_"

filtered_noVA2=separate(filtered_noVA2, Name, into="Sample_name", sep = "_", remove = F,
                       convert = FALSE, extra = "drop")
filtered_noVA2[which(filtered_noVA2$Sample_name%in%c("pt93", "pt93gelcleaned")), "Sample_name"]="pt_93_"
#filtered_noVA2=separate(filtered_noVA2, Fragment, into=c("nothing", "Fragment"), sep="F")
filtered_noVA2%<>%mutate(Experiment="2")


filtered_all=rbind(filtered_noVA1[,c("Sample_name", "APD", "Month", "Fragment", "Experiment", "GoodAAs")], filtered_noVA2[,c("Sample_name", "APD", "Month", "Fragment", "Experiment", "Good.AAs")])

pdf("O:/JOLabShared/SaraDrescher/Molecular_clock/some_plots.pdf")
ggplot(filtered_noVA1, aes(y=APD._1, x=Month)) + geom_point() +  geom_smooth(method = "lm") +  stat_regline_equation(label.x = 3, label.y = .015) + ggtitle ("Exp 1--filtered")

ggplot() + geom_point(data=filtered_noVA1, aes(y=APD._1, x=Month), color="red") +  geom_point(data=filtered_noVA1, aes(y=APD_5, x=Month), color="green")  + geom_point(data=filtered_noVA1, aes(y=APD_10, x=Month), color="blue") + geom_smooth(method = "lm")



ggplot(filtered_noVA2, aes(y=APD, x=Month)) + geom_point() +  geom_smooth(method = "lm") +  stat_regline_equation(label.x = 3, label.y = .015) + ggtitle ("Exp 2--filtered")
ggplot(filtered_all, aes(y=APD, x=Month)) + geom_point() +  geom_smooth(method = "lm") +  stat_regline_equation(label.x = 3, label.y = .015) + ggtitle ("All--filtered 10%")


ggplot(filtered_noVA1, aes(y=APD, x=Month, color=Fragment, fill=Fragment)) + geom_point(aes(color=Fragment)) + geom_smooth(method = "lm", se=T) + stat_regline_equation(label.x = 3, label.y = .025) + ggtitle ("Exp 1--filtered") + facet_grid(col=vars(Fragment))
ggplot(filtered_noVA2, aes(y=APD, x=Month, color=Fragment, fill=Fragment)) + geom_point(aes(color=Fragment)) + geom_smooth(method = "lm", se=T) + stat_regline_equation(label.x = 3, label.y = .025) + ggtitle ("Exp 2--filtered") + facet_grid(col=vars(Fragment))
ggplot(filtered_all, aes(y=APD, x=Month, color=Fragment, fill=Fragment)) + geom_point(aes(color=Fragment)) + geom_smooth(method = "lm", se=T) + stat_regline_equation(label.x = 3, label.y = .025) + ggtitle ("All--filtered 10%") + facet_grid(col=vars(Fragment))


ggplot(filtered_all, aes(y=APD, x=Month, color=Experiment, fill=Experiment)) + geom_point(aes(color=Experiment)) + geom_smooth(method = "lm", se=T) + ggtitle ("All--filtered 10%") + facet_grid(col=vars(Fragment))



ggplot(filtered_noVA1, aes(y=APD, x=Month, color=Sample_name, fill=Sample_name)) + geom_point(aes(color=Sample_name)) + geom_smooth(method = "lm", se=F) + ggtitle ("Exp 1, by Fragment") + facet_grid(col=vars(Fragment))
ggplot(filtered_noVA2, aes(y=APD, x=Month, color=Sample_name, fill=Sample_name)) + geom_point(aes(color=Sample_name)) + geom_smooth(method = "lm", se=F) + ggtitle ("Exp 2, by Fragment") + facet_grid(col=vars(Fragment))
ggplot(filtered_all, aes(y=APD, x=Month, color=Sample_name, fill=Sample_name)) + geom_point(aes(color=Sample_name)) + geom_smooth(method = "lm", se=F) + ggtitle ("All, by Fragment 10%") + facet_grid(col=vars(Fragment))

dev.off()


filtered_keep=filtered_all[-which(is.nan(filtered_all$APD)),]
glm.fit=glm(APD~Month, data=filtered_keep)

#LOOCV
cv.glm(filtered_keep,glm.fit)$delta 


lm(APD~Month, filtered_noVA[which(filtered_noVA$Fragment=="1_"),])
lm(APD~Month, filtered_noVA[which(filtered_noVA$Fragment=="2_"),])
lm(APD~Month, filtered_noVA[which(filtered_noVA$Fragment=="3_"),])
lm(APD~Month, filtered_noVA)