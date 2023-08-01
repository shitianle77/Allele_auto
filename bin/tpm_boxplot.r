#!Rscript
# Usage: Rscript myscript.R
##################################################################################################################
#Box plotting
##################################################################################################################

args<-commandArgs(TRUE)

#################################################
# R package
myPaths <- .libPaths()
library("ggplot2")

#################################################
setwd(args[1])

#################################################
#Inputdata
data = read.table(args[2],header=F)
colnames(data) <- c("TPM","Group")

#################################################
#run
p <- ggplot(data, aes(x=Group,y=TPM,fill=Group)) + 
  geom_boxplot(outlier.colour="white", outlier.size  = 0)+
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0,100))+ ##The range of the Y-axis is adjusted according to your own data.
  theme_bw() +theme(panel.grid=element_blank())+
  theme(axis.title.x = element_text(face="plain", size=10),axis.text.x = element_text(vjust=1, size=10)) + theme(axis.title.y = element_text(face="plain", size=10),axis.text.y  = element_text(vjust=1, size=10)) +
  theme(legend.title = element_text(face = "plain",size=9),legend.text = element_text(face="plain",size = 9)) +
  theme(text=element_text(size=11,  family="serif")) +
  scale_fill_manual(values = c("#6DA8C9","#32709E","#294F7C"))
# p
pdf("tpm_boxplot.pdf",width = 3, height = 4)
p
dev.off()

data = read.table(args[2],header=F)
head(data)
data1 <- data[data[,2]=="diff0",]
#head(data1)
data2 <- data[data[,2]=="diff2",]
#head(data2)
data3 <- data[data[,2]=="diff8",]
#head(data3)

wilcox.test(data1$V1,data2$V1)
wilcox.test(data2$V1,data3$V1)
wilcox.test(data1$V1,data3$V1)
