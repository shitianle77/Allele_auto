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
data = read.table(args[2],header=F,sep = "\t")
colnames(data) <- c("Allelepair","Ka","Ks","Ka/Ks","Group")

#################################################
#run
##Ka
p1 <- ggplot(data, aes(x=Group,y=Ka,fill=Group)) + 
  geom_boxplot(outlier.colour="white", outlier.size  = 0)+
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0,0.15))+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(axis.title.x = element_text(face="plain", size=10),axis.text.x = element_text(vjust=1, size=10)) + theme(axis.title.y = element_text(face="plain", size=10),axis.text.y  = element_text(vjust=1, size=10)) +
  theme(legend.title = element_text(face = "plain",size=9),legend.text = element_text(face="plain",size = 9)) +
  theme(text=element_text(size=11,  family="serif"))+
  scale_fill_manual(values = c("#6DA8C9","#32709E","#294F7C"))

pdf("ka.pdf", width = 2.5, height = 3)
p1
dev.off()

##Ks
p2 <- ggplot(data, aes(x=Group,y=Ks,fill=Group)) + 
  geom_boxplot(outlier.colour="white", outlier.size  = 0)+
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0,0.22))+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(axis.title.x = element_text(face="plain", size=10),axis.text.x = element_text(vjust=1, size=10)) + theme(axis.title.y = element_text(face="plain", size=10),axis.text.y  = element_text(vjust=1, size=10)) +
  theme(legend.title = element_text(face = "plain",size=9),legend.text = element_text(face="plain",size = 9)) +
  theme(text=element_text(size=11,  family="serif"))+
  scale_fill_manual(values = c("#6DA8C9","#32709E","#294F7C"))

pdf("ks.pdf", width = 2.5, height = 3)
p2
dev.off()

##KaKs
p3 <- ggplot(data, aes(x=Group,y=Ka/Ks,fill=Group)) + 
  geom_boxplot(outlier.colour="white", outlier.size  = 0)+
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0,1.85))+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(axis.title.x = element_text(face="plain", size=10),axis.text.x = element_text(vjust=1, size=10)) + theme(axis.title.y = element_text(face="plain", size=10),axis.text.y  = element_text(vjust=1, size=10)) +
  theme(legend.title = element_text(face = "plain",size=9),legend.text = element_text(face="plain",size = 9)) +
  theme(text=element_text(size=11,  family="serif")) +
  scale_fill_manual(values = c("#6DA8C9","#32709E","#294F7C"))

pdf("kaks.pdf", width = 2.5, height = 3)
p3
dev.off()

##Mann–Whitney–Wilcoxon test
data = read.table(args[2],header=F,sep = "\t")
head(data)
SA <- data[data[,5]=="diff0",]
SB <- data[data[,5]=="diff2",]
SC <- data[data[,5]=="diff8",]
wilcox.test(SA$V2,SB$V2)
wilcox.test(SB$V2,SC$V2)
wilcox.test(SA$V2,SC$V2)

wilcox.test(SA$V3,SB$V3)
wilcox.test(SB$V3,SC$V3)
wilcox.test(SA$V3,SC$V3)

wilcox.test(SA$V4,SB$V4)
wilcox.test(SB$V4,SC$V4)
wilcox.test(SA$V4,SC$V4)
<<<<<<< HEAD
=======






















>>>>>>> 8cc22e0 (all updata)
