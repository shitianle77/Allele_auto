#!Rscript
# Usage: Rscript myscript.R
##################################################################################################################
# Comparison on number of highly expressed genes from different homologous chromosomes
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
df <- read.table(args[2],header = T)

#################################################
#run
mean <- aggregate(df$Number,by=list(df$Group,df$Chromosomes),FUN=mean)
sd <- aggregate(df$Number,by=list(df$Group,df$Chromosomes),FUN=sd)
N <- aggregate(df$Number,by=list(df$Group,df$Chromosomes),FUN=length)

data <- data.frame(mean,sd=sd$x,N=N$x)
colnames(data)=c("Group","Chromosomes","Number","sd","N")
data$se <- data$sd / sqrt(data$N)

p <- ggplot(data, aes(x=Chromosomes, y=Number, fill=Group))+
  geom_bar(position=position_dodge(), 
           color="black",
           stat="identity",
           width=.6)+
  geom_errorbar(aes(ymin=Number-se,
                    ymax=Number+se),
                width=.2, 
                size=0.6,
                position=position_dodge(.6))+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        text=element_text(size=14, family="serif"))+
  scale_fill_manual(values = c("#84C3B9","#A976AA"))

pdf("Highexpression_chr.pdf",width = 8, height = 6)
p
dev.off()


