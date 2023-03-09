#Rscript
# Usage: Rscript myscript.R
##################################################################################################################
# Grouping of allele-specific gene (ASE) expression profile probably among samples from different tissues and treatments
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
data <- read.table(args[2],header=F)
list <- read.table(args[3],header=F)
list <- rev(as.character(list$V1))

#################################################
#run
data$V2 <- factor(data$V2,levels=list)
data$V3 <- factor(data$V3, levels=c('diff8','diff2','diff0','diff00','Noexpression'))
colnames(data) <- c("Percentage","Treatment","Group")
#head(data)

cols <-c("#294F7C","#32709E","#6DA8C9","#BFD3E6","#E5E1E4")	
p <- ggplot(data,aes(Treatment,Percentage,fill=Group))+
  geom_bar(stat="identity",position="fill")+ 
  theme(legend.position="none") +
  theme_bw() +theme(panel.grid=element_blank()) +
  scale_fill_manual(values = cols) +
  coord_flip()
pdf("duidietu.pdf",width = 5, height = 8)
p
dev.off()

