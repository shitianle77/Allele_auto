#!Rscript
# Usage: Rscript myscript.R
##################################################################################################################
#Identification and classification of differentially expressed allele pairs
##################################################################################################################

args<-commandArgs(TRUE)

#################################################
# R package
myPaths <- .libPaths()
library(DESeq2)
#################################################
setwd(args[1])

#################################################
#DEseq2
##All
#################################################
#Inputdata
inputdata = read.table(args[2],header=F,row.names=1)
s = args[3]

names(inputdata)
options(stringsAsFactors = FALSE)
data = round(inputdata,digits = 0)  ##Rounding
data = as.matrix(data)
condition = factor(c(rep(c("A"),s),rep(c("B"),s)))
names=colnames(data)
coldata = data.frame(row.names=colnames(data),condition,names)
dds1 <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition ) 
dds <- estimateSizeFactors(dds1)
dds <- DESeq(dds)
AvsB <-  results(dds, contrast=c("condition","A","B"),cooksCutoff=TRUE,alpha = 0.5)
write.table(AvsB,file="AvsB.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##up
AvsB_up = subset(AvsB,padj <= 0.05 & (log2FoldChange > 0))
write.table(AvsB_up,file="AvsB.up.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##down
AvsB_down = subset(AvsB,padj <= 0.05 & (log2FoldChange < 0))
write.table(AvsB_down,file="AvsB.down.txt",sep = "\t",quote = F,row.names = T,col.names = T)

#-------------------
#00-all
#AvsB <- as.data.frame(AvsB[order(AvsB$padj),])
#all = subset(AvsB,padj > 0.05)
#write.table(all,file="diff00.txt",sep = "\t",quote = F,row.names = T,col.names = T)
#-------------------
#0-all
#all = subset(AvsB,padj <= 0.05 & (-1 <= log2FoldChange & log2FoldChange <= 1))
#write.table(all,file="diff0.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##0-up
#all = subset(AvsB,padj <= 0.05 & ((log2FoldChange > 0 & log2FoldChange <= 1)))
#write.table(all,file="diff0.up.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##0-down
#all = subset(AvsB,padj <= 0.05 & ((-1 <= log2FoldChange & log2FoldChange < 0)))
#write.table(all,file="diff0.down.txt",sep = "\t",quote = F,row.names = T,col.names = T)
#------------------
#2-all
#all = subset(AvsB,padj <= 0.05 & ((3 > log2FoldChange & log2FoldChange > 1) | (-3 < log2FoldChange & log2FoldChange < -1)))
#write.table(all,file="diff2.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##2-up
#all = subset(AvsB,padj <= 0.05 & ((3 > log2FoldChange & log2FoldChange > 1)))
#write.table(all,file="diff2.up.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##2-down
#all = subset(AvsB,padj <= 0.05 & ((-3 < log2FoldChange & log2FoldChange < -1)))
#write.table(all,file="diff2.down.txt",sep = "\t",quote = F,row.names = T,col.names = T)
#-------------------
#8-all
#all = subset(AvsB,padj <= 0.05 & (log2FoldChange >= 3 | log2FoldChange <= -3))
#write.table(all,file="diff8.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##8-up
#all = subset(AvsB,padj <= 0.05 & (log2FoldChange >= 3))
#write.table(all,file="diff8.up.txt",sep = "\t",quote = F,row.names = T,col.names = T)
##8-down
#all = subset(AvsB,padj <= 0.05 & (log2FoldChange <= -3))
#write.table(all,file="diff8.down.txt",sep = "\t",quote = F,row.names = T,col.names = T)

#################################################
#DEseq2
##Each treatments
#################################################
#Inputdata
inputdata = read.table(args[2],header=F,row.names=1)

l <- ncol(inputdata)
g <- seq(from=1,by=3,length=ncol(inputdata)/6)

#################################################
n <- 1
for(i in g){
a <- i
b <- i+2
c <- l/2+i
d <- l/2+i+2
  if(n == 1){
print(c(a,b,c,d))
data = inputdata[,c(a:b,c:d)]
  }
#head(data)
names(data)
options(stringsAsFactors = FALSE)
data = round(data,digits = 0)
data = as.matrix(data)
condition = factor(c(rep(c("A"),3),rep(c("B"),3)))
names=colnames(data)
coldata = data.frame(row.names=colnames(data),condition,names)
dds1 <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition )
dds <- estimateSizeFactors(dds1)
dds <- DESeq(dds)
AvsB <-  results(dds, contrast=c("condition","A","B"),cooksCutoff=F,alpha = 0.5)
AvsB <- as.data.frame(AvsB[order(AvsB$padj),])
all = subset(AvsB,padj > 0.05 )
write.table(all,paste("a", i, ".diff00.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & (-1 <= log2FoldChange & log2FoldChange <= 1))
write.table(all,paste("a", i, ".diff0.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & ((3 > log2FoldChange & log2FoldChange > 1) | (-3 < log2FoldChange & log2FoldChange < -1)))
write.table(all,paste("a", i, ".diff2.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & (log2FoldChange >= 3 | log2FoldChange <= -3))
write.table(all,paste("a", i, ".diff8.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
#up
all = subset(AvsB,padj <= 0.05 & (log2FoldChange > 0 & log2FoldChange <= 1))
write.table(all,paste("a", i, ".diff0.sup.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & ((-3 < log2FoldChange & log2FoldChange < -1)))
write.table(all,paste("a", i, ".diff2.sup.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & (log2FoldChange <= -3))
write.table(all,paste("a", i, ".diff8.sup.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
#down
all = subset(AvsB,padj <= 0.05 & ((-1 <= log2FoldChange & log2FoldChange < 0)))
write.table(all,paste("a", i, ".diff0.dom.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & ((3 > log2FoldChange & log2FoldChange > 1)))
write.table(all,paste("a", i, ".diff2.dom.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
all = subset(AvsB,padj <= 0.05 & (log2FoldChange >= 3))
write.table(all,paste("a", i, ".diff8.dom.txt", sep = ""), sep = "\t",quote = F,row.names = T,col.names = T)
}

