#Rscript
# Usage: Rscript myscript.R
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

#################################################
#run
p <- ggplot(data, aes(x = (V2+V3)/2/1000000, y = (V4+V5)/2/1000000)) +
  geom_point(size = .75) +
  facet_wrap(~ V1, scales = "free") +
  labs(x = "ChromosomeA (Mb)", y = "ChromosomeB (Mb)") +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line(color = "black"))

pdf("pairs.allele.pdf") 
p
dev.off()  

