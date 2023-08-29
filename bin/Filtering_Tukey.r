#!Rscript Usage: Rscript myscript.R
##################################################################################################################
#Outlier filtering
##################################################################################################################

args<-commandArgs(TRUE)

#################################################
# R package
myPaths <- .libPaths()
library(ggplot2)
library(dplyr)
#################################################
setwd(args[1])

#################################################
#Tukey method
#################################################
#Inputdata
inputdata = read.table(args[2],header=T) 
s = as.numeric(args[3])
f = rlang::sym(args[4])

find_outlier <- function(x) {
  return(x < quantile(x, .25) - s*IQR(x) | x > quantile(x, .75) + s*IQR(x))
}

df <- inputdata %>%
  mutate(outlier = ifelse(find_outlier(!!f), NA, !!f))

write.table(df,paste(rlang::as_string(f),"-", s, ".genepairs_info.tsv", sep = ""), sep = "\t",quote = F,row.names = F,col.names = T)
