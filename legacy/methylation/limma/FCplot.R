library(limma)

setwd("/home/jiande/prjs/BreastTumor/methylation/limma/batch/")

list1 <- read.csv(file = "Methy_RNAseq_all_FC.txt", row.names = 1, sep = "\t")

plot(list1$Methy.logFC, list1$RNAseq.logFC, main="Scatterplot Example", 
     xlab="Methylation logFC ", ylab="RNAseq logFc", pch=19)
