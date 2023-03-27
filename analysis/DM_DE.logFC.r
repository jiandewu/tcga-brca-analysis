# install.packages("readr")
library(readr)
data1 <- read_delim("/home/jiande/prjs/Analysis/BreastCancer/paper_methy/DM_DE.logFC.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

View(data1)

data1$Colour="black"
data1$Colour[data1$Methylation <= -1]="green"
data1$Colour[data1$Methylation >= 1]="orange"
data1$Colour[data1$RNAseq <= -2]="aquamarine1"
data1$Colour[data1$RNAseq >= 2]="red"



data1$Colour[data1$Methylation >= 1 & data1$RNAseq >= 2]="brown"
data1$Colour[data1$Methylation <= -1 & data1$RNAseq >= 2]="blue"
data1$Colour[data1$Methylation >= 1 & data1$RNAseq <= -2]="blue"
data1$Colour[data1$Methylation <= -1 & data1$RNAseq <= -2]="brown"




plot(data1$Methylation, data1$RNAseq, col=data1$Colour , pch=19, cex=.3, xlab = "", ylab = "")

title("TNBC vs Normal", cex.main = 2, col.main= "black",
      xlab="Methylation (LogFC)", ylab="RNAseq (LogFC)", cex.lab=1.5)

abline(v= -1, col="red")
abline(v= 1, col="red")
abline(h= 2, col="red")
abline(h= -2, col="red")


legend1 <- c(rep("red", 1), rep("brown", 1), rep("green", 1), rep("aquamarine1", 1), rep("blue",1), rep("orange", 1))
legend2 <- c(paste("Up-regulated (",sum(data1$Colour=="red"),")",sep = ""),
             paste("Hypo-methyl, Down-regulated +\nHyper-methyl,Up-regulated (",sum(data1$Colour=="brown"),")",sep = ""),
             paste("Hypo-methylated (",sum(data1$Colour=="green"),")",sep = ""),
             paste("Down-regulated (",sum(data1$Colour=="aquamarine1"),")",sep = ""),
             paste("Hypo-methyl, Up-regulated +\nHyper-methyl,Down-regulated (",sum(data1$Colour=="blue"),")",sep = ""),
             paste("Hyper-methylated (",sum(data1$Colour=="orange"),")",sep = ""))
legend("bottomright", inset=.02, title="Legend",
      legend=legend2, fill=legend1, cex=1.5)

# colors()
