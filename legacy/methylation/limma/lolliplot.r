# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("trackViewer", version = "3.8")
library(trackViewer)
library(readxl)

# gene <- read_excel("/home/jiande/nas/Tarun_projects/prostate_cancer_new/All_info/Plot/KLK2.xlsx")
setwd("~/prjs/BreastTumor/methylation/limma/result")
hdat <- read.csv(file = "all_intersection_position.txt", sep = "\t")

genelist <- unique(hdat$Gene)
genespec <- 4

# gene <- hdat[hdat$Gene=="COL5A2",c(2,4)]
gene <- hdat[hdat$Gene==genelist[genespec],c(2,4)]

# sample.groups <- c(rep("red", 4), rep("green", 5), rep("blue",1))
sample.groups <- as.character(gene$Mutation)
sample.groups <- gsub("Germline", "red", sample.groups)
sample.groups <- gsub("Somatic", "blue", sample.groups)
sample.groups <- gsub("Methylation", "green", sample.groups)

chromo <- as.character(hdat[hdat$Gene==genelist[genespec],3])[1]
start <- min(gene$Start)
end <- max(gene$Start)

# KLK2 <- GRanges("chr17", IRanges(gene$Position, width=1))
gene_name <- GRanges(chromo, IRanges(gene$Start, width=1))
# KLK2$color <- sample.groups
gene_name$color <- sample.groups
# KLK2$score <- sample.int(20, length(KLK2), replace=TRUE)
gene_name$score <- sample.int(20, length(gene_name), replace=TRUE)
# features <- GRanges("chr17:51364600- 51393200")
features <- GRanges(paste(chromo,":",start,"- ",end,sep=""))

# legend <- c(rep("red", 1), rep("blue",1), rep("green", 1))
legend <- gsub("Methylation", "green", gsub("Somatic", "blue", gsub("Germline", "red", as.character(unique(gene$Mutation)))))
# names(legend) <- c("Germline","Somatic", "Methylation")
names(legend) <- as.character(unique(gene$Mutation))

# lolliplot(KLK2, features, ranges = GRanges("chr17", IRanges(51364600, 51393200))
#           , yaxis=FALSE, legend=legend)
lolliplot(gene_name, features, ranges = GRanges(chromo, IRanges(start-2, end+2))
          , yaxis=FALSE, ylab=as.character(genelist[genespec]), rescale=TRUE, legend=legend)

