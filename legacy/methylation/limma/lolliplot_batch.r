# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("trackViewer", version = "3.8")
library(trackViewer)
library(readxl)

# gene <- read_excel("/home/jiande/nas/Tarun_projects/prostate_cancer_new/All_info/Plot/KLK2.xlsx")
setwd("~/prjs/BreastTumor/methylation/limma/result")
hdat <- read.csv(file = "all_intersection_position_geneloc.txt", sep = "\t")

genelist <- unique(hdat$Gene)

gene_name <- list()
features <- list()
# legend <- list()
legended <- c("red", "blue", "green")
names(legended) <- c("Germline","Somatic","Methylation")
legendlist <- list()

plotlist <- 1:length(genelist)
# plotlist <- c(2)

# genespec <- 2
for (genespec in plotlist) {

# gene <- hdat[hdat$Gene=="COL5A2",c(2,4)]
gene_syb <- as.character(genelist[genespec])
gene <- hdat[hdat$Gene==gene_syb,c(2,4)]

chromo <- as.character(hdat[hdat$Gene==genelist[genespec],3])[1]
loli_start <- min(gene$Position)
loli_end <- max(gene$Position)
gene_start <- hdat[hdat$Gene==genelist[genespec],5][1]
gene_end <- hdat[hdat$Gene==genelist[genespec],6][1]
start <- min(loli_start,gene_start)
end <- max(loli_end,gene_end)
sample.groups <- as.character(gene$Mutation)
legendlist <- c(legendlist,sample.groups)

# sample.groups <- c(rep("red", 4), rep("green", 5), rep("blue",1))
sample.groups <- gsub("Germline", "red", sample.groups)
sample.groups <- gsub("Somatic", "blue", sample.groups)
sample.groups <- gsub("Methylation", "green", sample.groups)

# KLK2 <- GRanges("chr17", IRanges(gene$Position, width=1))
gene_name[[gene_syb]] <- GRanges(chromo, IRanges(gene$Position, width=1))
# KLK2$color <- sample.groups
gene_name[[gene_syb]]$color <- sample.groups
# KLK2$score <- sample.int(20, length(KLK2), replace=TRUE)
# gene_name[[gene_syb]]$score <- sample.int(20, length(gene_name[[gene_syb]]), replace=TRUE)
gene_name[[gene_syb]]$score <- rep(5.2,length(gene_name[[gene_syb]]))
# features <- GRanges("chr17:51364600- 51393200")
# features[[gene_syb]] <- GRanges(paste(chromo,":",start,"- ",end,sep=""))
features[[gene_syb]] <- GRanges(chromo,IRanges(c(start,gene_start,gene_end),width=c(gene_start-start, gene_end-gene_start, end-gene_end)))
features[[gene_syb]]$fill <- c("#FFFFFF", "#51C6E6", "#FFFFFF")

# legend <- c(rep("red", 1), rep("blue",1), rep("green", 1))
# legend[[gene_syb]] <- gsub("Methylation", "green", gsub("Somatic", "blue", gsub("Germline", "red", as.character(unique(gene$Mutation)))))
# names(legend) <- c("Germline","Somatic", "Methylation")
# names(legend[[gene_syb]]) <- as.character(unique(gene$Mutation))
}

# lolliplot(KLK2, features, ranges = GRanges("chr17", IRanges(51364600, 51393200))
#           , yaxis=FALSE, legend=legend)
legendlist <- unique(legendlist)
png(filename="lolli.png",width=800,height=125*length(plotlist)+40)
lolliplot(gene_name, features, yaxis=FALSE, rescale=TRUE, legend=legended[unlist(legendlist)])
dev.off()
