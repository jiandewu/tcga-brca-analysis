library(limma)
library(Glimma)
library(edgeR)
# source("http://bioconductor.org/biocLite.R")
# biocLite("Homo.sapiens")
# biocLite("biomaRt")
library(biomaRt)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RNAseq123", version = "3.8")

options(stringsAsFactors = FALSE)
setwd("/home/jiande/prjs/BreastTumor/analyse_data/count")

files <- dir(path="TNBC_all.counts", pattern="*\\.txt$", full.names = TRUE)
x <- readDGE(files, columns=c(1,2), header=FALSE)
colnames(x) = gsub("TNBC_all\\.counts/", "", colnames(x))
colnames(x) = gsub("\\.htseq\\.counts", "", colnames(x))
TNBC.samples <- dim(x)[2]

files <- dir(path="nonTNBC_all.counts", pattern="*\\.txt$", full.names = TRUE)
y <- readDGE(files, columns=c(1,2), header=FALSE)
colnames(y) = gsub("nonTNBC_all\\.counts/", "", colnames(y))
colnames(y) = gsub("\\.htseq\\.counts", "", colnames(y))
nonTNBC.samples <- dim(y)[2]

x <- cbind(x,y)

files <- dir(path="TNBC.ctrl.counts", pattern="*\\.txt$", full.names = TRUE)
files <- c(files,dir(path="nonTNBC.ctrl.counts", pattern="*\\.txt$", full.names = TRUE))
y <- readDGE(files, columns=c(1,2), header=FALSE)
colnames(y) = gsub("^.*TNBC\\.ctrl\\.counts/", "", colnames(y))
colnames(y) = gsub("\\.htseq\\.counts", "", colnames(y))
Ctrl.samples <- dim(y)[2]

x <- cbind(x,y)
x <- x[1:(dim(x)[1]-5),]
ENSEMBL <- gsub("(.*?)\\.\\d*?(_\\d)*", "\\1", rownames(x))
row.names(x) <- ENSEMBL

group <- as.factor(rep(c("TNBC","nonTNBC","Ctrl"), c(TNBC.samples,nonTNBC.samples,Ctrl.samples)))
x$samples$group <- group
lane <- as.factor(rep(c("Patient", "Normal"), c(TNBC.samples+nonTNBC.samples,Ctrl.samples)))
x$samples$lane <- lane

geneid <- rownames(x)
genes_ensembl <- read.table("/home/jiande/prjs/RefGene/hg37ref/GRCh37_gene-symbols.txt", sep="\t", header = TRUE)
genes_ensembl <- genes_ensembl[, c(1,2)]
row.names(genes_ensembl) <- genes_ensembl[,1]
genes <- genes_ensembl[geneid,]
genes <- genes[!duplicated(genes$ensembl_gene_id),]

# x$counts <- merge(x$counts,genes,by.x="row.names",by.y="ensembl_gene_id")
x$genes <- genes

cpm <- cpm(x)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
# c(L, M)

table(rowSums(x$counts==0)==ncol(x))

lcpm <- cpm(x, log=TRUE)
# summary(lcpm)
lcpm.cutoff <- log2(10/M + 2/L)

library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
# legend("topright", samplenames, text.col=col, bty="n")


keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
# dim(x)

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.3), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
# legend("topright", samplenames, text.col=col, bty="n")


x <- calcNormFactors(x, method = "TMM")
# x$samples$norm.factors

# x2 <- x
# x2$samples$norm.factors <- 1
# x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
# x2$counts[,2] <- x2$counts[,2]*5
# 
# par(mfrow=c(1,2))
# lcpm <- cpm(x2, log=TRUE)
# boxplot(lcpm, las=2, col=col, main="")
# title(main="A. Example: Unnormalised data",ylab="Log-cpm")
# x2 <- calcNormFactors(x2)  
# x2$samples$norm.factors
# 
# lcpm <- cpm(x2, log=TRUE)
# boxplot(lcpm, las=2, col=col, main="")
# title(main="B. Example: Normalised data",ylab="Log-cpm")


lcpm <- cpm(x, log=TRUE)

# par(mfrow=c(1,2))
# col.group <- group
# levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
# col.group <- as.character(col.group)
# col.lane <- lane
# levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
# col.lane <- as.character(col.lane)
# plotMDS(lcpm, labels=group, col=col.group)
# title(main="A. Sample groups")
# plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
# title(main="B. Sequencing lanes")
# 
# glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), groups=x$samples[,c(2,5)], launch=FALSE)

design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("lane", "", colnames(design))
# design

contr.matrix <- makeContrasts(
  TNBCvsCtrl = TNBC-Ctrl, 
  nonTNBCvsCtrl = nonTNBC-Ctrl, 
  TNBCvsnonTNBC = TNBC-nonTNBC, 
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
# v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

## Significance is defined using an adjusted p-value cutoff that is set at 5% by default.
summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

tefit <- treat(efit, lfc=1)
det <- decideTests(tfit)
summary(det)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$gene_name[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

write.fit(tfit, dt, file="../limmaedgeR/results.txt")

TNBC.vs.Ctrl <- topTreat(tfit, coef=1, n=Inf)
nonTNBC.vs.Ctrl <- topTreat(tfit, coef=2, n=Inf)
head(TNBC.vs.Ctrl)
head(nonTNBC.vs.Ctrl)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)

library(gplots)
TNBC.vs.Ctrl.topgenes <- TNBC.vs.Ctrl$ensembl_gene_id[1:100]
i <- which(v$genes$ensembl_gene_id %in% TNBC.vs.Ctrl.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$gene_name[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

# 
# load("/home/jiande/prjs/RefGene/human_c2_v5p2.rdata")
# idx <- ids2indices(Hs.c2,v$genes$ensembl_gene_id)
# cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
# head(cam.BasalvsLP,5)
