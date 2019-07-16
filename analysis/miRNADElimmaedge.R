library(limma)
library(Glimma)
library(edgeR)
library(data.table)
library(edgeR)
library(stringi)
# source("http://bioconductor.org/biocLite.R")
# biocLite("Homo.sapiens")
# biocLite("biomaRt")
library(Homo.sapiens)
library(biomaRt)

setwd("/home/jiande/prjs/Analysis/BreastCancer/miRNA")
miRNA_normal <- read.csv(file = "Normal104.miRNA.txt", sep = "\t")
miRNA_tnbc <- read.csv(file = "TNBC114.miRNA.txt", sep = "\t")
miRNA_exp <- merge(miRNA_tnbc,miRNA_normal,by = intersect(names(miRNA_tnbc), names(miRNA_normal)))

miRNA_names <- miRNA_exp[,"miRNA_ID"]
miRNA_exp$miRNA_ID <- NULL
rownames(miRNA_exp)<-miRNA_names

miRNA_exp_2<-copy(miRNA_exp)
TNBC.samples <- dim(miRNA_tnbc)[2]-1
Normal.samples <- dim(miRNA_normal)[2]-1
grp_chr <- rep(c("T","N"), c(TNBC.samples,Normal.samples))
grp <- as.factor(grp_chr)
colnames(miRNA_exp_2)<-paste(grp,1:length(grp),sep="")

pairs(log2(1+miRNA_exp_2[,c(1:4,216:218)]), pch=".",lower.panel=NULL)

d <- DGEList(counts=miRNA_exp_2, group=grp)
d <- calcNormFactors(d,method = "TMM")
cols <- as.numeric(d$samples$group)
plotMDS(d,col=cols)

mm <- model.matrix(~-1+grp)
d <- estimateGLMCommonDisp(d,mm)
d <- estimateGLMTrendedDisp(d,mm)
d <- estimateGLMTagwiseDisp(d,mm)
plotBCV(d)
plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
plotSmear(d, pair=c("N","T"), ylim=c(-5,5))

f <- glmFit(d,mm)
con <- makeContrasts("T-N"=grpT-grpN,levels=colnames(mm))
lrt <- glmLRT(f,contrast=con)
topTags(lrt,20)

miRNA_DE_result<-topTags(lrt,Inf)$table
miRNA_DE_result<-cbind(miRNAs=rownames(miRNA_DE_result),miRNA_DE_result)
write.table(miRNA_DE_result,file="./edgeR_miRNA_DE.tsv",quote = F,col.names = T,row.names = F,sep="\t")
nc <- cpm(d,log=T)
nc<-cbind(miRNAs=rownames(nc),nc)
write.table(nc,file="./miRNA_Normalized_logcpm.tsv",quote = F,col.names = T,row.names = F,sep="\t")

d<-dist(t(miRNA_exp_2))
hc1<-hclust(d)
plot(hc1,cex=0.4,hang=-1)

# kmeans聚类
km<-kmeans(t(miRNA_exp_2),2)
km_cluster<-as.data.frame(km$cluster)

# Normal里的分类:
table(km_cluster[which(substr(rownames(km_cluster),1,1)=="N"),])

# Tumor里的分类:
table(km_cluster[which(substr(rownames(km_cluster),1,1)=="T"),])

# library(pheatmap)
# pheatmap(miRNA_exp_2,
#          scale = "none",
#          color = colorRampPalette(colors = c("blue","white","red"))(100)
# )
