#!/usr/bin/Rscript

# source("http://bioconductor.org/biocLite.R")
# biocLite("rgl")
# biocLite("limma")
# biocLite("edgeR")
library(limma)
# library("rgl")
library(gplots)
library(edgeR)
# library(preprocessCore)

setwd("/home/jiande/prjs/BreastTumor/analyse_data/")
treat <- read.csv(file = "TNBC.Tumor.counts.txt", row.names = 1, sep = "\t")
ctrl <- read.csv(file = "TNBC.ctrl.counts.txt", row.names = 1, sep = "\t")

names(treat) <- paste("treat.", names(treat), sep = "")
names(ctrl) <- paste("ctrl.", names(ctrl), sep = "")

pheno <- cbind(c(names(treat), names(ctrl)), c(rep("Treat",length(treat)), rep("Ctrl",length(ctrl))))
colnames(pheno) <- c("sample", "class")
pheno <- as.data.frame(pheno)

data <- cbind(treat,ctrl)

# boxplot(data)

# DESeq2
calc_sf <- function (expr_mat, spikes=NULL){
  geomeans <- exp(rowMeans(log(expr_mat)))
  SF <- function(cnts){
    median((cnts/geomeans)[(is.finite(geomeans) & geomeans >0)])
  }
  norm_factor <- apply(expr_mat,2,SF)
  return(t(t(expr_mat)/norm_factor))
}

data_norm <- data[!rowSums(data==0) > ncol(data)*.3,]
data_norm <- calc_sf(data_norm)
# boxplot(data_norm)

#  TMM
# f <- calcNormFactors(data_norm, method = "TMM")
# data_norm <- t(t(data_norm) /f )


data_norm <- log2(data_norm+1)

rownames(pheno) <- pheno$sample
colnames(data_norm) <- pheno$sample
# rownames(pheno) == colnames(data_norm)

# data_mat <- data.matrix(hdat) 
# data_norm <- normalize.quantiles(data_mat, copy = TRUE)
# colnames(data_norm) <- rownames(pheno)
# rownames(data_norm) <- rownames(hdat)
# hdat <- as.data.frame(data_norm)

##### PCA
# sample.groups <- c(rep("red", length(treat)), rep("blue",length(ctrl)))
# pc <- princomp(data_norm)
# plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
# plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)
# rgl.postscript("results/TNBC_Normal3d.pdf","pdf")
# 
# pc$loadings

samples <- factor(pheno$class)
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)
# design

fit <- lmFit(data_norm, design)
# names(fit)

fit1 <- contrasts.fit(fit, c(-1,1))


fit1 <- eBayes(fit1)
#volcanoplot(fit1,  main = "TNBC vs Normal", xlim=c(-3,3), ylim=c(0,30))
volcanoplot(fit1, main = "TNBC vs Ctrl.", highlight = 0)

# nrow(topTable(fit1, coef=1, number=27767, lfc=2, p.value = 0.001))
probeset.list1 <- topTable(fit1, coef=1, number=30000, lfc=2, p.value = 0.001)
# head(probeset.list1)
write.table(probeset.list1, file = "results/TNBC_Ctrl_p_fc2.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list1)
# head(genes)
sign_exp <- data_norm[genes,]
# head(sign_exp)
write.table(sign_exp, file = "results/TNBC_Ctrl_exp_fc2.txt", sep="\t", quote = F, row.names=T)

