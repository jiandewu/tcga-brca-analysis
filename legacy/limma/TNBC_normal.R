# source("http://bioconductor.org/biocLite.R")
# biocLite("rgl")
# biocLite("limma")
# biocLite("edgeR")
library(limma)
library("rgl")
library(gplots)
# library(edgeR)
# library(preprocessCore)


setwd("/home/jiande/prjs/BreastTumor/limma")
hdat1 <- read.csv(file = "TNBC110_Normal113_remove.txt", row.names = 1, sep = "\t")
pheno <- read.csv(file = "TNBC110_Normal113_pheno.txt", sep = "\t")

hdat <- apply(hdat1, 2, function(x) log2(x+1))

#hdat <- sweep(hdat1, 2, colSums(hdat1,na.rm=TRUE,dims = 1), "/")*1000

boxplot(hdat)
# 


# boxplot(data_norm)
# colm <- colMeans(hdat)
# diff <- mean(colm[1:110])-mean(colm[111:223])
# hdat[,c(1:110)] <- hdat[,c(1:110)]-diff
# boxplot(hdat)


rownames(pheno) <- pheno$sample
colnames(hdat) <- rownames(pheno)
rownames(pheno) == colnames(hdat)

# data_mat <- data.matrix(hdat) 
# data_norm <- normalize.quantiles(data_mat, copy = TRUE)
# colnames(data_norm) <- rownames(pheno)
# rownames(data_norm) <- rownames(hdat)
# hdat <- as.data.frame(data_norm)

##### PCA
sample.groups <- c(rep("red", 110), rep("blue",113))
pc <- princomp(hdat)
plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)
rgl.postscript("results/TNBC_Normal3d.pdf","pdf")

pc$loadings

samples <- factor(pheno$class)
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)
design

fit <- lmFit(hdat, design)
names(fit)

fit1 <- contrasts.fit(fit, c(-1,1))


fit1 <- eBayes(fit1)
#volcanoplot(fit1,  main = "TNBC vs Normal", xlim=c(-3,3), ylim=c(0,30))
volcanoplot(fit1,  main = "TNBC vs Normal")

nrow(topTable(fit1, coef=1, number=36450, lfc=1, p.value = 0.001))
probeset.list1 <- topTable(fit1, coef=1, number=36450, lfc=1, p.value = 0.001)
head(probeset.list1)
write.table(probeset.list1, file = "results/TNBC_normal_p_fc1.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list1)
head(genes)
sign_exp <- hdat[genes,]
head(sign_exp)
write.table(sign_exp, file = "results/TNBC_normal_exp_fc1.txt", sep="\t", quote = F, row.names=T)

