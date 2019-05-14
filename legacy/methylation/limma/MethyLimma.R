# source("http://bioconductor.org/biocLite.R")
# biocLite("rgl")
# biocLite("limma")
library(limma)
library("rgl")
library(gplots)
library(preprocessCore)

setwd("~/prjs/BreastTumor/methylation/limma")
hdat1 <- read.csv(file = "normal_TNBC.txt", row.names = 1, sep = "\t")
pheno <- read.csv(file = "normal_TNBC_pheno.txt", sep = "\t")

# hdat=removeBatchEffect(hdat,batch=c(rep("0",96),rep("1",80)))
# hdat <- hdat1/colSums(hdat1,na.rm=TRUE,dims = 1)*200000
# hdat <- sweep(hdat1, 2, colSums(hdat1,na.rm=TRUE,dims = 1), "/")*200000
data_mat <- data.matrix(hdat1)
data_norm <- normalize.quantiles(data_mat, copy = TRUE)
colnames(data_norm) <- colnames(hdat1)
rownames(data_norm) <- rownames(hdat1)
hdat <- as.data.frame(data_norm)

boxplot(hdat)


rownames(pheno) <- pheno$sample
colnames(hdat) <- rownames(pheno)
rownames(pheno) == colnames(hdat)

##### PCA
sample.groups <- c(rep("red", 96), rep("blue",79))
pc <- princomp(hdat)
plot(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", main = "PCA")
plot3d(pc$loadings, col=sample.groups, xlab= "Component 1", ylab="Component 2", zlab="Component 3", type = "s", size = 1.5)

pc$loadings

samples <- factor(pheno$class)
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)
design

fit <- lmFit(hdat, design)
names(fit)

fit1 <- contrasts.fit(fit, c(-1,1))


fit1 <- eBayes(fit1)
volcanoplot(fit1,  main = "Normal VS TNBC")


write.table(topTable(fit1, number=dim(fit1)[1]), file = "result/normal_TNBC_p.all.txt", sep="\t", quote = F, row.names=T)
write.table(hdat, file = "result/normal_TNBC_exp.all.txt", sep="\t", quote = F, row.names=T)


nrow(topTable(fit1, coef=1, number=394400, lfc=0.25, p.value = 0.001))
probeset.list1 <- topTable(fit1, coef=1, number=394400, lfc=0.25, p.value = 0.001)
# probeset.list2 <- topTable(fit1, coef=1, number=394399, lfc=2, p.value = 0.05)
head(probeset.list1)

write.table(probeset.list1, file = "result/normal_TNBC_p.txt", sep="\t", quote = F, row.names=T)
# write.table(probeset.list2, file = "result/BLIA_BLIS_p_lfc2.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list1)
head(genes)
sign_exp <- hdat[genes,]
head(sign_exp)
write.table(sign_exp, file = "result/normal_TNBC_exp.txt", sep="\t", quote = F, row.names=T)
# 
# genes <- rownames(probeset.list2)
# head(genes)
# sign_exp <- hdat[genes,]
# head(sign_exp)
# write.table(sign_exp, file = "result/BLIA_BLIS_exp_lfc2.txt", sep="\t", quote = F, row.names=T)

