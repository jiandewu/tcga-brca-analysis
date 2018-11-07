# source("http://bioconductor.org/biocLite.R")
# biocLite("rgl")
# biocLite("limma")
library(limma)
library("rgl")
library(gplots)


setwd("~/prjs/BreastTumor/HicksRaw")
hdat <- read.csv(file = "BLIA_BLIS.txt", row.names = 1, sep = "\t")
pheno <- read.csv(file = "BLIA_BLIS_pheno.txt", sep = "\t")

# hdat=removeBatchEffect(hdat,batch=c("0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"
# ))
boxplot(hdat)


rownames(pheno) <- pheno$sample
colnames(hdat) <- rownames(pheno)
rownames(pheno) == colnames(hdat)

##### PCA
sample.groups <- c(rep("red", 54), rep("blue",60))
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
volcanoplot(fit1,  main = "BLIA VS BLIS")

nrow(topTable(fit1, coef=1, number=54675, lfc=1, p.value = 0.05))
probeset.list1 <- topTable(fit1, coef=1, number=54675, lfc=1, p.value = 0.05)
probeset.list2 <- topTable(fit1, coef=1, number=54675, lfc=2, p.value = 0.05)
head(probeset.list1)

write.table(probeset.list1, file = "result/BLIA_BLIS_p.txt", sep="\t", quote = F, row.names=T)
write.table(probeset.list2, file = "result/BLIA_BLIS_p_lfc2.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list1)
head(genes)
sign_exp <- hdat[genes,]
head(sign_exp)
write.table(sign_exp, file = "result/BLIA_BLIS_exp.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list2)
head(genes)
sign_exp <- hdat[genes,]
head(sign_exp)
write.table(sign_exp, file = "result/BLIA_BLIS_exp_lfc2.txt", sep="\t", quote = F, row.names=T)

