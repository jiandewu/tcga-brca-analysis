# source("http://bioconductor.org/biocLite.R")
# biocLite("rgl")
# biocLite("limma")
# biocLite("edgeR")
library(limma)
library("rgl")
library(gplots)
# library(edgeR)


setwd("~/prjs/BreastTumor/HicksRaw")
hdat <- read.csv(file = "BLIS_normal.txt", row.names = 1, sep = "\t")
pheno <- read.csv(file = "BLIS_normal_pheno.txt", sep = "\t")

# hdat=removeBatchEffect(hdat,batch=c("0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"
# ))
colm <- colMeans(hdat)
diff <- mean(colm[1:60])-mean(colm[61:203])
hdat[,c(1:60)] <- hdat[,c(1:60)]-diff
boxplot(hdat)


rownames(pheno) <- pheno$sample
colnames(hdat) <- rownames(pheno)
rownames(pheno) == colnames(hdat)

##### PCA
sample.groups <- c(rep("red", 60), rep("blue",143))
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
volcanoplot(fit1,  main = "BLIA vs normal")

nrow(topTable(fit1, coef=1, number=54675, lfc=2, p.value = 0.05))
probeset.list1 <- topTable(fit1, coef=1, number=54675, lfc=2, p.value = 0.05)
head(probeset.list1)
write.table(probeset.list1, file = "result/BLIS_normal_p.txt", sep="\t", quote = F, row.names=T)

genes <- rownames(probeset.list1)
head(genes)
sign_exp <- hdat[genes,]
head(sign_exp)
write.table(sign_exp, file = "result/BLIS_normal_exp.txt", sep="\t", quote = F, row.names=T)

