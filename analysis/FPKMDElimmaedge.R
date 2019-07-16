# library(edgeR)
# library(stringi)
# source("http://bioconductor.org/biocLite.R")
# biocLite("Homo.sapiens")
# biocLite("biomaRt")
library(limma)
library(foreach)
library(doParallel)
library(pheatmap)
library(ggplot2)


setwd("/home/jiande/prjs/Analysis/BreastCancer/GeneExpression/FPKM.Methy_CASE")

#差异化 tnbc vs normal
rm(list=ls())

ctrlinput <- "GeneExpression.Normal83.txt"
treatinput <- "GeneExpression.TNBC83.txt"
ctrl.type <- "Normal"
treat.type <- "TNBC"

parameterout <- paste("Rscript.para",treat.type,ctrl.type,"txt", sep = ".")

adjust.pValue.threshold <- 0.05
logFC.threshold <- 2

ctrl <- read.csv(file = ctrlinput, sep = "\t")
treat <- read.csv(file = treatinput, sep = "\t")
ctrl.count <- dim(ctrl)[2]-1
treat.count <- dim(treat)[2]-1
treat.ctrl.exp <- merge(treat,ctrl,by = intersect(names(treat), names(ctrl)))
treat.ctrl.exp <- treat.ctrl.exp[!rowSums(treat.ctrl.exp==0) > ncol(treat.ctrl.exp)*.3,]

id_gene <- read.csv("/home/jiande/prjs/RefGene/hg38ref/GRCh38_gene-symbols.txt", sep = "\t")
id_gene[,3:6] <- NULL
treat.ctrl.exp[,1] <- gsub("\\..*$","",treat.ctrl.exp[,1])
treat.ctrl.exp <- merge(id_gene,treat.ctrl.exp, by.x="ensembl_gene_id", by.y="ID")

treat.ctrl.exp[,1] <- NULL
treat.ctrl.exp <- avereps(treat.ctrl.exp,ID=treat.ctrl.exp[,1])
rowname <- treat.ctrl.exp[,1]
treat.ctrl.exp <- treat.ctrl.exp[,2:dim(treat.ctrl.exp)[2]]
colname <- colnames(treat.ctrl.exp)

treat.ctrl.exp <- scale(matrix(as.numeric(treat.ctrl.exp),nrow=nrow(treat.ctrl.exp)), center=FALSE, scale=colSums(matrix(as.numeric(treat.ctrl.exp),nrow=nrow(treat.ctrl.exp)))/1000000)
colnames(treat.ctrl.exp) <- colname
rownames(treat.ctrl.exp) <- rowname
treat.ctrl.exp.log <- log2(treat.ctrl.exp+1)

save(treat.count,ctrl.count,ctrl.type,treat.type,adjust.pValue.threshold,logFC.threshold,
     file=parameterout,ascii=TRUE)
normalizeout <- paste("TPM",treat.type,ctrl.type,"txt", sep = ".")
normalizelogout <- paste("TPMlog2",treat.type,ctrl.type,"txt", sep = ".")
treat.ctrl.exp <- cbind(id=row.names(treat.ctrl.exp),treat.ctrl.exp) 
write.table(treat.ctrl.exp,file=normalizeout,sep="\t",row.names=F,quote=F)
treat.ctrl.exp.log <- cbind(id=row.names(treat.ctrl.exp.log),treat.ctrl.exp.log) 
write.table(treat.ctrl.exp.log,file=normalizelogout,sep="\t",row.names=F,quote=F)



load(parameterout)
normalizeout <- paste("TPM",treat.type,ctrl.type,"txt", sep = ".")
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)

treat.ctrl.grade <- c(rep(1,treat.count),rep(2,ctrl.count))
treat.ctrl.DE <- data.frame()

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

treat.ctrl.DE <- foreach(i = 1:nrow(treat.ctrl.exp), .combine=rbind) %dopar% {
  rt <- rbind(expression=treat.ctrl.exp[i,],grade=treat.ctrl.grade)
  rt <- as.matrix(t(rt))
  wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
  
  tumorGeneMeans <- mean(as.numeric(treat.ctrl.exp[i,1:treat.count]))
  normalGeneMeans <- mean(as.numeric(treat.ctrl.exp[i,(treat.count+1):ncol(treat.ctrl.exp)]))
  pvalue <- wilcoxTest$p.value
  tumorMed <- median(as.numeric(treat.ctrl.exp[i,1:treat.count]))
  normalMed <- median(as.numeric(treat.ctrl.exp[i,(treat.count+1):ncol(treat.ctrl.exp)]))
  diffMed <- tumorMed - normalMed
  diffMean <- tumorGeneMeans - normalGeneMeans
  logFC <- log2(tumorMed) - log2(normalMed)
  tempMatrix <- cbind(id=row.names(treat.ctrl.exp)[i],normalMean=normalGeneMeans,TumorMean=tumorGeneMeans,normalMed=normalMed,tumorMed=tumorMed,logFC=logFC,diffMean=diffMean,diffMed=diffMed,pValue=pvalue)
  
  tempMatrix #Equivalent to tnbc.normal.DE = rbind(tnbc.normal.DE, tempMatrix)
}
#stop cluster
stopCluster(cl)

pValue <- treat.ctrl.DE[,"pValue"]
fdr <- p.adjust(as.numeric(as.vector(pValue)),method="fdr")
treat.ctrl.DE <- cbind(treat.ctrl.DE,adjust.pValue=fdr)

write.table(treat.ctrl.DE,file=DEout,sep="\t",row.names=F,quote=F)



load(parameterout)
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
volcano.title <- paste(treat.type,"vs",ctrl.type, sep = " ")
volcano.out <- paste("volcano",treat.type,ctrl.type,"pdf", sep = ".")

treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
data <- treat.ctrl.DE[,c("id","logFC","adjust.pValue")]
data$color <- as.factor(ifelse(data$adjust.pValue<adjust.pValue.threshold & 
                                 abs(data$logFC)>=logFC.threshold,ifelse(data$logFC > logFC.threshold,'Up','Down'),'Not'))
p <- ggplot(data=data, aes(x=logFC, y=-log10(adjust.pValue), col=color)) +
  geom_point(alpha=0.8, size = 0.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adjust.pValue)",title = volcano.title) +
  geom_hline(yintercept = -log10(adjust.pValue.threshold), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.threshold, logFC.threshold), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("red", "blue", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(volcano.out,plot = p,dpi = 300)



load(parameterout)
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
normalizelogout <- paste("TPMlog2",treat.type,ctrl.type,"txt", sep = ".")
DEfilterout <- paste("limma.DE",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
outDiff <- treat.ctrl.DE[(abs(as.numeric(as.vector(treat.ctrl.DE[,"logFC"]))) > logFC.threshold & 
                            as.numeric(as.vector(treat.ctrl.DE[,"adjust.pValue"])) < adjust.pValue.threshold),]
treat.ctrl.exp.log <- read.table(normalizelogout,sep="\t",header=T,row.names=1,check.names=F)
heatmap <- treat.ctrl.exp.log[as.vector(outDiff[,1]),]
heatmap <- cbind(id=row.names(heatmap),heatmap)
write.table(outDiff,file=DEfilterout,sep="\t",row.names=F,quote=F)
write.table(heatmap,file=heatmaptxt,sep="\t",row.names=F,quote=F)



load(parameterout)
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
heatmaptiff <- paste("heatmap",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"tiff", sep = ".")
rt <- read.table(heatmaptxt,sep="\t",header=T,row.names=1,check.names=F)
Type <- c(rep(treat.type,treat.count),rep(ctrl.type,ctrl.count))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)
tiff(file=heatmaptiff,width=45,height=70,units="cm",compression="lzw",bg="white",res=300)
pheatmap(rt,annotation=Type,color=colorRampPalette(c("blue","white","red"))(50),
         cluster_cols=F,fontsize_row=5,fontsize_col=4,cutree_rows=2,show_rownames=F,show_colnames=F)
dev.off()

#差异化 nontnbc vs normal
rm(list=ls())

ctrlinput <- "GeneExpression.Normal83.txt"
treatinput <- "GeneExpression.nonTNBC597.605.txt"
ctrl.type <- "Normal"
treat.type <- "nonTNBC"

parameterout <- paste("Rscript.para",treat.type,ctrl.type,"txt", sep = ".")

adjust.pValue.threshold <- 0.05
logFC.threshold <- 2

ctrl <- read.csv(file = ctrlinput, sep = "\t")
treat <- read.csv(file = treatinput, sep = "\t")
ctrl.count <- dim(ctrl)[2]-1
treat.count <- dim(treat)[2]-1
treat.ctrl.exp <- merge(treat,ctrl,by = intersect(names(treat), names(ctrl)))
treat.ctrl.exp <- treat.ctrl.exp[!rowSums(treat.ctrl.exp==0) > ncol(treat.ctrl.exp)*.3,]

id_gene <- read.csv("/home/jiande/prjs/RefGene/hg38ref/GRCh38_gene-symbols.txt", sep = "\t")
id_gene[,3:6] <- NULL
treat.ctrl.exp[,1] <- gsub("\\..*$","",treat.ctrl.exp[,1])
treat.ctrl.exp <- merge(id_gene,treat.ctrl.exp, by.x="ensembl_gene_id", by.y="ID")

treat.ctrl.exp[,1] <- NULL
treat.ctrl.exp <- avereps(treat.ctrl.exp,ID=treat.ctrl.exp[,1])
rowname <- treat.ctrl.exp[,1]
treat.ctrl.exp <- treat.ctrl.exp[,2:dim(treat.ctrl.exp)[2]]
colname <- colnames(treat.ctrl.exp)

treat.ctrl.exp <- scale(matrix(as.numeric(treat.ctrl.exp),nrow=nrow(treat.ctrl.exp)), center=FALSE, scale=colSums(matrix(as.numeric(treat.ctrl.exp),nrow=nrow(treat.ctrl.exp)))/1000000)
colnames(treat.ctrl.exp) <- colname
rownames(treat.ctrl.exp) <- rowname
treat.ctrl.exp.log <- log2(treat.ctrl.exp+1)

save(treat.count,ctrl.count,ctrl.type,treat.type,adjust.pValue.threshold,logFC.threshold,
     file=parameterout,ascii=TRUE)
normalizeout <- paste("TPM",treat.type,ctrl.type,"txt", sep = ".")
normalizelogout <- paste("TPMlog2",treat.type,ctrl.type,"txt", sep = ".")
treat.ctrl.exp <- cbind(id=row.names(treat.ctrl.exp),treat.ctrl.exp) 
write.table(treat.ctrl.exp,file=normalizeout,sep="\t",row.names=F,quote=F)
treat.ctrl.exp.log <- cbind(id=row.names(treat.ctrl.exp.log),treat.ctrl.exp.log) 
write.table(treat.ctrl.exp.log,file=normalizelogout,sep="\t",row.names=F,quote=F)



load(parameterout)
normalizeout <- paste("TPM",treat.type,ctrl.type,"txt", sep = ".")
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)

treat.ctrl.grade <- c(rep(1,treat.count),rep(2,ctrl.count))
treat.ctrl.DE <- data.frame()

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

treat.ctrl.DE <- foreach(i = 1:nrow(treat.ctrl.exp), .combine=rbind) %dopar% {
  rt <- rbind(expression=treat.ctrl.exp[i,],grade=treat.ctrl.grade)
  rt <- as.matrix(t(rt))
  wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
  
  tumorGeneMeans <- mean(as.numeric(treat.ctrl.exp[i,1:treat.count]))
  normalGeneMeans <- mean(as.numeric(treat.ctrl.exp[i,(treat.count+1):ncol(treat.ctrl.exp)]))
  pvalue <- wilcoxTest$p.value
  tumorMed <- median(as.numeric(treat.ctrl.exp[i,1:treat.count]))
  normalMed <- median(as.numeric(treat.ctrl.exp[i,(treat.count+1):ncol(treat.ctrl.exp)]))
  diffMed <- tumorMed - normalMed
  diffMean <- tumorGeneMeans - normalGeneMeans
  logFC <- log2(tumorMed) - log2(normalMed)
  tempMatrix <- cbind(id=row.names(treat.ctrl.exp)[i],normalMean=normalGeneMeans,TumorMean=tumorGeneMeans,normalMed=normalMed,tumorMed=tumorMed,logFC=logFC,diffMean=diffMean,diffMed=diffMed,pValue=pvalue)
  
  tempMatrix #Equivalent to tnbc.normal.DE = rbind(tnbc.normal.DE, tempMatrix)
}
#stop cluster
stopCluster(cl)

pValue <- treat.ctrl.DE[,"pValue"]
fdr <- p.adjust(as.numeric(as.vector(pValue)),method="fdr")
treat.ctrl.DE <- cbind(treat.ctrl.DE,adjust.pValue=fdr)

write.table(treat.ctrl.DE,file=DEout,sep="\t",row.names=F,quote=F)



load(parameterout)
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
volcano.title <- paste(treat.type,"vs",ctrl.type, sep = " ")
volcano.out <- paste("volcano",treat.type,ctrl.type,"pdf", sep = ".")

treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
data <- treat.ctrl.DE[,c("id","logFC","adjust.pValue")]
data$color <- as.factor(ifelse(data$adjust.pValue<adjust.pValue.threshold & 
                                 abs(data$logFC)>=logFC.threshold,ifelse(data$logFC > logFC.threshold,'Up','Down'),'Not'))
p <- ggplot(data=data, aes(x=logFC, y=-log10(adjust.pValue), col=color)) +
  geom_point(alpha=0.8, size = 0.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adjust.pValue)",title = volcano.title) +
  geom_hline(yintercept = -log10(adjust.pValue.threshold), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.threshold, logFC.threshold), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("red", "blue", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(volcano.out,plot = p,dpi = 300)



load(parameterout)
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
normalizelogout <- paste("TPMlog2",treat.type,ctrl.type,"txt", sep = ".")
DEfilterout <- paste("limma.DE",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
outDiff <- treat.ctrl.DE[(abs(as.numeric(as.vector(treat.ctrl.DE[,"logFC"]))) > logFC.threshold & 
                            as.numeric(as.vector(treat.ctrl.DE[,"adjust.pValue"])) < adjust.pValue.threshold),]
treat.ctrl.exp.log <- read.table(normalizelogout,sep="\t",header=T,row.names=1,check.names=F)
heatmap <- treat.ctrl.exp.log[as.vector(outDiff[,1]),]
heatmap <- cbind(id=row.names(heatmap),heatmap)
write.table(outDiff,file=DEfilterout,sep="\t",row.names=F,quote=F)
write.table(heatmap,file=heatmaptxt,sep="\t",row.names=F,quote=F)



load(parameterout)
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
heatmaptiff <- paste("heatmap",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"tiff", sep = ".")
rt <- read.table(heatmaptxt,sep="\t",header=T,row.names=1,check.names=F)
Type <- c(rep(treat.type,treat.count),rep(ctrl.type,ctrl.count))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)
tiff(file=heatmaptiff,width=45,height=70,units="cm",compression="lzw",bg="white",res=300)
pheatmap(rt,annotation=Type,color=colorRampPalette(c("blue","white","red"))(50),
         cluster_cols=F,fontsize_row=5,fontsize_col=4,cutree_rows=2,show_rownames=F,show_colnames=F)
dev.off()


#差异化 tnbc vs nontnbc
rm(list=ls())

ctrlinput <- "GeneExpression.nonTNBC597.605.txt"
treatinput <- "GeneExpression.TNBC83.txt"
ctrl.type <- "nonTNBC"
treat.type <- "TNBC"

parameterout <- paste("Rscript.para",treat.type,ctrl.type,"txt", sep = ".")

adjust.pValue.threshold <- 0.05
logFC.threshold <- 2

ctrl <- read.csv(file = ctrlinput, sep = "\t")
treat <- read.csv(file = treatinput, sep = "\t")
ctrl.count <- dim(ctrl)[2]-1
treat.count <- dim(treat)[2]-1
treat.ctrl.exp <- merge(treat,ctrl,by = intersect(names(treat), names(ctrl)))
treat.ctrl.exp <- treat.ctrl.exp[!rowSums(treat.ctrl.exp==0) > ncol(treat.ctrl.exp)*.3,]

id_gene <- read.csv("/home/jiande/prjs/RefGene/hg38ref/GRCh38_gene-symbols.txt", sep = "\t")
id_gene[,3:6] <- NULL
treat.ctrl.exp[,1] <- gsub("\\..*$","",treat.ctrl.exp[,1])
treat.ctrl.exp <- merge(id_gene,treat.ctrl.exp, by.x="ensembl_gene_id", by.y="ID")

treat.ctrl.exp[,1] <- NULL
treat.ctrl.exp <- avereps(treat.ctrl.exp,ID=treat.ctrl.exp[,1])
rowname <- treat.ctrl.exp[,1]
treat.ctrl.exp <- treat.ctrl.exp[,2:dim(treat.ctrl.exp)[2]]
colname <- colnames(treat.ctrl.exp)

treat.ctrl.exp <- scale(matrix(as.numeric(treat.ctrl.exp),nrow=nrow(treat.ctrl.exp)), center=FALSE, scale=colSums(matrix(as.numeric(treat.ctrl.exp),nrow=nrow(treat.ctrl.exp)))/1000000)
colnames(treat.ctrl.exp) <- colname
rownames(treat.ctrl.exp) <- rowname
treat.ctrl.exp.log <- log2(treat.ctrl.exp+1)

save(treat.count,ctrl.count,ctrl.type,treat.type,adjust.pValue.threshold,logFC.threshold,
     file=parameterout,ascii=TRUE)
normalizeout <- paste("TPM",treat.type,ctrl.type,"txt", sep = ".")
normalizelogout <- paste("TPMlog2",treat.type,ctrl.type,"txt", sep = ".")
treat.ctrl.exp <- cbind(id=row.names(treat.ctrl.exp),treat.ctrl.exp) 
write.table(treat.ctrl.exp,file=normalizeout,sep="\t",row.names=F,quote=F)
treat.ctrl.exp.log <- cbind(id=row.names(treat.ctrl.exp.log),treat.ctrl.exp.log) 
write.table(treat.ctrl.exp.log,file=normalizelogout,sep="\t",row.names=F,quote=F)




load(parameterout)
normalizeout <- paste("TPM",treat.type,ctrl.type,"txt", sep = ".")
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)

treat.ctrl.grade <- c(rep(1,treat.count),rep(2,ctrl.count))
treat.ctrl.DE <- data.frame()

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

treat.ctrl.DE <- foreach(i = 1:nrow(treat.ctrl.exp), .combine=rbind) %dopar% {
  rt <- rbind(expression=treat.ctrl.exp[i,],grade=treat.ctrl.grade)
  rt <- as.matrix(t(rt))
  wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
  
  tumorGeneMeans <- mean(as.numeric(treat.ctrl.exp[i,1:treat.count]))
  normalGeneMeans <- mean(as.numeric(treat.ctrl.exp[i,(treat.count+1):ncol(treat.ctrl.exp)]))
  pvalue <- wilcoxTest$p.value
  tumorMed <- median(as.numeric(treat.ctrl.exp[i,1:treat.count]))
  normalMed <- median(as.numeric(treat.ctrl.exp[i,(treat.count+1):ncol(treat.ctrl.exp)]))
  diffMed <- tumorMed - normalMed
  diffMean <- tumorGeneMeans - normalGeneMeans
  logFC <- log2(tumorMed) - log2(normalMed)
  tempMatrix <- cbind(id=row.names(treat.ctrl.exp)[i],normalMean=normalGeneMeans,TumorMean=tumorGeneMeans,normalMed=normalMed,tumorMed=tumorMed,logFC=logFC,diffMean=diffMean,diffMed=diffMed,pValue=pvalue)
  
  tempMatrix #Equivalent to tnbc.normal.DE = rbind(tnbc.normal.DE, tempMatrix)
}
#stop cluster
stopCluster(cl)

pValue <- treat.ctrl.DE[,"pValue"]
fdr <- p.adjust(as.numeric(as.vector(pValue)),method="fdr")
treat.ctrl.DE <- cbind(treat.ctrl.DE,adjust.pValue=fdr)

write.table(treat.ctrl.DE,file=DEout,sep="\t",row.names=F,quote=F)



load(parameterout)
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
volcano.title <- paste(treat.type,"vs",ctrl.type, sep = " ")
volcano.out <- paste("volcano",treat.type,ctrl.type,"pdf", sep = ".")

treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
data <- treat.ctrl.DE[,c("id","logFC","adjust.pValue")]
data$color <- as.factor(ifelse(data$adjust.pValue<adjust.pValue.threshold & 
                                 abs(data$logFC)>=logFC.threshold,ifelse(data$logFC > logFC.threshold,'Up','Down'),'Not'))
p <- ggplot(data=data, aes(x=logFC, y=-log10(adjust.pValue), col=color)) +
  geom_point(alpha=0.8, size = 0.8) +
  theme_bw(base_size = 15) +
  labs(x="log2 (fold change)",y="-log10 (adjust.pValue)",title = volcano.title) +
  geom_hline(yintercept = -log10(adjust.pValue.threshold), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-logFC.threshold, logFC.threshold), lty=4,col="black",lwd=0.6) +
  scale_color_manual(name = "", values = c("red", "blue", "gray"), limits = c("Up", "Down", "Not")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(volcano.out,plot = p,dpi = 300)



load(parameterout)
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
normalizelogout <- paste("TPMlog2",treat.type,ctrl.type,"txt", sep = ".")
DEfilterout <- paste("limma.DE",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
outDiff <- treat.ctrl.DE[(abs(as.numeric(as.vector(treat.ctrl.DE[,"logFC"]))) > logFC.threshold & 
                            as.numeric(as.vector(treat.ctrl.DE[,"adjust.pValue"])) < adjust.pValue.threshold),]
treat.ctrl.exp.log <- read.table(normalizelogout,sep="\t",header=T,row.names=1,check.names=F)
heatmap <- treat.ctrl.exp.log[as.vector(outDiff[,1]),]
heatmap <- cbind(id=row.names(heatmap),heatmap)
write.table(outDiff,file=DEfilterout,sep="\t",row.names=F,quote=F)
write.table(heatmap,file=heatmaptxt,sep="\t",row.names=F,quote=F)



load(parameterout)
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"txt", sep = ".")
heatmaptiff <- paste("heatmap",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold, sep = "_"),"tiff", sep = ".")
rt <- read.table(heatmaptxt,sep="\t",header=T,row.names=1,check.names=F)
Type <- c(rep(treat.type,treat.count),rep(ctrl.type,ctrl.count))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)
tiff(file=heatmaptiff,width=45,height=70,units="cm",compression="lzw",bg="white",res=300)
pheatmap(rt,annotation=Type,color=colorRampPalette(c("blue","white","red"))(50),
         cluster_cols=F,fontsize_row=5,fontsize_col=4,cutree_rows=2,show_rownames=F,show_colnames=F)
dev.off()


