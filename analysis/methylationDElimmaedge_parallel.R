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


setwd("/home/jiande/prjs/Analysis/BreastCancer/Methylation/Methy.FPKM_CASE")

#差异甲基化 tnbc vs normal
ctrlinput <- "Methylation.Normal83.txt"
treatinput <- "Methylation.TNBC83.84.txt"
treat.type <- "TNBC"
ctrl.type <- "Normal"

normalizeout <- paste("normalizeMethy",treat.type,ctrl.type,"txt", sep = ".")
parameterout <- paste("Rscript.para",treat.type,ctrl.type,"txt", sep = ".")
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
volcano.title <- paste(treat.type,"vs",ctrl.type, sep = " ")
volcano.out <- paste("volcano",treat.type,ctrl.type,"pdf", sep = ".")

ctrl <- read.csv(file = ctrlinput, sep = "\t")
treat <- read.csv(file = treatinput, sep = "\t")
ctrl.count <- dim(ctrl)[2]-1
treat.count <- dim(treat)[2]-1
treat.ctrl.exp <- merge(treat,ctrl,by = intersect(names(treat), names(ctrl)))
rownames(treat.ctrl.exp) <- treat.ctrl.exp[,1]
treat.ctrl.exp[,1] <- NULL

treat.ctrl.exp <- avereps(treat.ctrl.exp)
treat.ctrl.exp <- treat.ctrl.exp[rowMeans(treat.ctrl.exp)>0,]
treat.ctrl.exp <- na.omit(treat.ctrl.exp)
treat.ctrl.exp <- normalizeBetweenArrays(treat.ctrl.exp)

treat.ctrl.exp <- cbind(id=row.names(treat.ctrl.exp),treat.ctrl.exp) #为了让id号写入文件
write.table(treat.ctrl.exp,file=normalizeout,sep="\t",row.names=F,quote=F)
save(treat.count,ctrl.count,file=parameterout,ascii=TRUE)

treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)
load(parameterout)
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

adjust.pValue.threshold <- 0.05
logFC.threshold <- 1
diffMean.threshold = 0.2

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

DEfilterout <- paste("limma.DE",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"txt", sep = ".")
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"txt", sep = ".")
heatmaptiff <- paste("heatmap",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"tiff", sep = ".")

treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
outDiff <- treat.ctrl.DE[(abs(as.numeric(as.vector(treat.ctrl.DE[,"diffMean"]))) > diffMean.threshold & 
                            abs(as.numeric(as.vector(treat.ctrl.DE[,"logFC"]))) > logFC.threshold & 
                            as.numeric(as.vector(treat.ctrl.DE[,"adjust.pValue"])) < adjust.pValue.threshold),]
treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)
heatmap <- treat.ctrl.exp[as.vector(outDiff[,1]),]
heatmap <- cbind(id=row.names(heatmap),heatmap)
write.table(outDiff,file=DEfilterout,sep="\t",row.names=F,quote=F)
write.table(heatmap,file=heatmaptxt,sep="\t",row.names=F,quote=F)

rt <- read.table(heatmaptxt,sep="\t",header=T,row.names=1,check.names=F)
Type <- c(rep(treat.type,treat.count),rep(ctrl.type,ctrl.count))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)
tiff(file=heatmaptiff,width=45,height=70,units="cm",compression="lzw",bg="white",res=300)
pheatmap(rt,annotation=Type,color=colorRampPalette(c("blue","white","red"))(50),
         cluster_cols=F,fontsize_row=5,fontsize_col=4,cutree_rows=2,show_rownames=F,show_colnames=F)
dev.off()


#差异甲基化 tnbc vs nontnbc
rm(list=ls())
ctrlinput <- "Methylation.nonTNBC597.604.txt"
treatinput <- "Methylation.TNBC83.84.txt"
treat.type <- "TNBC"
ctrl.type <- "nonTNBC"

normalizeout <- paste("normalizeMethy",treat.type,ctrl.type,"txt", sep = ".")
parameterout <- paste("Rscript.para",treat.type,ctrl.type,"txt", sep = ".")
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
volcano.title <- paste(treat.type,"vs",ctrl.type, sep = " ")
volcano.out <- paste("volcano",treat.type,ctrl.type,"pdf", sep = ".")

ctrl <- read.csv(file = ctrlinput, sep = "\t")
treat <- read.csv(file = treatinput, sep = "\t")
ctrl.count <- dim(ctrl)[2]-1
treat.count <- dim(treat)[2]-1
treat.ctrl.exp <- merge(treat,ctrl,by = intersect(names(treat), names(ctrl)))
rownames(treat.ctrl.exp) <- treat.ctrl.exp[,1]
treat.ctrl.exp[,1] <- NULL

treat.ctrl.exp <- avereps(treat.ctrl.exp)
treat.ctrl.exp <- treat.ctrl.exp[rowMeans(treat.ctrl.exp)>0,]
treat.ctrl.exp <- na.omit(treat.ctrl.exp)
treat.ctrl.exp <- normalizeBetweenArrays(treat.ctrl.exp)

treat.ctrl.exp <- cbind(id=row.names(treat.ctrl.exp),treat.ctrl.exp) #为了让id号写入文件
write.table(treat.ctrl.exp,file=normalizeout,sep="\t",row.names=F,quote=F)
save(treat.count,ctrl.count,file=parameterout,ascii=TRUE)

treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)
load(parameterout)
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

adjust.pValue.threshold <- 0.05
logFC.threshold <- 1
diffMean.threshold = 0.2

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

DEfilterout <- paste("limma.DE",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"txt", sep = ".")
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"txt", sep = ".")
heatmaptiff <- paste("heatmap",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"tiff", sep = ".")

treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
outDiff <- treat.ctrl.DE[(abs(as.numeric(as.vector(treat.ctrl.DE[,"diffMean"]))) > diffMean.threshold & 
                            abs(as.numeric(as.vector(treat.ctrl.DE[,"logFC"]))) > logFC.threshold & 
                            as.numeric(as.vector(treat.ctrl.DE[,"adjust.pValue"])) < adjust.pValue.threshold),]
treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)
heatmap <- treat.ctrl.exp[as.vector(outDiff[,1]),]
heatmap <- cbind(id=row.names(heatmap),heatmap)
write.table(outDiff,file=DEfilterout,sep="\t",row.names=F,quote=F)
write.table(heatmap,file=heatmaptxt,sep="\t",row.names=F,quote=F)

rt <- read.table(heatmaptxt,sep="\t",header=T,row.names=1,check.names=F)
Type <- c(rep(treat.type,treat.count),rep(ctrl.type,ctrl.count))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)
tiff(file=heatmaptiff,width=45,height=70,units="cm",compression="lzw",bg="white",res=300)
pheatmap(rt,annotation=Type,color=colorRampPalette(c("blue","white","red"))(50),
         cluster_cols=F,fontsize_row=5,fontsize_col=4,cutree_rows=2,show_rownames=F,show_colnames=F)
dev.off()

#差异甲基化 nontnbc vs normal
rm(list=ls())
ctrlinput <- "Methylation.Normal83.txt"
treatinput <- "Methylation.nonTNBC597.604.txt"
treat.type <- "nonTNBC"
ctrl.type <- "Normal"

normalizeout <- paste("normalizeMethy",treat.type,ctrl.type,"txt", sep = ".")
parameterout <- paste("Rscript.para",treat.type,ctrl.type,"txt", sep = ".")
DEout <- paste("limma.DE",treat.type,ctrl.type,"txt", sep = ".")
volcano.title <- paste(treat.type,"vs",ctrl.type, sep = " ")
volcano.out <- paste("volcano",treat.type,ctrl.type,"pdf", sep = ".")

ctrl <- read.csv(file = ctrlinput, sep = "\t")
treat <- read.csv(file = treatinput, sep = "\t")
ctrl.count <- dim(ctrl)[2]-1
treat.count <- dim(treat)[2]-1
treat.ctrl.exp <- merge(treat,ctrl,by = intersect(names(treat), names(ctrl)))
rownames(treat.ctrl.exp) <- treat.ctrl.exp[,1]
treat.ctrl.exp[,1] <- NULL

treat.ctrl.exp <- avereps(treat.ctrl.exp)
treat.ctrl.exp <- treat.ctrl.exp[rowMeans(treat.ctrl.exp)>0,]
treat.ctrl.exp <- na.omit(treat.ctrl.exp)
treat.ctrl.exp <- normalizeBetweenArrays(treat.ctrl.exp)

treat.ctrl.exp <- cbind(id=row.names(treat.ctrl.exp),treat.ctrl.exp) #为了让id号写入文件
write.table(treat.ctrl.exp,file=normalizeout,sep="\t",row.names=F,quote=F)
save(treat.count,ctrl.count,file=parameterout,ascii=TRUE)

treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)
load(parameterout)
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

adjust.pValue.threshold <- 0.05
logFC.threshold <- 1
diffMean.threshold = 0.2

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

DEfilterout <- paste("limma.DE",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"txt", sep = ".")
heatmaptxt <- paste("heatmap",treat.type,ctrl.type,
                    paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"txt", sep = ".")
heatmaptiff <- paste("heatmap",treat.type,ctrl.type,
                     paste("filt",logFC.threshold,adjust.pValue.threshold,diffMean.threshold, sep = "_"),"tiff", sep = ".")

treat.ctrl.DE <- read.table(DEout,sep="\t",header=T,check.names=F)
outDiff <- treat.ctrl.DE[(abs(as.numeric(as.vector(treat.ctrl.DE[,"diffMean"]))) > diffMean.threshold & 
                            abs(as.numeric(as.vector(treat.ctrl.DE[,"logFC"]))) > logFC.threshold & 
                            as.numeric(as.vector(treat.ctrl.DE[,"adjust.pValue"])) < adjust.pValue.threshold),]
treat.ctrl.exp <- read.table(normalizeout,sep="\t",header=T,row.names=1,check.names=F)
heatmap <- treat.ctrl.exp[as.vector(outDiff[,1]),]
heatmap <- cbind(id=row.names(heatmap),heatmap)
write.table(outDiff,file=DEfilterout,sep="\t",row.names=F,quote=F)
write.table(heatmap,file=heatmaptxt,sep="\t",row.names=F,quote=F)

rt <- read.table(heatmaptxt,sep="\t",header=T,row.names=1,check.names=F)
Type <- c(rep(treat.type,treat.count),rep(ctrl.type,ctrl.count))
names(Type) <- colnames(rt)
Type <- as.data.frame(Type)
tiff(file=heatmaptiff,width=45,height=70,units="cm",compression="lzw",bg="white",res=300)
pheatmap(rt,annotation=Type,color=colorRampPalette(c("blue","white","red"))(50),
         cluster_cols=F,fontsize_row=5,fontsize_col=4,cutree_rows=2,show_rownames=F,show_colnames=F)
dev.off()





###################################################################
#甲基化与转录组数据融合
methy <- read.table(file="normalizeMethy.txt",row.names=1,header=T,sep="\t",check.names=F)
rna <- read.table(file="normalizeExp.txt",row.name=1,header=T,sep="\t",check.names=F)

colnames(methy) <- gsub("(.*?)-(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3-\\4",colnames(methy))
colnames(rna) <- gsub("(.*?)-(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3-\\4",colnames(rna))
sameSample <- intersect(colnames(methy),colnames(rna))
rownames(methy) <- paste(rownames(methy),"methy",sep="|")
rownames(rna) <- paste(rownames(rna),"exp",sep="|")
merge <- rbind(id=sameSample,methy[,sameSample],rna[,sameSample])
write.table(merge,file="merge.txt",sep="\t",quote=F,col.names=F)

#基因甲基化与转录组数据相关关系
inputFile <- "merge.txt"
methyGene <- "NCAPD3|methy"
expGene <- "NCAPD3|exp"
rt <- read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
i <- rt[methyGene,]
j <- rt[expGene,]
x <- as.numeric(i)
y <- log2(as.numeric(j)+1)
corT <- cor.test(x,y)
methyGeneName <- unlist(strsplit(methyGene,"\\|"))[[1]]
expGeneName <- unlist(strsplit(expGene,"\\|"))[[1]]
z <- lm(y~x)
estimate <- corT$estimate
cor <- round(estimate,3)
pvalue <- corT$p.value   #以小数的方式保留位数
pval <- signif(pvalue,4) #以科学计数法的方式保留位数
pval <- format(pval,scientific=T)

outTiff <- "cor.tiff"
tiff(file=outTiff,width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
plot(x,y,type="p",pch=16,main=paste("Cor=",cor,"(p-value=",pval,")",sep=""),xlab=paste(methyGeneName,"methylation"),ylab=paste(expGeneName,"expression"))
lines(x,fitted(z),col=2)
dev.off()

#甲基化位点提取
# nohup ./posget.pl metadata.cart.2018-08-01.json NCAPD3 > /home/train/Liuping/raw_data/TCGA_Methylation_ProstateCancer/NCAPD3_Methy_location.log 2>&1 &
  #[1] 40635
  
#甲基化位点与gene表达合并
methyFile="posMethy.txt"
expFile="normalizeExp.txt"
methy <- read.table(methyFile,row.names=1,header=T,sep="\t",check.names=F)
rna <- read.table(expFile,row.names=1,header=T,sep="\t",check.names=F)
colnames(methy) <- gsub("(.*?)-(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3-\\4",colnames(methy))
colnames(rna) <- gsub("(.*?)-(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3-\\4",colnames(rna))
sameSample <- intersect(colnames(methy),colnames(rna))
rownames(methy) <- paste(rownames(methy),"methy",sep="|")
rownames(rna) <- paste(rownames(rna),"exp",sep="|")
merge <- rbind(id=sameSample,methy[,sameSample],rna[,sameSample])

#甲基化位点与gene表达相关
inputFile <- "merge.txt"
posmethy <- "cg14934141|methy"
expGene <- "NCAPD3|exp"
rt <- read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
i <- rt[posmethy,]
j <- rt[expGene,]
x <- as.numeric(i)
y <- log2(as.numeric(j)+1)
corT <- cor.test(x,y)
posmethyid <- unlist(strsplit(posmethy,"\\|"))[[1]]
expgenename <- unlist(strsplit(expGene,"\\|"))[[1]]
z <- lm(y~x)
cor <- corT$estimate
cor <- round(cor,3)
pvalue <- corT$p.value
pval <- signif(pvalue,4)
pval <- format(pval,scientific=T)
outTiff <- "posCor.tiff"
tiff(file=outTiff,width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
plot(x,y,type="p",pch=16,main=paste("Cor=",cor,"(p-value=",pval,")",sep=""),xlab=paste(posmethyid,"methylation"),ylab=paste(expgenename,"expression"))
lines(x,fitted(z),col=2)
dev.off()
