library(limma)

setwd("/home/jiande/prjs/BreastTumor/analyse_data/")

list1 <- read.csv(file = "tnbc_ctrl/Treat_Ctrl_p_all.txt", row.names = 1, sep = "\t")
list2 <- read.csv(file = "tnbc.all_ctrl/Treat_Ctrl_p_all.txt", row.names = 1, sep = "\t")
gene <- unique(c(rownames(list1),rownames(list2)))
a <- cbind(rep(FALSE,length(gene)),rep(FALSE,length(gene)))
rownames(a) <- gene

gene1 <- list1[(list1$logFC>=2|list1$logFC<=-2)&list1$adj.P.Val<=0.001,]
gene2 <- list2[(list2$logFC>=2|list2$logFC<=-2)&list2$adj.P.Val<=0.001,]

a[rownames(gene1),1] <- TRUE
a[rownames(gene2),2] <- TRUE
aa <- vennCounts(a)
vennDiagram(aa, include = "both", names = c("Only Treat_Ctrl","All Treat_Ctrl"), counts.col = "red")

