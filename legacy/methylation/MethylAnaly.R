library(MethylMix)
library(doParallel)
## Loading required package: foreach
## Loading required package: iterators
## Loading required package: parallel
# data(METcancer)
# data(METnormal)
# data(GEcancer)
# head(METcancer[, 1:4])

setwd("~/prjs/BreastTumor/methylation")

METcancer <- read.csv(file = "./TNBC/combin.txt", row.names = 1, sep = "\t")
METnormal <- read.csv(file = "./Normal/combin.txt", row.names = 1, sep = "\t")
GEcancer <- read.csv(file = "./TNBC/GEtrip_syb.txt", row.names = 1, sep = "\t")

METcancer1 <- data.matrix(METcancer)
METnormal1 <- data.matrix(METnormal)
GEcancer1 <- data.matrix(GEcancer)

cl <- makeCluster(5)
registerDoParallel(cl)
MethylMixResults <- MethylMix(METcancer1, GEcancer1, METnormal1)
stopCluster(cl)

plots <- MethylMix_PlotModel("A2M", MethylMixResults, METcancer)
MethylMixResults$MethylationStates[, 1:5]
MethylMixResults$Classifications[, 1:5]
