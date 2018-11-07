
library(limma)
data1=removeBatchEffect(data,batch=c("0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"
))
row.names(data.matrix(data1))=rname
View(data1
     )
write.table(data1, file="C:/Users/tmamid/Desktop/all_data.txt", quote=F, sep="\t")
