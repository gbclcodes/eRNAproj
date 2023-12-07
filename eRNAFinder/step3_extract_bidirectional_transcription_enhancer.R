enh.expr.gsr <- read.csv(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_EPM.csv",header=TRUE,stringsAsFactors=FALSE)
td.enh.info <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/filtered_mm10_td.bed",sep=" ",header=FALSE,stringsAsFactors=FALSE)

total.enh.pos <- paste(enh.expr.gsr[,2],enh.expr.gsr[,3],enh.expr.gsr[,4],sep="_")
td.enh.pos <- paste(td.enh.info[,1],td.enh.info[,2],td.enh.info[,3],sep="_")

match.indexes <- match(td.enh.pos,total.enh.pos)
match.indexes <- match.indexes[which(!is.na(match.indexes))]
enh.expr.gsr <- enh.expr.gsr[match.indexes,]
enh.expr.gsr <- enh.expr.gsr[which(rowSums(enh.expr.gsr[6:ncol(enh.expr.gsr)])/9 > 0.01),]
write.table(enh.expr.gsr,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_clean_td_0.01EPM.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)


enh.expr.xw <- read.csv(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_EPM.csv",header=TRUE,stringsAsFactors=FALSE)
td.enh.info <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/filtered_mm10_td.bed",sep=" ",header=FALSE,stringsAsFactors=FALSE)

total.enh.pos <- paste(enh.expr.xw[,2],enh.expr.xw[,3],enh.expr.xw[,4],sep="_")
td.enh.pos <- paste(td.enh.info[,1],td.enh.info[,2],td.enh.info[,3],sep="_")

match.indexes <- match(td.enh.pos,total.enh.pos)
match.indexes <- match.indexes[which(!is.na(match.indexes))]
enh.expr.xw <- enh.expr.xw[match.indexes,]
enh.expr.xw <- enh.expr.xw[which(rowSums(enh.expr.xw[6:ncol(enh.expr.xw)])/7 > 0.01),]
write.table(enh.expr.xw,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_clean_td_0.01EPM.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)
