mapData <- read.table(file="/media/yuhua/yuhua_projects/genomeanno/gencode.vM17.genetype.tab",header=FALSE,stringsAsFactors=FALSE)
mapData[,1] <- gsub("\\.\\d+","",mapData[,1])

netDataGSR <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_GSR_final.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(netDataGSR[,2],mapData[,1])
netDataGSR <- cbind(netDataGSR[,1],mapData[matchIndexes,3],netDataGSR[,2],netDataGSR[,c(3,4,5,6)])
write.csv(netDataGSR,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_target_gene_regulatory_network_Wang_et_al_study.csv",row.names=FALSE,col.names=TRUE,quote=FALSE)

netDataXW <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_XW_final.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(netDataXW[,2],mapData[,1])
netDataXW <- cbind(netDataXW[,1],mapData[matchIndexes,3],netDataXW[,2],netDataXW[,c(3,4,5,6)])
write.csv(netDataXW,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_target_gene_regulatory_network_Wu_et_al_study.csv",row.names=FALSE,col.names=TRUE,quote=FALSE)