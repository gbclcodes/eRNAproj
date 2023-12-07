mapData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE_ID_to_NAME.txt",header=FALSE,stringsAsFactors=FALSE,sep="\t")
tf.enh.motif.pairs <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile/TF_enhancer_network.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(tf.enh.motif.pairs[,1],mapData[,1])
tf.enh.motif.pairs <- cbind(mapData[matchIndexes[which(!is.na(matchIndexes))],2],tf.enh.motif.pairs[which(!is.na(matchIndexes)),2:ncol(tf.enh.motif.pairs)])
tf.enh.motif.pairs.yes <- tf.enh.motif.pairs[which(tf.enh.motif.pairs[,5]=="Yes"),]
tf.enh.motif.pairs.no <- tf.enh.motif.pairs[which(tf.enh.motif.pairs[,5]=="No"),]  

tf.enh.co.pairs <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_GSR.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(paste(toupper(tf.enh.co.pairs[,5]),tf.enh.co.pairs[,1],sep="\t"),paste(toupper(tf.enh.motif.pairs.no[,1]),tf.enh.motif.pairs.no[,2],sep="\t"))
tf.enh.motif.pairs.no <- tf.enh.motif.pairs.no[matchIndexes[which(!is.na(matchIndexes))],]
tf.enh.motif.pairs <- rbind(tf.enh.motif.pairs.yes,tf.enh.motif.pairs.no)

tf.enh.co.pairs.all <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_GSR_all.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(paste(toupper(tf.enh.co.pairs.all[,5]),tf.enh.co.pairs.all[,1],sep="\t"),paste(toupper(tf.enh.motif.pairs[,1]),tf.enh.motif.pairs[,2],sep="\t"))
tf.enh.co.pairs.all <- tf.enh.co.pairs.all[which(!is.na(matchIndexes)),]
tf.enh.motif.pairs <- tf.enh.motif.pairs[matchIndexes[which(!is.na(matchIndexes))],]
tf.enh.co.pairs.all[which(tf.enh.co.pairs.all[,4]==0),4] <- min(tf.enh.co.pairs.all[which(tf.enh.co.pairs.all[,4]>0),4])
tf.enh.motif.pairs <- cbind(tf.enh.motif.pairs[,c(2,1,3,4)],tf.enh.motif.pairs[,1])
tf.enh.motif.pairs[,2] <- tf.enh.co.pairs.all[,2]

row.indexes <- c()
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
geneexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/stringtiefile_GSR/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(geneexpData) <- gsub("\\.\\d+","",geneexpData[,1])
geneexpData <- as.matrix(geneexpData[,3:ncol(geneexpData)])
for (rowindex in seq(1,nrow(tf.enh.motif.pairs))){
	rownum <- match(tf.enh.motif.pairs[rowindex,2],rownames(geneexpData))
	colnum <- match(tf.enh.motif.pairs[rowindex,4],colnames(geneexpData))
	if(geneexpData[tf.enh.motif.pairs[rowindex,2],tf.enh.motif.pairs[rowindex,4]] > 0){
		tf.enh.motif.pairs[rowindex,3] <- tf.enh.motif.pairs[rowindex,3]*geneexpData[tf.enh.motif.pairs[rowindex,2],tf.enh.motif.pairs[rowindex,4]]
		row.indexes <- c(row.indexes,rowindex)
	}
}
tf.enh.motif.pairs <- tf.enh.motif.pairs[row.indexes,]
write.table(tf.enh.motif.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile/TF_enhancer_network_final.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

mapData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE_ID_to_NAME.txt",header=FALSE,stringsAsFactors=FALSE,sep="\t")
tf.enh.motif.pairs <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile/TF_enhancer_network.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(tf.enh.motif.pairs[,1],mapData[,1])
tf.enh.motif.pairs <- cbind(mapData[matchIndexes[which(!is.na(matchIndexes))],2],tf.enh.motif.pairs[which(!is.na(matchIndexes)),2:ncol(tf.enh.motif.pairs)])
tf.enh.motif.pairs.yes <- tf.enh.motif.pairs[which(tf.enh.motif.pairs[,5]=="Yes"),]
tf.enh.motif.pairs.no <- tf.enh.motif.pairs[which(tf.enh.motif.pairs[,5]=="No"),]  

tf.enh.co.pairs <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_XW.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(paste(toupper(tf.enh.co.pairs[,5]),tf.enh.co.pairs[,1],sep="\t"),paste(toupper(tf.enh.motif.pairs.no[,1]),tf.enh.motif.pairs.no[,2],sep="\t"))
tf.enh.motif.pairs.no <- tf.enh.motif.pairs.no[matchIndexes[which(!is.na(matchIndexes))],]
tf.enh.motif.pairs <- rbind(tf.enh.motif.pairs.yes,tf.enh.motif.pairs.no)

tf.enh.co.pairs.all <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_XW_all.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(paste(toupper(tf.enh.co.pairs.all[,5]),tf.enh.co.pairs.all[,1],sep="\t"),paste(toupper(tf.enh.motif.pairs[,1]),tf.enh.motif.pairs[,2],sep="\t"))
tf.enh.co.pairs.all <- tf.enh.co.pairs.all[which(!is.na(matchIndexes)),]
tf.enh.motif.pairs <- tf.enh.motif.pairs[matchIndexes[which(!is.na(matchIndexes))],]
tf.enh.co.pairs.all[which(tf.enh.co.pairs.all[,4]==0),4] <- min(tf.enh.co.pairs.all[which(tf.enh.co.pairs.all[,4]>0),4])
tf.enh.motif.pairs <- cbind(tf.enh.motif.pairs[,c(2,1,3,4)],tf.enh.motif.pairs[,1])
tf.enh.motif.pairs[,2] <- tf.enh.co.pairs.all[,2]

row.indexes <- c()
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
geneexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/stringtiefile_XW/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(geneexpData) <- gsub("\\.\\d+","",geneexpData[,1])
geneexpData <- as.matrix(geneexpData[,3:ncol(geneexpData)])
for (rowindex in seq(1,nrow(tf.enh.motif.pairs))){
	rownum <- match(tf.enh.motif.pairs[rowindex,2],rownames(geneexpData))
	colnum <- match(tf.enh.motif.pairs[rowindex,4],colnames(geneexpData))
	if(geneexpData[tf.enh.motif.pairs[rowindex,2],tf.enh.motif.pairs[rowindex,4]] > 0){
		tf.enh.motif.pairs[rowindex,3] <- tf.enh.motif.pairs[rowindex,3]*geneexpData[tf.enh.motif.pairs[rowindex,2],tf.enh.motif.pairs[rowindex,4]]
		row.indexes <- c(row.indexes,rowindex)
	}
}
tf.enh.motif.pairs <- tf.enh.motif.pairs[row.indexes,]
write.table(tf.enh.motif.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile/TF_enhancer_network_final.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
