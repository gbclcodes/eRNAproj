load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_pvalues_GSR.RData")
load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_corrs_GSR.RData")
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)

pvalueMatrix[which(is.na(pvalueMatrix),arr.ind=TRUE)] <- 1
pvalueMatrix <- t(apply(pvalueMatrix,1,p.adjust,"fdr"))
corrMatrix[which(is.na(corrMatrix),arr.ind=TRUE)] <- 0  
row.col.indexes <- which(corrMatrix > 0.7, arr.ind=TRUE)
sign.pairs <- cbind(rownames(mapInfo)[row.col.indexes[,1]],colnames(corrMatrix)[row.col.indexes[,2]],corrMatrix[row.col.indexes],pvalueMatrix[row.col.indexes])
sign.pairs <- sign.pairs[which(sign.pairs[,4] < 0.01),]
sign.pairs[,2] <- gsub("\\.\\d+","",sign.pairs[,2])

tf.genes <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/AnimalTFDB_mus_musculus_TF.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
match.indexes <- match(sign.pairs[,2],tf.genes[,3])
sign.tf.enh.pairs <- sign.pairs[which(!is.na(match.indexes)),]
sign.tf.enh.pairs <- cbind(sign.tf.enh.pairs,tf.genes[match.indexes[which(!is.na(match.indexes))],2])
write.table(sign.tf.enh.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_GSR.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

row.col.indexes <- which(corrMatrix > 0, arr.ind=TRUE)
sign.pairs <- cbind(rownames(mapInfo)[row.col.indexes[,1]],colnames(corrMatrix)[row.col.indexes[,2]],corrMatrix[row.col.indexes],pvalueMatrix[row.col.indexes])
sign.pairs[,2] <- gsub("\\.\\d+","",sign.pairs[,2])

tf.genes <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/AnimalTFDB_mus_musculus_TF.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
match.indexes <- match(sign.pairs[,2],tf.genes[,3])
sign.tf.enh.pairs <- sign.pairs[which(!is.na(match.indexes)),]
sign.tf.enh.pairs <- cbind(sign.tf.enh.pairs,tf.genes[match.indexes[which(!is.na(match.indexes))],2])
write.table(sign.tf.enh.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_GSR_all.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_pvalues_XW.RData")
load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_corrs_XW.RData")
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)

pvalueMatrix[which(is.na(pvalueMatrix),arr.ind=TRUE)] <- 1
pvalueMatrix <- t(apply(pvalueMatrix,1,p.adjust,"fdr"))
corrMatrix[which(is.na(corrMatrix),arr.ind=TRUE)] <- 0  
row.col.indexes <- which(corrMatrix > 0.7, arr.ind=TRUE)
sign.pairs <- cbind(rownames(mapInfo)[row.col.indexes[,1]],colnames(corrMatrix)[row.col.indexes[,2]],corrMatrix[row.col.indexes],pvalueMatrix[row.col.indexes])
sign.pairs <- sign.pairs[which(sign.pairs[,4] < 0.01),]
sign.pairs[,2] <- gsub("\\.\\d+","",sign.pairs[,2])

tf.genes <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/AnimalTFDB_mus_musculus_TF.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
match.indexes <- match(sign.pairs[,2],tf.genes[,3])
sign.tf.enh.pairs <- sign.pairs[which(!is.na(match.indexes)),]
sign.tf.enh.pairs <- cbind(sign.tf.enh.pairs,tf.genes[match.indexes[which(!is.na(match.indexes))],2])
write.table(sign.tf.enh.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_XW.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

row.col.indexes <- which(corrMatrix > 0, arr.ind=TRUE)
sign.pairs <- cbind(rownames(mapInfo)[row.col.indexes[,1]],colnames(corrMatrix)[row.col.indexes[,2]],corrMatrix[row.col.indexes],pvalueMatrix[row.col.indexes])
sign.pairs[,2] <- gsub("\\.\\d+","",sign.pairs[,2])

tf.genes <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/AnimalTFDB_mus_musculus_TF.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
match.indexes <- match(sign.pairs[,2],tf.genes[,3])
sign.tf.enh.pairs <- sign.pairs[which(!is.na(match.indexes)),]
sign.tf.enh.pairs <- cbind(sign.tf.enh.pairs,tf.genes[match.indexes[which(!is.na(match.indexes))],2])
write.table(sign.tf.enh.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/tf_enh_pairs_XW_all.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)