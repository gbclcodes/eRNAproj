load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_pvalues_GSR.RData")
load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_corrs_GSR.RData")
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
enhexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GSRData/RNAseq_GSR/deepTools/enhancer_TPM.tab",sep=" ",header=TRUE,stringsAsFactors=FALSE)
rownames(enhexpData) <- rownames(mapInfo)
target.stages <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/gene_stage_specific_scores_GSR.txt",header=FALSE,stringsAsFactors=FALSE)

pvalueMatrix[which(is.na(pvalueMatrix),arr.ind=TRUE)] <- 1
pvalueMatrix <- t(apply(pvalueMatrix,1,p.adjust,"fdr"))
corrMatrix[which(is.na(corrMatrix),arr.ind=TRUE)] <- 0  
row.col.indexes <- which(corrMatrix > 0.3, arr.ind=TRUE)
matchIndexes <- match(colnames(corrMatrix)[row.col.indexes[,2]],target.stages[,1])
matchIndexes <- match(target.stages[matchIndexes,3],colnames(enhexpData))

sign.pairs <- cbind(rownames(mapInfo)[row.col.indexes[,1]],colnames(corrMatrix)[row.col.indexes[,2]],corrMatrix[row.col.indexes],pvalueMatrix[row.col.indexes],as.numeric(corrMatrix[row.col.indexes])*(-log10(as.numeric(pvalueMatrix[row.col.indexes])+1e-100)),as.numeric(corrMatrix[row.col.indexes])*enhexpData[cbind(row.col.indexes[,1],matchIndexes)])
sign.pairs <- sign.pairs[order(as.numeric(sign.pairs[,5]),decreasing=TRUE),]
sign.pairs <- sign.pairs[1:round((11442*53715)*0.05),]
sign.pairs[,2] <- gsub("\\.\\d+","",sign.pairs[,2])
write.table(sign.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_coexpr_pairs_GSR_all.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_pvalues_XW.RData")
load("/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_corrs_XW.RData")
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
enhexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/XWData/RNAseqData_XW/deepTools/enhancer_TPM.tab",sep=" ",header=TRUE,stringsAsFactors=FALSE)
rownames(enhexpData) <- rownames(mapInfo)
target.stages <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/gene_stage_specific_scores_XW.txt",header=FALSE,stringsAsFactors=FALSE)

pvalueMatrix[which(is.na(pvalueMatrix),arr.ind=TRUE)] <- 1
pvalueMatrix <- t(apply(pvalueMatrix,1,p.adjust,"fdr"))
corrMatrix[which(is.na(corrMatrix),arr.ind=TRUE)] <- 0  
row.col.indexes <- which(corrMatrix > 0.3, arr.ind=TRUE)
matchIndexes <- match(colnames(corrMatrix)[row.col.indexes[,2]],target.stages[,1])
matchIndexes <- match(target.stages[matchIndexes,3],colnames(enhexpData))

sign.pairs <- cbind(rownames(mapInfo)[row.col.indexes[,1]],colnames(corrMatrix)[row.col.indexes[,2]],corrMatrix[row.col.indexes],pvalueMatrix[row.col.indexes],as.numeric(corrMatrix[row.col.indexes])*(-log10(as.numeric(pvalueMatrix[row.col.indexes])+1e-100)),as.numeric(corrMatrix[row.col.indexes])*enhexpData[cbind(row.col.indexes[,1],matchIndexes)])
sign.pairs <- sign.pairs[order(as.numeric(sign.pairs[,5]),decreasing=TRUE),]
sign.pairs <- sign.pairs[1:round((8471*53715)*0.05),]
sign.pairs[,2] <- gsub("\\.\\d+","",sign.pairs[,2])
write.table(sign.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_coexpr_pairs_XW_all.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)