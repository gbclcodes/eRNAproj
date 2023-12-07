library(foreach)
library(doParallel)
cl <- makeCluster(as.numeric(4))
registerDoParallel(cl)

getPvalue <- function(x,y){
	res <- cor.test(x,y,method="pearson")
	return(res$p.value)
}

getCorrValue <- function(x,y){
	res <- cor.test(x,y,method="pearson")
	return(res$estimate)
}

evalPvalue <- function(vec, mat){
	return(apply(mat,1,getPvalue,vec))
}

evalCorr <- function(vec, mat){
	return(apply(mat,1,getCorrValue,vec))
}

geneexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/stringtiefile_GSR/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(geneexpData) <- geneexpData[,1]
geneexpData <- as.matrix(geneexpData[,3:ncol(geneexpData)])
geneexpData <- geneexpData[,c(1,2,3,6,7,8)]
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
enhexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GSRData/RNAseq_GSR/deepTools/enhancer_TPM.tab",sep=" ",header=TRUE,stringsAsFactors=FALSE)
rownames(enhexpData) <- rownames(mapInfo)
enhexpData <- enhexpData[,c(1,2,3,6,7,8)]

pvalueMatrix <- foreach(iterer=1:nrow(enhexpData),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% evalPvalue(as.numeric(enhexpData[iterer,]),geneexpData)
save(pvalueMatrix,file="/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_pvalues_GSR.RData")
corrMatrix <- foreach(iterer=1:nrow(enhexpData),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% evalCorr(as.numeric(enhexpData[iterer,]),geneexpData)
save(corrMatrix,file="/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_corrs_GSR.RData")

geneexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/stringtiefile_XW/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(geneexpData) <- geneexpData[,1]
geneexpData <- as.matrix(geneexpData[,3:ncol(geneexpData)])
mapInfo <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_clean_td_0.01EPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
enhexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/XWData/RNAseqData_XW/deepTools/enhancer_TPM.tab",sep=" ",header=TRUE,stringsAsFactors=FALSE)
rownames(enhexpData) <- rownames(mapInfo)

pvalueMatrix <- foreach(iterer=1:nrow(enhexpData),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% evalPvalue(as.numeric(enhexpData[iterer,]),geneexpData)
save(pvalueMatrix,file="/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_pvalues_XW.RData")
corrMatrix <- foreach(iterer=1:nrow(enhexpData),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% evalCorr(as.numeric(enhexpData[iterer,]),geneexpData)
save(corrMatrix,file="/media/yuhua/yuhua_projects/enhProj/ENHData/all_gene_enhancer_pair_corrs_XW.RData")