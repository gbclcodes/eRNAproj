library(foreach)
library(doParallel)
library(reshape2)
cl <- makeCluster(as.numeric(20))
registerDoParallel(cl)

expValNormalization <- function(expVec){
	expVec <- log2(expVec+1)/sum(log2(expVec+1))
	return(expVec)
}

JSscoreCal <- function(probVec,stdMat,stageNames){
	HscoreVec <- c()
	probVec[which(probVec < 2e-100)] <- runif(length(which(probVec < 2e-100)),min=1e-100,max=2e-100)
	for(colIndex in seq(1,ncol(stdMat))){
		meanProb <- (probVec+stdMat[,colIndex])/2
		meanH <- -sum(meanProb*log2(meanProb))
		pop1H <- -sum(probVec*log2(probVec))
		pop2H <- -sum(stdMat[,colIndex]*log2(stdMat[,colIndex]))
		HscoreVec <- c(HscoreVec,(1-sqrt(meanH-(pop1H+pop2H)/2)))
	}
	JSscore <- max(HscoreVec)
    stageName <- stageNames[which(HscoreVec==max(HscoreVec))]
	return(c(JSscore,stageName))
}

tf.genes <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/AnimalTFDB_mus_musculus_TF.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
tf.genes <- tf.genes[which(!duplicated(tf.genes[,2])),]
mapData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE_ID_to_NAME.txt",header=FALSE,stringsAsFactors=FALSE)
mapData <- mapData[which(!duplicated(mapData[,2])),]
matchIndexes <- match(toupper(mapData[,2]),toupper(tf.genes[,2]))
mapData <- mapData[which(!is.na(matchIndexes)),]
tf.genes <- tf.genes[matchIndexes[which(!is.na(matchIndexes))],]

geneexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/stringtiefile_XW/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(geneexpData) <- geneexpData[,1]
geneexpData <- as.matrix(geneexpData[,3:ncol(geneexpData)])
rownames(geneexpData) <- gsub("\\.\\d+","",rownames(geneexpData))
matchIndexes <- match(rownames(geneexpData),tf.genes[,3])
geneexpData <- geneexpData[which(!is.na(matchIndexes)),]
tf.genes <- tf.genes[matchIndexes[which(!is.na(matchIndexes))],]
mapData <- mapData[matchIndexes[which(!is.na(matchIndexes))],]
rownames(geneexpData) <- toupper(tf.genes[,2])

files <- Sys.glob(file.path("/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile","*","*","ame.tsv"))
enData <- read.delim(file=files[1],header=TRUE,stringsAsFactors=FALSE)
posIndex <- regexpr('tiffile.*ame.tsv',files[1])
sampleid <- substring(files[1],posIndex+8,posIndex+attr(posIndex,'match.length')-9)
enscoreMat <- c(unlist(strsplit(sampleid,"\\/")),(-log10(enData[1,7])))

for (filename in files[2:length(files)]){
	if(file.info(filename)[1] > 0){
		enData <- read.delim(file=filename,header=TRUE,stringsAsFactors=FALSE)
		posIndex <- regexpr('tiffile.*ame.tsv',filename)
		sampleid <- substring(filename,posIndex+8,posIndex+attr(posIndex,'match.length')-9)
		if(!is.null(enData[1,6])){
			enscoreVec <- c(unlist(strsplit(sampleid,"\\/")),(-log10(enData[1,7])))
		}else{
			enscoreVec <- c(unlist(strsplit(sampleid,"\\/")),0)
		}
		enscoreMat <- rbind(enscoreMat,enscoreVec)
	}
}
enscoreMat <- as.data.frame(enscoreMat)
colnames(enscoreMat) <- c("stageName","motifName","enScore")
enscoreMat <- acast(enscoreMat, motifName~stageName, value.var="enScore")
matchIndexes <- match(mapData[,1],rownames(enscoreMat))

enscoreMat <- enscoreMat[matchIndexes[which(!is.na(matchIndexes))],]
mapData <- mapData[which(!is.na(matchIndexes)),]
tf.genes <- tf.genes[which(!is.na(matchIndexes)),]
geneexpData <- geneexpData[which(!is.na(matchIndexes)),]

stageNames <- colnames(enscoreMat)
enscoreMat <- matrix(as.numeric(enscoreMat), ncol = ncol(enscoreMat))
rownames(enscoreMat) <- toupper(mapData[,2])
colnames(enscoreMat) <- stageNames

matchIndexes <- which(apply(enscoreMat,1,max) > 3 & rowSums(geneexpData) > 0)
enscoreMat <- enscoreMat[c(matchIndexes,214),]
geneexpData <- geneexpData[c(matchIndexes,214),]

stageNames <- c("MIIOocyte","Zygote","E2C","L2C","M4C","M8C","ICM")
plotMat <- enscoreMat[,stageNames]*geneexpData[,stageNames]
plotMat <- plotMat[which(rowSums(plotMat) > 0),]
tf.enh.motif.pairs.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile/TF_enhancer_network_final.txt",header=FALSE,stringsAsFactors=FALSE)
tf.enh.motif.pairs.gsr <- tf.enh.motif.pairs.gsr[,c(5,1,3,4)]
tf.enh.motif.pairs.gsr[,1] <- toupper(tf.enh.motif.pairs.gsr[,1])
tf.names <- unique(tf.enh.motif.pairs.gsr[,1])
matchIndexes <- match(rownames(plotMat),tf.names)
plotMat <- plotMat[which(!is.na(matchIndexes)),]
plotMat <- plotMat[which(rowSums(plotMat) > 0),]

bgMatrix <- matrix(0,nrow=length(stageNames),ncol=length(stageNames))
cat("The stage number:",length(stageNames),"\n")
for(stageIter in seq(1,length(stageNames))){
	randomNum <- runif(length(stageNames),min=1e-100,max=2e-100)
	bgMatrix[stageIter,stageIter] <- 1
	bgMatrix[,stageIter] <- bgMatrix[,stageIter] + randomNum
}

# JS score calculation
plotMatNorm <- t(apply(plotMat,1,expValNormalization))
tfJSscoreVec <- foreach(iter=1:nrow(plotMatNorm),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% JSscoreCal(plotMatNorm[iter,],bgMatrix,stageNames)
rownames(tfJSscoreVec) <- rownames(plotMatNorm)
enhnonaIndexes <- which(!is.na(tfJSscoreVec[,1]))
tfJSscoreVec <- tfJSscoreVec[enhnonaIndexes,]
tfJSscoreVec <- as.data.frame(tfJSscoreVec)
colnames(tfJSscoreVec) <- c("JSscore","StageName")
tfJSscoreVec[,2] <- factor(tfJSscoreVec[,2],levels=rev(c("MIIOocyte","Zygote","E2C","L2C","M4C","M8C","ICM")),ordered=TRUE)
tfJSscoreVec <- tfJSscoreVec[order(tfJSscoreVec[,2],tfJSscoreVec[,1],decreasing=TRUE),]
matchIndexes <- match(rownames(tfJSscoreVec),rownames(plotMat))
plotMat <- plotMat[matchIndexes,]

library(pheatmap)
library(wesanderson)
library(dendextend)
pdf(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_tf_simi_clustering_heatmap.pdf",width=40,height=4)
p <- pheatmap(t(plotMat),show_rownames=TRUE,scale="column",clustering_distance_cols="manhattan",cluster_rows=FALSE,cluster_cols=TRUE,fontsize=10,border_color="NA",color = colorRampPalette(c("#1E90FF","#FFFFF0","#FF00FF"))(100))
dev.off()

pdf(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_tf_simi_clustering_tree.pdf",width=40,height=4)
p.dend <- as.dendrogram(p$tree_col)
p.dend <- reorder(p.dend,seq(1,2160000,10000),agglo.FUN="mean")
p.dend %>% set("branches_k_col", k = 7) %>% plot()
dev.off()

pdf(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_tf_simi_clustering_heatmap.pdf",width=40,height=4)
p <- pheatmap(t(plotMat[labels(p.dend),]),show_rownames=TRUE,scale="column",cluster_rows=FALSE,cluster_cols=FALSE,fontsize=10,border_color="NA",color = colorRampPalette(c("#1E90FF","#FFFFF0","#FF00FF"))(100))
dev.off()

maxStage <- function(x,stageNames){
	return(stageNames[which(x==max(x))])
}
maxStages <- apply(plotMat[labels(p.dend),],1,maxStage,stageNames)
write.table(cbind(plotMat[labels(p.dend),],maxStages)[seq(nrow(plotMat),1,-1),],file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_reg_strength_mat.txt",row.names=TRUE,col.names=TRUE)

tfCorMat <- cor(t(plotMat))
pdf(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_tf_comodule_heatmap.pdf",width=40,height=40)
pheatmap(tfCorMat[labels(p.dend),labels(p.dend)],show_rownames=TRUE,scale="column",cluster_rows=FALSE,cluster_cols=FALSE,fontsize=10,border_color="NA",color = colorRampPalette(c("#1E90FF","#FFFFF0","#FF00FF"))(100))
dev.off()
write.table(rev(labels(p.dend)),file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tfnames_in_modules.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
