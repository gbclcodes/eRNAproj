library(foreach)
library(doParallel)
library(reshape2)
library(Matrix)
library(igraph)
library(bigSCale)
cl <- makeCluster(as.numeric(20))
registerDoParallel(cl)

evalCorr <- function(vec, mat){
	return(apply(mat,1,getCorrValue,vec))
}

getCorrValue <- function(x,y){
	res <- cor.test(x,y,method="pearson")
	return(res$estimate)
}

evalPvalue <- function(vec, mat){
	return(apply(mat,1,getPvalue,vec))
}

getPvalue <- function(x,y){
	res <- cor.test(x,y,method="pearson")
	return(res$p.value)
}

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

coRePvalMatrix <- foreach(iterer=1:nrow(plotMat),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% evalPvalue(as.numeric(plotMat[iterer,]),plotMat)
coRePvalMatrix  <- t(apply(coRePvalMatrix,1,p.adjust,"fdr"))
coReCorrMatrix <- foreach(iterer=1:nrow(plotMat),.combine=rbind,.multicombine=TRUE,.verbose=TRUE) %dopar% evalCorr(as.numeric(plotMat[iterer,]),plotMat)
coRePvalMatrix[lower.tri(coRePvalMatrix,diag=TRUE)] <- 1
coReCorrMatrix[lower.tri(coReCorrMatrix,diag=TRUE)] <- 0

row.col.indexes <- which(coRePvalMatrix < 0.05 & coReCorrMatrix > 0.3, arr.ind=TRUE)
sign.pairs <- cbind(rownames(plotMat)[row.col.indexes[,1]],rownames(plotMat)[row.col.indexes[,2]],-log10(coRePvalMatrix[row.col.indexes])*coReCorrMatrix[row.col.indexes])
colnames(sign.pairs) <- c("TF1","TF2","coreg.score")
write.table(sign.pairs,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_tf_motif_coregnetwork.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)


knownTFs <- as.character(read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/knownTFs.csv",sep=",",header=FALSE,stringsAsFactors=FALSE)[,1])
commonTFs <- as.character(read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/commonTFs.csv",sep=",",header=FALSE,stringsAsFactors=FALSE)[,1])

tf.tf.pairs.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tf_tf_motif_coregnetwork.txt",header=TRUE,stringsAsFactors=FALSE)
tf.tf.pairs.gsr[which(is.infinite(tf.tf.pairs.gsr[,3])),3] <- 60
tf.module.names <- unique(c(tf.tf.pairs.gsr[,1],tf.tf.pairs.gsr[,2]))

tf.class <- rep("Other",length(tf.module.names))
match.known <- match(knownTFs,tf.module.names)
match.common <- match(commonTFs,tf.module.names)
tf.class[match.common] <- "Common"
tf.class[match.known] <- "Known"

tf.to.stage <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tfnames_to_stages.txt",header=FALSE,stringsAsFactors=FALSE)
match.indexes <- match(tf.module.names,tf.to.stage[,1])
tf.stage <- tf.to.stage[match.indexes,2]

tf.row.indexes <- match(tf.tf.pairs.gsr[,1],tf.module.names)
tf.col.indexes <- match(tf.tf.pairs.gsr[,2],tf.module.names)

tf.tf.sp.mat <- sparseMatrix(i=tf.row.indexes,j=tf.col.indexes,x=as.numeric(tf.tf.pairs.gsr[,3]),symmetric=TRUE,dims=c(length(tf.module.names),length(tf.module.names)),dimnames=list(c(tf.module.names),c(tf.module.names)))
tf.tf.net.graph <- graph.adjacency(tf.tf.sp.mat,mode="undirected",weighted=TRUE)

tf.tf.net.graph <- set_vertex_attr(tf.tf.net.graph, V(tf.tf.net.graph), name = "stage", value = tf.stage)
tf.tf.net.graph <- set_vertex_attr(tf.tf.net.graph, V(tf.tf.net.graph), name = "class", value = tf.class)
toCytoscape(G = tf.tf.net.graph,file.name = '/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/TF_TF_motif_coregnetwork.json')