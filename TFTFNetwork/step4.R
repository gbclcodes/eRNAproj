library(foreach)
library(doParallel)
library(reshape2)
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
maxVec <- apply(plotMat,1,which.max)
write.table(cbind(names(maxVec),colnames(plotMat)[maxVec]),file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tfnames_to_stages.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

TFtoStageMapDataGSR <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/tfnames_to_stages.txt",stringsAsFactors=FALSE,header=FALSE)
TFtoStageMapDataXW <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tfnames_to_stages.txt",stringsAsFactors=FALSE,header=FALSE)

OocyteCommon <- intersect(TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1],TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1])
X2cellCommon <- intersect(TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C" | TFtoStageMapDataXW[,2]=="L2C"),1],TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1])
X4cellCommon <- intersect(TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1],TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1])
X8cellCommon <- intersect(TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1],TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1])
ICMCommon <- intersect(TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1],TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1])
TFcommon <- c(OocyteCommon,X2cellCommon,X4cellCommon,X8cellCommon,ICMCommon)
write.csv(TFcommon,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/commonTFs.csv",quote=FALSE,row.names=FALSE,col.names=FALSE)

TFtoModuleMapDataXW <- read.csv(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tfnames_in_modules.csv",stringsAsFactors=FALSE,header=FALSE)
maxVec <- apply(plotMat,1,max)

matchIndexesM1XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M1"),1],names(maxVec))
stageVecM1XW <- maxVec[matchIndexesM1XW]
stageVecM1XW <- sort(stageVecM1XW,decreasing=TRUE)

M1XWCommon <- intersect(names(stageVecM1XW),TFcommon)
M1XWSpecific <- setdiff(names(stageVecM1XW),TFcommon)

length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM1XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))

matchIndexesM2XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M2"),1],names(maxVec))
stageVecM2XW <- maxVec[matchIndexesM2XW]
stageVecM2XW <- sort(stageVecM2XW,decreasing=TRUE)

M2XWCommon <- intersect(names(stageVecM2XW),TFcommon)
M2XWSpecific <- setdiff(names(stageVecM2XW),TFcommon)

length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM2XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))

matchIndexesM3XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M3"),1],names(maxVec))
stageVecM3XW <- maxVec[matchIndexesM3XW]
stageVecM3XW <- sort(stageVecM3XW,decreasing=TRUE)

M3XWCommon <- intersect(names(stageVecM3XW),TFcommon)
M3XWSpecific <- setdiff(names(stageVecM3XW),TFcommon)

length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM3XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))


matchIndexesM4XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M4"),1],names(maxVec))
stageVecM4XW <- maxVec[matchIndexesM4XW]
stageVecM4XW <- sort(stageVecM4XW,decreasing=TRUE)

M4XWCommon <- intersect(names(stageVecM4XW),TFcommon)
M4XWSpecific <- setdiff(names(stageVecM4XW),TFcommon)

length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM4XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))


matchIndexesM5XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M5"),1],names(maxVec))
stageVecM5XW <- maxVec[matchIndexesM5XW]
stageVecM5XW <- sort(stageVecM5XW,decreasing=TRUE)

M5XWCommon <- intersect(names(stageVecM5XW),TFcommon)
M5XWSpecific <- setdiff(names(stageVecM5XW),TFcommon)

length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM5XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))


matchIndexesM6XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M6"),1],names(maxVec))
stageVecM6XW <- maxVec[matchIndexesM6XW]
stageVecM6XW <- sort(stageVecM6XW,decreasing=TRUE)

M6XWCommon <- intersect(names(stageVecM6XW),TFcommon)
M6XWSpecific <- setdiff(names(stageVecM6XW),TFcommon)

length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM6XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))


matchIndexesM7XW <- match(TFtoModuleMapDataXW[which(TFtoModuleMapDataXW[,2]=="M7"),1],names(maxVec))
stageVecM7XW <- maxVec[matchIndexesM7XW]
stageVecM7XW <- sort(stageVecM7XW,decreasing=TRUE)

M7XWCommon <- intersect(names(stageVecM7XW),TFcommon)
M7XWSpecific <- setdiff(names(stageVecM7XW),TFcommon)

length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="Zygote"),1]))
length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C"),1]))
length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="L2C"),1]))
length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1]))
length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1]))
length(intersect(names(stageVecM7XW),TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1]))
