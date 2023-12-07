library(foreach)
library(doParallel)
library(reshape2)
cl <- makeCluster(as.numeric(20))
registerDoParallel(cl)

expValNormalization <- function(expVec){
	expVec <- log2(expVec+1)/sum(log2(expVec+1))
	return(expVec)
}

maxvalue <- function(x){
    max_value <- max(x)
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

expValNormalization <- function(expVec){
	expVec <- log2(expVec+1)/sum(log2(expVec+1))
	return(expVec)
}

tf.genes <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/AnimalTFDB_mus_musculus_TF.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
tf.genes <- tf.genes[which(!duplicated(tf.genes[,2])),]
mapData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE_ID_to_NAME.txt",header=FALSE,stringsAsFactors=FALSE)
mapData <- mapData[which(!duplicated(mapData[,2])),]
matchIndexes <- match(toupper(mapData[,2]),toupper(tf.genes[,2]))
mapData <- mapData[which(!is.na(matchIndexes)),]
tf.genes <- tf.genes[matchIndexes[which(!is.na(matchIndexes))],]

geneexpData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/GENEData/stringtiefile_GSR/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(geneexpData) <- geneexpData[,1]
geneexpData <- as.matrix(geneexpData[,3:ncol(geneexpData)])
rownames(geneexpData) <- gsub("\\.\\d+","",rownames(geneexpData))
matchIndexes <- match(rownames(geneexpData),tf.genes[,3])
geneexpData <- geneexpData[which(!is.na(matchIndexes)),]
tf.genes <- tf.genes[matchIndexes[which(!is.na(matchIndexes))],]
mapData <- mapData[matchIndexes[which(!is.na(matchIndexes))],]
rownames(geneexpData) <- toupper(tf.genes[,2])

files <- Sys.glob(file.path("/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile","*","*","ame.tsv"))
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
enscoreMat <- enscoreMat[c(matchIndexes,234),]
geneexpData <- geneexpData[c(matchIndexes,234),]

stageNames <- c("MIIOocyte","X2cell","X4cell","X8cell","morula","ICM")
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
write.table(cbind(names(maxVec),colnames(plotMat)[maxVec]),file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/tfnames_to_stages.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

TFtoStageMapDataGSR <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/tfnames_to_stages.txt",stringsAsFactors=FALSE,header=FALSE)
TFtoStageMapDataXW <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/tfnames_to_stages.txt",stringsAsFactors=FALSE,header=FALSE)

OocyteCommon <- intersect(TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1],TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="MIIOocyte"),1])
X2cellCommon <- intersect(TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1],TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="E2C" | TFtoStageMapDataXW[,2]=="L2C"),1])
X4cellCommon <- intersect(TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1],TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M4C"),1])
X8cellCommon <- intersect(TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1],TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="M8C"),1])
ICMCommon <- intersect(TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1],TFtoStageMapDataXW[which(TFtoStageMapDataXW[,2]=="ICM"),1])
TFcommon <- c(OocyteCommon,X2cellCommon,X4cellCommon,X8cellCommon,ICMCommon)
write.csv(TFcommon,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/commonTFs.csv",quote=FALSE,row.names=FALSE,col.names=FALSE)

TFtoModuleMapDataGSR <- read.csv(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/tfnames_in_modules.csv",stringsAsFactors=FALSE,header=FALSE)
maxVec <- apply(plotMat,1,max)

matchIndexesM1GSR <- match(TFtoModuleMapDataGSR[which(TFtoModuleMapDataGSR[,2]=="M1"),1],names(maxVec))
stageVecM1GSR <- maxVec[matchIndexesM1GSR]
stageVecM1GSR <- sort(stageVecM1GSR,decreasing=TRUE)

M1GSRCommon <- intersect(names(stageVecM1GSR),TFcommon)
M1GSRSpecific <- setdiff(names(stageVecM1GSR),TFcommon)

length(intersect(names(stageVecM1GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM1GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1]))
length(intersect(names(stageVecM1GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1]))
length(intersect(names(stageVecM1GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1]))
length(intersect(names(stageVecM1GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="morula"),1]))
length(intersect(names(stageVecM1GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1]))

matchIndexesM2GSR <- match(TFtoModuleMapDataGSR[which(TFtoModuleMapDataGSR[,2]=="M2"),1],names(maxVec))
stageVecM2GSR <- maxVec[matchIndexesM2GSR]
stageVecM2GSR <- sort(stageVecM2GSR,decreasing=TRUE)

M2GSRCommon <- intersect(names(stageVecM2GSR),TFcommon)
M2GSRSpecific <- setdiff(names(stageVecM2GSR),TFcommon)

length(intersect(names(stageVecM2GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM2GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1]))
length(intersect(names(stageVecM2GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1]))
length(intersect(names(stageVecM2GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1]))
length(intersect(names(stageVecM2GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="morula"),1]))
length(intersect(names(stageVecM2GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1]))

matchIndexesM3GSR <- match(TFtoModuleMapDataGSR[which(TFtoModuleMapDataGSR[,2]=="M3"),1],names(maxVec))
stageVecM3GSR <- maxVec[matchIndexesM3GSR]
stageVecM3GSR <- sort(stageVecM3GSR,decreasing=TRUE)

M3GSRCommon <- intersect(names(stageVecM3GSR),TFcommon)
M3GSRSpecific <- setdiff(names(stageVecM3GSR),TFcommon)

length(intersect(names(stageVecM3GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM3GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1]))
length(intersect(names(stageVecM3GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1]))
length(intersect(names(stageVecM3GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1]))
length(intersect(names(stageVecM3GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="morula"),1]))
length(intersect(names(stageVecM3GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1]))

matchIndexesM4GSR <- match(TFtoModuleMapDataGSR[which(TFtoModuleMapDataGSR[,2]=="M4"),1],names(maxVec))
stageVecM4GSR <- maxVec[matchIndexesM4GSR]
stageVecM4GSR <- sort(stageVecM4GSR,decreasing=TRUE)

M4GSRCommon <- intersect(names(stageVecM4GSR),TFcommon)
M4GSRSpecific <- setdiff(names(stageVecM4GSR),TFcommon)

length(intersect(names(stageVecM4GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM4GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1]))
length(intersect(names(stageVecM4GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1]))
length(intersect(names(stageVecM4GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1]))
length(intersect(names(stageVecM4GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="morula"),1]))
length(intersect(names(stageVecM4GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1]))

matchIndexesM5GSR <- match(TFtoModuleMapDataGSR[which(TFtoModuleMapDataGSR[,2]=="M5"),1],names(maxVec))
stageVecM5GSR <- maxVec[matchIndexesM5GSR]
stageVecM5GSR <- sort(stageVecM5GSR,decreasing=TRUE)

M5GSRCommon <- intersect(names(stageVecM5GSR),TFcommon)
M5GSRSpecific <- setdiff(names(stageVecM5GSR),TFcommon)

length(intersect(names(stageVecM5GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="MIIOocyte"),1]))
length(intersect(names(stageVecM5GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X2cell"),1]))
length(intersect(names(stageVecM5GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X4cell"),1]))
length(intersect(names(stageVecM5GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="X8cell"),1]))
length(intersect(names(stageVecM5GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="morula"),1]))
length(intersect(names(stageVecM5GSR),TFtoStageMapDataGSR[which(TFtoStageMapDataGSR[,2]=="ICM"),1]))
