files <- Sys.glob(file.path("/media/yuhua/yuhua_projects/enhProj/ENHData/depthfiles_GSR","*","depth.txt"))
expData <- read.delim(file=files[1],header=FALSE,stringsAsFactors=FALSE)
posIndex <- regexpr('GSR\\/.*\\/depth.txt',files[1])
sampleid <- substring(files[1],posIndex+4,posIndex+attr(posIndex,'match.length')-11)
columnNames <- c("chr","pos",sampleid)
for (filename in files[2:length(files)]){
	if(file.info(filename)[1] > 0){
		sampleMat <- read.delim(file=filename,header=FALSE,stringsAsFactors=FALSE)
		expData <- cbind(expData,sampleMat[,3])
		posIndex <- regexpr('GSR\\/.*\\/depth.txt',filename)
		sampleid <- substring(filename,posIndex+4,posIndex+attr(posIndex,'match.length')-11)
		cat(sampleid,"\n")
		columnNames <- append(columnNames,sampleid)
	}
}
colnames(expData) <- columnNames
write.csv(expData,file="/media/yuhua/yuhua_projects/enhProj/ENHData/depthfiles_GSR/enhancer_position_count.csv",sep=',',quote=FALSE,row.names=F,col.names=T)


files <- Sys.glob(file.path("/media/yuhua/yuhua_projects/enhProj/GENEData/depthfiles_GSR","*","depth.txt"))
expData <- read.delim(file=files[1],header=FALSE,stringsAsFactors=FALSE)
posIndex <- regexpr('\\GSR/.*\\/depth.txt',files[1])
sampleid <- substring(files[1],posIndex+4,posIndex+attr(posIndex,'match.length')-11)
columnNames <- c("chr","pos",sampleid)
for (filename in files[2:length(files)]){
	if(file.info(filename)[1] > 0){
		sampleMat <- read.delim(file=filename,header=FALSE,stringsAsFactors=FALSE)
		expData <- cbind(expData,sampleMat[,3])
		posIndex <- regexpr('\\GSR/.*\\/depth.txt',filename)
		sampleid <- substring(filename,posIndex+4,posIndex+attr(posIndex,'match.length')-11)
		cat(sampleid,"\n")
		columnNames <- append(columnNames,sampleid)
	}
}
colnames(expData) <- columnNames
write.csv(expData,file="/media/yuhua/yuhua_projects/enhProj/GENEData/depthfiles_GSR/gene_position_count.csv",sep=',',quote=FALSE,row.names=F,col.names=T) 
