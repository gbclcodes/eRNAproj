files <- Sys.glob(file.path("/media/yuhua/yuhua_projects/enhProj/ENHData/htseqcountfile_GSR","*","htseq_count.txt"))
expData <- read.delim(file=files[1],header=FALSE,stringsAsFactors=FALSE)
expData <- expData[,c(1,2)]
posIndex <- regexpr('htseqcountfile.*htseq_count',files[1])
sampleid <- substring(files[1],posIndex+19,posIndex+attr(posIndex,'match.length')-13)
uniqIndex <- which(!duplicated(as.character(expData[,1])))
uniqEnsemids <- as.character(expData[uniqIndex,1])
expData <- expData[uniqIndex,]
columnNames <- c("gene_id",sampleid)
for (filename in files[2:length(files)]){
	if(file.info(filename)[1] > 0){
		sampleMat <- read.delim(file=filename,header=FALSE,stringsAsFactors=FALSE)
		matchIndex <- match(uniqEnsemids,sampleMat[,1])
		expData <- cbind(expData,sampleMat[matchIndex,2])
		posIndex <- regexpr('htseqcountfile.*htseq_count',filename)
		sampleid <- substring(filename,posIndex+19,posIndex+attr(posIndex,'match.length')-13)
		cat(sampleid,"\n")
		columnNames <- append(columnNames,sampleid)
	}
}
colnames(expData) <- columnNames
write.table(expData,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_GSR_counts.txt",sep='\t',quote=FALSE,row.names=T,col.names=T) 