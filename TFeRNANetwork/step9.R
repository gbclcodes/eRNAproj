netData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile/TF_enhancer_network_final.txt",sep=" ",header=FALSE,stringsAsFactors=FALSE)
outData <- c()
tfIDs <- unique(netData[,2])
for (tfID in tfIDs){
	tfData <- netData[which(netData[,2]==tfID),]
	tfData <- tfData[order(tfData[,3],decreasing=TRUE),]
	if(nrow(tfData) >= 100){
		if(length(outData)==0){
			outData <- tfData[1:round(nrow(tfData)*0.05),]
		}else{
			outData <- rbind(outData,tfData[1:round(nrow(tfData)*0.05),])
		}
	}else if(nrow(tfData) < 100 && nrow(tfData) >= 5){
		if(length(outData)==0){
			outData <- tfData[1:5,]
		}else{
			outData <- rbind(outData,tfData[1:5,])
		}
	}else{
		if(length(outData)==0){
			outData <- tfData
		}else{
			outData <- rbind(outData,tfData)
		}
	}
}
write.table(outData,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile/TF_enhancer_network_core.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

netData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile/TF_enhancer_network_final.txt",sep=" ",header=FALSE,stringsAsFactors=FALSE)
outData <- c()
tfIDs <- unique(netData[,2])
for (tfID in tfIDs){
	tfData <- netData[which(netData[,2]==tfID),]
	tfData <- tfData[order(tfData[,3],decreasing=TRUE),]
	if(nrow(tfData) >= 100){
		if(length(outData)==0){
			outData <- tfData[1:round(nrow(tfData)*0.05),]
		}else{
			outData <- rbind(outData,tfData[1:round(nrow(tfData)*0.05),])
		}
	}else if(nrow(tfData) < 100 && nrow(tfData) >= 5){
		if(length(outData)==0){
			outData <- tfData[1:5,]
		}else{
			outData <- rbind(outData,tfData[1:5,])
		}
	}else{
		if(length(outData)==0){
			outData <- tfData
		}else{
			outData <- rbind(outData,tfData)
		}
	}
}
write.table(outData,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile/TF_enhancer_network_core.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
