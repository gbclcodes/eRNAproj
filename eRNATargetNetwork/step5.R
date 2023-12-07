high.enhancers.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_highscores_GSR_subset.txt",header=FALSE,stringsAsFactors=FALSE)
median.enhancers.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_medianscores_GSR_subset.txt",header=FALSE,,stringsAsFactors=FALSE)
low.enhancers.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_lowscores_GSR_subset.txt",header=FALSE,,stringsAsFactors=FALSE)

netData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_GSR_all.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(netData[,1],enhancer.stages[,1])
netData <- cbind(netData[,c(1,2,3,6,7)],enhancer.stages[matchIndexes,3])
write.table(netData,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_GSR_final.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

high.enhancers.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_highscores_XW.txt",header=FALSE,stringsAsFactors=FALSE)
median.enhancers.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_medianscores_XW.txt",header=FALSE,,stringsAsFactors=FALSE)
low.enhancers.gsr <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_lowscores_XW.txt",header=FALSE,,stringsAsFactors=FALSE)
enhancer.stages <- rbind(high.enhancers.gsr,median.enhancers.gsr,low.enhancers.gsr)

netData <- read.table(file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_XW_all.txt",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(netData[,1],enhancer.stages[,1])
netData <- cbind(netData[,c(1,2,3,6,7)],enhancer.stages[matchIndexes,3])
write.table(netData,file="/media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_XW_final.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)