args <- commandArgs(trailingOnly = TRUE)
message(args[1],"\n")

library(foreach)
library(doParallel)
cl <- makeCluster(as.numeric(4))
registerDoParallel(cl)

cat.td.signal <- function (pos.vec, depth.mat){
	match.indexes <- match(seq(as.numeric(pos.vec[2])+1,as.numeric(pos.vec[3])),depth.mat[,2])
	relength <- (round(length(match.indexes)/3))*3
	pos.reads.count <- approx(depth.mat[match.indexes,3], n=relength)$y
	td.signal <- ((2*sum(pos.reads.count[((relength/3)+1):(2*(relength/3))])) + 1e-10)/((sum(pos.reads.count[1:(relength/3)])+sum(pos.reads.count[(2*(relength/3)+1):relength])) + 1e-10)
	return(td.signal)
}

gene.pos.data <- read.table(file="/public/ZhangJinLab/project_enhancer/known_genes.bed",header=FALSE,stringsAsFactors=FALSE)
gene.depth.data <- read.csv(file="/public/ZhangJinLab/project_enhancer/gene_position_count.csv",header=TRUE,stringsAsFactors=FALSE)
gene.depth.data[,3] <- rowSums(data.matrix(gene.depth.data[,3:ncol(gene.depth.data)]))
gene.depth.data[,4:ncol(gene.depth.data)] <- NULL
save(gene.pos.data,gene.depth.data,file="/public/ZhangJinLab/project_enhancer/known_genes_info.RData")
load("/public/ZhangJinLab/project_enhancer/known_genes_info.RData")

bg.td.signal.vec <- foreach(iterer=1:nrow(gene.pos.data),.combine=c,.multicombine=TRUE,.verbose=TRUE) %dopar% cat.td.signal(as.character(gene.pos.data[iterer,]),gene.depth.data)
bg.td.signal.thre <- quantile(bg.td.signal.vec,probs=0.95)
enh.pos.data <- read.table(file=args[1],header=FALSE,stringsAsFactors=FALSE)
enh.depth.data <- read.csv(file="/public/ZhangJinLab/project_enhancer/enhancer_position_count.csv",header=TRUE,stringsAsFactors=FALSE)
enh.depth.data[,3] <- rowSums(data.matrix(enh.depth.data[,3:ncol(enh.depth.data)]))
enh.depth.data[,4:ncol(enh.depth.data)] <- NULL
save(enh.depth.data,file="/public/ZhangJinLab/project_enhancer/enh_depth_info.RData")
enh.td.signal.vec <- foreach(iterer=1:nrow(enh.pos.data),.combine=c,.multicombine=TRUE,.verbose=TRUE) %dopar% cat.td.signal(as.character(enh.pos.data[iterer,]),enh.depth.data)
write.table(enh.pos.data[which(enh.td.signal.vec > bg.td.signal.thre),],file=paste(args[1],".bed",sep=""),row.names=FALSE,col.names=FALSE)

# Rscript --vanilla fig1a_step7_idnetify_two_direction_transcription.R /media/yuhua/yuhua_projects/enhProj/ENHData/filtered_mm10_td