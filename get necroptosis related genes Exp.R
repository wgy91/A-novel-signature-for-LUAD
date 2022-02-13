
library(limma)      
setwd("E:\\necroptosis")   


rt=read.table("symbol.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


gene=read.table("gene.txt", header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]


out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="tcga.necroptosisExp.txt",sep="\t",quote=F,col.names=F)


