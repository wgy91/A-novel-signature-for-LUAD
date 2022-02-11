

library(Rtsne)
library(ggplot2)
setwd("E:\\文章\\文章\\坏死+肺腺癌\\132pyroptosis\\25.PCA")   

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){

	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=rt[c(3:(ncol(rt)-2))]
	risk=rt[,"risk"]


	data.pca=prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)
	PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	

	pdf(file=pcaFile, height=4.5, width=5.5)      
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
		scale_colour_manual(name="Risk",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
	


	tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
	tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	

	pdf(file=tsneFile, height=4.5, width=5.5)       #保存输入出文件
	p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = risk)) +
		scale_colour_manual(name="Risk",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
}
bioPCA(inputFile="trainRisk.txt", pcaFile="train.PCA.pdf", tsneFile="train.t-SNE.pdf")
bioPCA(inputFile="testRisk.txt", pcaFile="test.PCA.pdf", tsneFile="test.t-SNE.pdf")

