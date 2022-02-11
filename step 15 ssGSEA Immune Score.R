
library(limma)
library(reshape2)
library(ggpubr)
setwd("E:\\ÎÄÕÂ\\ÎÄÕÂ\\»µËÀ+·ÎÏÙ°©\\132pyroptosis\\33.scoreCor\\test\\22")      


scoreCor=function(riskFile=null, scoreFile=null, project=null){

	data=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(data)

	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	

	sameSample=intersect(row.names(data),row.names(risk))
	data=data[sameSample,,drop=F]
	risk=risk[sameSample,,drop=F]
	rt=cbind(data,risk[,c("riskScore","risk")])
	rt=rt[,-(ncol(rt)-1)]
	

	immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages M1","Macrophages M2",
	          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
	          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
	rt1=rt[,c(immCell,"risk")]
	data=melt(rt1,id.vars=c("risk"))
	colnames(data)=c("Risk","Type","Score")
	data$Risk=factor(data$Risk, levels=c("low","high"))
	p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     	xlab="",ylab="Score",add = "none",palette = c("blue","red") )
	p=p+rotate_x_text(50)
	p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

	pdf(file=paste0(project,".immCell.pdf"), width=7, height=6)
	print(p)
	dev.off()
	

	immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
	          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
	          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
	          "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
	rt1=rt[,c(immFunction,"risk")]
	data=melt(rt1,id.vars=c("risk"))
	colnames(data)=c("Risk","Type","Score")
	data$Risk=factor(data$Risk, levels=c("low","high"))
	p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     	xlab="",ylab="Score",add = "none",palette = c("blue","red") )
	p=p+rotate_x_text(50)
	p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

	pdf(file=paste0(project,".immFunction.pdf"), width=7, height=6)
	print(p)
	dev.off()
}

scoreCor(riskFile="trainRisk.txt", scoreFile="TCGA.score.txt", project="TCGA")
scoreCor(riskFile="testRisk.txt", scoreFile="GEO.score.txt", project="GEO")


