rm (list=ls())

setwd("/psc/home/zhaocheng/My_project/Dongmeng_SRDX")

wt1=read.table("tmp_data/W1_RMAP/htseq_out.R")
wt2=read.table("tmp_data/W2_RMAP/htseq_out.R")
mt1=read.table("tmp_data/T1_RMAP/htseq_out.R")
mt2=read.table("tmp_data/T2_RMAP/htseq_out.R")

all=data.frame(wt1=wt1$V2,wt2=wt2$V2,mt1=mt1$V2,mt2=mt2$V2)
rownames(all)=wt1$V1

group= c(rep("col0",2),rep("srdx",2))
library(edgeR)
cds = DGEList(all,group =group)
cds = calcNormFactors(cds)
keep = rowSums(cpm(cds)>1)>=2
cds =cds[keep,]
log2counts = log2(cds$counts+1)
pdf ("results/col0_srdx_correlation.pdf")
panel.cor <- function(x,y, ...)
{ 
par(usr=c(0,1,0,1))
txt <- as.character(format(cor(x,y,method=c("spearman")),digits=4))
text(0.5,0.5,txt,cex =2*abs(cor(x,y)))
}
pairs(log2counts[,1:4],upper.panel=panel.cor,,main="Relationship of col0 replicates and srdx replicates in Arabidopsis")
dev.off()

cds=estimateCommonDisp(cds,verbose=TRUE)
cds = estimateTagwiseDisp(cds)
cpmcounts = cpm(cds$counts+1)
tcounts = t(cpmcounts)
pca.total = prcomp(log2(tcounts), retx=TRUE)
pdf("results/col0_srdx_pca.pdf")
c=round(100*summary(pca.total)$importance[2,1],digits=2)
d=round(100*summary(pca.total)$importance[2,2],digits=2)
tt=c("col0_1","col0_2","srdx_1","srdx_2")
plot(pca.total$x[,1:2], pch=c(15:16), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),col=c(rep("black",2),rep("red",2)),main="PCA Plot of Samples")
legend("bottomright",cex=0.6,border=F, legend=tt,pch=c(15:16), col=c(rep("black",2),rep("red",2)),bty="n")
dev.off()
de.out=exactTest(cds,pair=c("col0","srdx"))
rawcount = as.matrix(all)[(rownames(topTags(de.out,n=100000))),]
rawcount=rawcount[order(rownames(rawcount)),]
table=topTags(de.out,n=100000)$table
table=table[order(rownames(table)),]
sam_cpm=cpm(cds$count)[row.names(rawcount),]
all=cbind(rawcount,sam_cpm,table)

all=all[order(all$PValue),]
write.table(all,"results/col0_and_srdx.edgeR.out",sep="\t",row.names=T,col.names=F,quote=F)

de.up=subset(all,((all$FDR<0.05) & (all$logFC) >=1))
de.down=subset(all,((all$FDR<0.05) & (all$logFC) <= -1))
write.table(de.up,"results/DEG.up",sep="\t",row.names=T,col.names=F,quote=F)
write.table(de.down,"results/DEG.down",sep="\t",row.names=T,col.names=F,quote=F)


c1up=nrow(subset(all,((all$FDR<0.05) & (all$logFC) >=1)))
c1down=nrow(subset(all,((all$FDR<0.05) & (all$logFC) <= -1)))
c0.58up=nrow(subset(all,((all$FDR<0.05) & (all$logFC) >=0.58)))
c0.58down=nrow(subset(all,((all$FDR<0.05) & (all$logFC) <= -0.58)))
de.log=matrix(c(c1up,c1down,c0.58up,c0.58down),nrow=2)
write.table(de.log,"results/Cutoff.log",sep="\t",row.names=F,col.names=F,quote=F)

library(topGO)

TG = function (ID) {
	geneID2GO <- readMappings(file = "/psc/home/zhaocheng/species/ATH/GO/ATH_tair10_geneid2go.map")
	ALL_gene=rownames(read.table("/psc/home/zhaocheng/species/ATH/Gene.ID",row.names=1))
	geneList <- factor(as.integer(ALL_gene %in% ID))
	names(geneList)=ALL_gene
	GOdata=new("topGOdata",ontology="MF",allGenes=geneList,annot=annFUN.gene2GO,gene2GO=geneID2GO)
	resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
	allRes=GenTable(GOdata,pvalue=resultFisher,topNodes=100)
	allRes$pvalue[grep('<',allRes$pvalue)] = "1e-30"
	allRes$pvalue=as.numeric(allRes$pvalue)
	allRes=allRes[order(as.numeric(allRes$pvalue),decreasing=T),]
	allRes$catagory="MF"
	SM=allRes


	GOdata=new("topGOdata",ontology="BP",allGenes=geneList,annot=annFUN.gene2GO,gene2GO=geneID2GO)
	resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
	allRes=GenTable(GOdata,pvalue=resultFisher,topNodes=100)
	allRes$pvalue[grep('<',allRes$pvalue)] = "1e-30"
	allRes$pvalue=as.numeric(allRes$pvalue)
	allRes=allRes[order(as.numeric(allRes$pvalue),decreasing=T),]
	allRes$catagory="BP"
	BP=allRes


	GOdata=new("topGOdata",ontology="CC",allGenes=geneList,annot=annFUN.gene2GO,gene2GO=geneID2GO)
	resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
	allRes=GenTable(GOdata,pvalue=resultFisher,topNodes=100)
	allRes$pvalue[grep('<',allRes$pvalue)] = "1e-30"
	allRes$pvalue=as.numeric(allRes$pvalue)
	allRes=allRes[order(as.numeric(allRes$pvalue),decreasing=T),]
	allRes$catagory="CC"
	CC=allRes

	OUT=rbind(SM,BP)
	OUT=rbind(OUT,CC)
	OUT
}

up=rownames(subset(all,((all$FDR<0.05) & (all$logFC) >=1)))
down=rownames(subset(all,((all$FDR<0.05) & (all$logFC) <= -1)))
GO_up=TG(up)
GO_down=TG(down)
write.table(GO_up,"results/DE_up.go",sep="\t",col.names=T,row.names=F,quote=F)
write.table(GO_down,"results/DE_down.go",sep="\t",col.names=T,row.names=F,quote=F)
system("perl ~/PC/code/Topgo_buqi.pl /psc/home/zhaocheng/species/GO.terms_and_ids results/DE_up.go > results/DE_up.go.tmp ")
system("mv results/DE_up.go.tmp results/DE_up.go")
system("perl ~/PC/code/Topgo_buqi.pl /psc/home/zhaocheng/species/GO.terms_and_ids results/DE_down.go > results/DE_down.go.tmp ")
system("mv results/DE_down.go.tmp results/DE_down.go")

pdf("results/DE_up_go.pdf")
GR=read.table("results/DE_up.go",header=T,sep="\t",quote="")
GS=subset(GR,GR$pvalue<0.001)
write.table(GS,"results/DE_up.go.sig",sep="\t",col.names=T,row.names=F,quote=F)
PV=as.numeric(GS$pvalue)
names(PV)= GS$Term
func=GS$catagory
GS$catagory=as.vector(GS$catagory)
library(RColorBrewer)
col=brewer.pal(9,"Set1") 
GS$col[GS$catagory=="BP"] <- col[3]
GS$col[GS$catagory=="MF"] <- col[2]
GS$col[GS$catagory=="CC"] <- col[1]
par(mar=c(1,22,5,1))
barplot(-(log10(PV)),cex.names=0.3,width=0.5,space=0.8,las=1,horiz=T,col=GS$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(PV)))+1.5),main=expression(paste("-Log10(",italic(P),"-value of GO ),Up")))
#mtext(expression(paste("-Log10(",italic(P),"-value of GO enrichment)")),at=0,side=3,line=1.5,cex=0.6,adj=0)

legend("right",legend=c("GOCC","GOBP","GOMF"),fill=col[c(1,3,2)],bty="n",cex=0.6,border=F)
axis(3,cex.axis=0.6,mgp=c(1,0.5,0),tck=-0.01)
dev.off()

pdf("results/DE_down_go.pdf")
GR=read.table("results/DE_down.go",header=T,sep="\t",quote="")
GS=subset(GR,GR$pvalue<0.001)
write.table(GS,"results/DE_down.go.sig",sep="\t",col.names=T,row.names=F,quote=F)
PV=as.numeric(GS$pvalue)
names(PV)= GS$Term
func=GS$catagory
GS$catagory=as.vector(GS$catagory)
library(RColorBrewer)
col=brewer.pal(9,"Set1") 
GS$col[GS$catagory=="BP"] <- col[3]
GS$col[GS$catagory=="MF"] <- col[2]
GS$col[GS$catagory=="CC"] <- col[1]
par(mar=c(1,22,5,1))
barplot(-(log10(PV)),cex.names=0.3,width=0.5,space=0.8,las=1,horiz=T,col=GS$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(PV)))+1),main=expression(paste("-Log10(",italic(P),"-value of GO ),Down")))
#mtext(expression(paste("-Log10(",italic(P),"-value of GO enrichment)")),at=0,side=3,line=1.5,cex=0.6,adj=0)

legend("right",legend=c("GOCC","GOBP","GOMF"),fill=col[c(1,3,2)],bty="n",cex=0.6,border=F)
axis(3,cex.axis=0.6,mgp=c(1,0.5,0),tck=-0.01)
dev.off()
