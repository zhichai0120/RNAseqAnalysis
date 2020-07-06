#GO富集分析

source("http://Bioconductor.org/biocLite.R")
#进入bioconductor第三方网站

setwd("c:/R_work/res_R11/");
#设置工作目录


library(clusterProfiler)
#biocLite("org.Hs.eg.db")
#library(org.Hs.eg.db)
##bioconductor website, Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.

#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
#Genome wide annotation for Mouse, primarily based on mapping using Entrez Gene identifiers.

allgene3<-read.csv("20181208allgene3.csv",header=T,row.names=1)
#读取所有基因检验表，赋值给allgene
gene_dif<-read.csv("20181208gene_dif.csv",header=T,row.names=1)
#读取差异基因检验表，赋值给gene_dif
gene_up<-read.csv("20181208gene_up.csv",header=T,row.names=1)
#读取上调基因检验表，赋值给gene_up
gene_down<-read.csv("20181208gene_down.csv",header=T,row.names=1)
#读取下调基因检验表，赋值给gene_down
EGU<-allgene3$entrez
#选取所有基因的EntrezID，赋值给EGU
EG_dif<-gene_dif$entrez
#选取差异基因的EntrezID，赋值给EGS
EG_up<-gene_up$entrez
#选取上调基因gl<-allgene$logFC,赋值给EGS_up
EG_down<-gene_down$entrez
#选取下调基因的EntrezID，赋值给EGS_up

mode(EGU)
mode(EG_dif)

EGU<-as.character(EGU)
#将EGU转化为字符串类型
EG_dif<-as.character(EG_dif)
#将EGU转化为字符串类型
EG_up<-as.character(EG_up)
#将EGU转化为字符串类型
EG_down<-as.character(EG_down)
#将EGU转化为字符串类型
mode(EG_dif)

length(EGU)
length(EG_up)
length(EG_down)
length(EG_dif)

#给allgene3增加一列FoldChange，并给数据赋值为2^log2FoldChange
allgene3$FoldChange<-2^(allgene3$log2FoldChange)

genelist<-allgene3$log2FoldChange
#提取所有基因的Log2FoldChange，赋值给genelist
names(genelist)<-allgene3$entrez
#用EntrezID对变化倍率进行命名
genelist<-sort(genelist,decreasing=T)
#将数据从大到小进行排列

#=======================================up

up_ego_CC<-enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#将差异上调基因在cc细胞组成上进行GO富集分析。
up_ego_MF<- enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#将差异上调基因在MF分子功能上进行GO富集分析。
up_ego_BP<- enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#将差异上调基因在BP生物过程上进行GO富集分析。

write.table(up_ego_CC,file="up_ego_CC.txt",sep="\t",quote=F,row.names=F)
#将差异基因在CC上的GO富集分析结果表进行存储。
write.table(up_ego_MF,file="up_ego_MF.txt",sep="\t",quote=F,row.names=F)
#将差异基因MF上的GO富集分析结果表进行存储。
write.table(up_ego_BP,file="up_ego_BP.txt",sep="\t",quote=F,row.names=F)
#将差异基因在BP上的GO富集分析结果表进行存储。



tiff(file="up_BP_barplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
barplot(up_ego_BP, showCategory=10)
#GO富集条形图
dev.off()
#关闭设备，保存图片。

tiff(file="up_BP_dotplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
dotplot(up_ego_BP,showCategory=10)
#GO富集散点图
dev.off()
#关闭设备，保存图片。

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ‘pvalue’, ‘p.adjust’ or ‘qvalue’.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"


tiff(file="up_BP_enrichMap.tif",res=100,units='in',width=20,height=20)
#打开tiff绘图设备，命名图片。
enrichMap(up_ego_BP)
#富集的GO term之间关系图
dev.off()
#关闭设备，保存图片。



tiff(file="up_BP_cnetplot2.tif",res=150,units='in',width=20,height=15)
#打开tiff绘图设备，命名图片。
cnetplot(up_ego_BP,showCategory =5,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1,vertext.label.font=20)
#GO与基因关系的网络图。
dev.off()
#关闭设备，保存图片。
## categorySize can be scaled by 'pvalue' or 'geneNum'


biocLite("topGO")
library(topGO)
tiff(file="dif_CC_topGO.tif",res=400,units='in',width=10,height=10)
#打开tiff绘图设备
plotGOgraph(dif_ego_CC,firstSigNodes =5)
#topGO有向无环图DAG。
dev.off()
#关闭设备，保存图片。

#============================down

down_ego_CC<-enrichGO(gene=EG_down,universe=EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#将差异下调基因在cc细胞组成上进行GO富集分析。
down_ego_MF<-enrichGO(gene=EG_down,universe=EGU,OrgDb= org.Mm.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#将差异下调基因在MF分子功能上进行GO富集分析。
down_ego_BP<-enrichGO(gene=EG_down,universe=EGU,OrgDb= org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#将差异下调基因在BP生物过程上进行GO富集分析。

write.table(down_ego_CC,file="down_ego_CC.txt",sep="\t",quote=F,row.names=F)
#将差异基因在CC上的GO富集分析结果表进行存储。
write.table(down_ego_MF,file="down_ego_MF.txt",sep="\t",quote=F,row.names=F)
#将差异基因MF上的GO富集分析结果表进行存储。
write.table(down_ego_BP,file="down_ego_BP.txt",sep="\t",quote=F,row.names=F)
#将差异基因在BP上的GO富集分析结果表进行存储。

tiff(file="down_MF_barplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
barplot(down_ego_MF, showCategory=10)
#GO富集条形图
dev.off()
#关闭设备，保存图片。

tiff(file="down_MF_dotplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
dotplot(down_ego_MF,showCategory=10)
#GO富集散点图
dev.off()
#关闭设备，保存图片。

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ‘pvalue’, ‘p.adjust’ or ‘qvalue’.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"


tiff(file="down_MF_enrichMap.tif",res=100,units='in',width=20,height=20)
#打开tiff绘图设备，命名图片。
enrichMap(down_ego_MF)
#富集的GO term之间关系图
dev.off()
#关闭设备，保存图片。



tiff(file="down_MF_cnetplot.tif",res=150,units='in',width=20,height=15)
#打开tiff绘图设备，命名图片。
cnetplot(down_ego_MF,showCategory =5,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#GO与基因关系的网络图。
dev.off()
#关闭设备，保存图片。
## categorySize can be scaled by 'pvalue' or 'geneNum'


biocLite("topGO")
library(topGO)
tiff(file="down_MF_topGO.tif",res=400,units='in',width=10,height=10)
#打开tiff绘图设备
plotGOgraph(down_ego_MF,firstSigNodes =5)
#topGO有向无环图DAG。
dev.off()
#关闭设备，保存图片。


#====================both up and down====
dif_ego_CC<-enrichGO(gene=EG_dif,universe= EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
#将差异基因在cc细胞组成上进行GO富集分析。
dif_ego_MF<- enrichGO(gene=EG_dif,universe= EGU,OrgDb= org.Mm.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
#将差异基因在MF分子功能上进行GO富集分析。
dif_ego_BP<- enrichGO(gene=EG_dif,universe= EGU,OrgDb= org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
#将差异基因在BP生物过程上进行GO富集分析。

write.table(dif_ego_CC,file="dif_ego_CC.txt",sep="\t",quote=F,row.names=F)
#将差异基因在CC上的GO富集分析结果表进行存储。
write.table(dif_ego_MF,file="dif_ego_MF.txt",sep="\t",quote=F,row.names=F)
#将差异基因MF上的GO富集分析结果表进行存储。
write.table(dif_ego_BP,file="dif_ego_BP.txt",sep="\t",quote=F,row.names=F)
#将差异基因在BP上的GO富集分析结果表进行存储。


tiff(file="dif_MF_barplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
barplot(dif_ego_MF, showCategory=10)
#font.size=26
#GO富集条形图
dev.off()
#关闭设备，保存图片。

tiff(file="dif_MF_dotplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
dotplot(dif_ego_MF,showCategory=10,font.size =32)
#GO富集散点图
dev.off()
#关闭设备，保存图片。

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ‘pvalue’, ‘p.adjust’ or ‘qvalue’.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"




tiff(file="dif_MF_enrichMap.tif",res=100,units='in',width=20,height=20)
#打开tiff绘图设备，命名图片。
enrichMap(dif_ego_MF)
#富集的GO term之间关系图
dev.off()
#关闭设备，保存图片。



tiff(file="dif_MF_cnetplot2.tif",res=150,units='in',width=20,height=15)
#打开tiff绘图设备，命名图片。
cnetplot(dif_ego_MF,showCategory=25 ,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#GO与基因关系的网络图。
dev.off()
#关闭设备，保存图片。
## categorySize can be scaled by 'pvalue' or 'geneNum'


biocLite("topGO")
library(topGO)
tiff(file="dif_CC_topGO.tif",res=400,units='in',width=10,height=10)
#打开tiff绘图设备
plotGOgraph(dif_ego_CC,firstSigNodes =5)
#topGO有向无环图DAG。
dev.off()
#关闭设备，保存图片。


#=====dif, BP==
  
tiff(file="dif_BP_barplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
barplot(dif_ego_BP, showCategory=10)
#GO富集条形图
dev.off()
#关闭设备，保存图片。

tiff(file="dif_BP_dotplot.tif",res=100,units='in',width=10,height=10)
#打开tiff绘图设备，命名图片。
dotplot(dif_ego_BP,showCategory=10)
#GO富集散点图
dev.off()
#关闭设备，保存图片。

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ‘pvalue’, ‘p.adjust’ or ‘qvalue’.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"




tiff(file="dif_BP_enrichMap.tif",res=100,units='in',width=20,height=20)
#打开tiff绘图设备，命名图片。
enrichMap(dif_ego_BP)
#富集的GO term之间关系图
dev.off()
#关闭设备，保存图片。



tiff(file="dif_BP_cnetplot2.tif",res=150,units='in',width=20,height=15)
#打开tiff绘图设备，命名图片。
cnetplot(dif_ego_BP,showCategory =15,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#GO与基因关系的网络图。
dev.off()
#关闭设备，保存图片。
## categorySize can be scaled by 'pvalue' or 'geneNum'

biocLite("topGO")
library(topGO)
tiff(file="dif_BP_topGO.tif",res=400,units='in',width=10,height=10)
#打开tiff绘图设备
plotGOgraph(dif_ego_BP,firstSigNodes =5)
#topGO有向无环图DAG。
dev.off()
#关闭设备，保存图片






