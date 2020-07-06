#GO��������

source("http://Bioconductor.org/biocLite.R")
#����bioconductor��������վ

setwd("c:/R_work/res_R11/");
#���ù���Ŀ¼


library(clusterProfiler)
#biocLite("org.Hs.eg.db")
#library(org.Hs.eg.db)
##bioconductor website, Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.

#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
#Genome wide annotation for Mouse, primarily based on mapping using Entrez Gene identifiers.

allgene3<-read.csv("20181208allgene3.csv",header=T,row.names=1)
#��ȡ���л�����������ֵ��allgene
gene_dif<-read.csv("20181208gene_dif.csv",header=T,row.names=1)
#��ȡ���������������ֵ��gene_dif
gene_up<-read.csv("20181208gene_up.csv",header=T,row.names=1)
#��ȡ�ϵ�������������ֵ��gene_up
gene_down<-read.csv("20181208gene_down.csv",header=T,row.names=1)
#��ȡ�µ�������������ֵ��gene_down
EGU<-allgene3$entrez
#ѡȡ���л����EntrezID����ֵ��EGU
EG_dif<-gene_dif$entrez
#ѡȡ��������EntrezID����ֵ��EGS
EG_up<-gene_up$entrez
#ѡȡ�ϵ�����gl<-allgene$logFC,��ֵ��EGS_up
EG_down<-gene_down$entrez
#ѡȡ�µ������EntrezID����ֵ��EGS_up

mode(EGU)
mode(EG_dif)

EGU<-as.character(EGU)
#��EGUת��Ϊ�ַ�������
EG_dif<-as.character(EG_dif)
#��EGUת��Ϊ�ַ�������
EG_up<-as.character(EG_up)
#��EGUת��Ϊ�ַ�������
EG_down<-as.character(EG_down)
#��EGUת��Ϊ�ַ�������
mode(EG_dif)

length(EGU)
length(EG_up)
length(EG_down)
length(EG_dif)

#��allgene3����һ��FoldChange���������ݸ�ֵΪ2^log2FoldChange
allgene3$FoldChange<-2^(allgene3$log2FoldChange)

genelist<-allgene3$log2FoldChange
#��ȡ���л����Log2FoldChange����ֵ��genelist
names(genelist)<-allgene3$entrez
#��EntrezID�Ա仯���ʽ�������
genelist<-sort(genelist,decreasing=T)
#�����ݴӴ�С��������

#=======================================up

up_ego_CC<-enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#�������ϵ�������ccϸ������Ͻ���GO����������
up_ego_MF<- enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#�������ϵ�������MF���ӹ����Ͻ���GO����������
up_ego_BP<- enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#�������ϵ�������BP��������Ͻ���GO����������

write.table(up_ego_CC,file="up_ego_CC.txt",sep="\t",quote=F,row.names=F)
#�����������CC�ϵ�GO����������������д洢��
write.table(up_ego_MF,file="up_ego_MF.txt",sep="\t",quote=F,row.names=F)
#���������MF�ϵ�GO����������������д洢��
write.table(up_ego_BP,file="up_ego_BP.txt",sep="\t",quote=F,row.names=F)
#�����������BP�ϵ�GO����������������д洢��



tiff(file="up_BP_barplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
barplot(up_ego_BP, showCategory=10)
#GO��������ͼ
dev.off()
#�ر��豸������ͼƬ��

tiff(file="up_BP_dotplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
dotplot(up_ego_BP,showCategory=10)
#GO����ɢ��ͼ
dev.off()
#�ر��豸������ͼƬ��

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ��pvalue��, ��p.adjust�� or ��qvalue��.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"


tiff(file="up_BP_enrichMap.tif",res=100,units='in',width=20,height=20)
#��tiff��ͼ�豸������ͼƬ��
enrichMap(up_ego_BP)
#������GO term֮���ϵͼ
dev.off()
#�ر��豸������ͼƬ��



tiff(file="up_BP_cnetplot2.tif",res=150,units='in',width=20,height=15)
#��tiff��ͼ�豸������ͼƬ��
cnetplot(up_ego_BP,showCategory =5,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1,vertext.label.font=20)
#GO������ϵ������ͼ��
dev.off()
#�ر��豸������ͼƬ��
## categorySize can be scaled by 'pvalue' or 'geneNum'


biocLite("topGO")
library(topGO)
tiff(file="dif_CC_topGO.tif",res=400,units='in',width=10,height=10)
#��tiff��ͼ�豸
plotGOgraph(dif_ego_CC,firstSigNodes =5)
#topGO�����޻�ͼDAG��
dev.off()
#�ر��豸������ͼƬ��

#============================down

down_ego_CC<-enrichGO(gene=EG_down,universe=EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#�������µ�������ccϸ������Ͻ���GO����������
down_ego_MF<-enrichGO(gene=EG_down,universe=EGU,OrgDb= org.Mm.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#�������µ�������MF���ӹ����Ͻ���GO����������
down_ego_BP<-enrichGO(gene=EG_down,universe=EGU,OrgDb= org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#�������µ�������BP��������Ͻ���GO����������

write.table(down_ego_CC,file="down_ego_CC.txt",sep="\t",quote=F,row.names=F)
#�����������CC�ϵ�GO����������������д洢��
write.table(down_ego_MF,file="down_ego_MF.txt",sep="\t",quote=F,row.names=F)
#���������MF�ϵ�GO����������������д洢��
write.table(down_ego_BP,file="down_ego_BP.txt",sep="\t",quote=F,row.names=F)
#�����������BP�ϵ�GO����������������д洢��

tiff(file="down_MF_barplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
barplot(down_ego_MF, showCategory=10)
#GO��������ͼ
dev.off()
#�ر��豸������ͼƬ��

tiff(file="down_MF_dotplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
dotplot(down_ego_MF,showCategory=10)
#GO����ɢ��ͼ
dev.off()
#�ر��豸������ͼƬ��

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ��pvalue��, ��p.adjust�� or ��qvalue��.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"


tiff(file="down_MF_enrichMap.tif",res=100,units='in',width=20,height=20)
#��tiff��ͼ�豸������ͼƬ��
enrichMap(down_ego_MF)
#������GO term֮���ϵͼ
dev.off()
#�ر��豸������ͼƬ��



tiff(file="down_MF_cnetplot.tif",res=150,units='in',width=20,height=15)
#��tiff��ͼ�豸������ͼƬ��
cnetplot(down_ego_MF,showCategory =5,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#GO������ϵ������ͼ��
dev.off()
#�ر��豸������ͼƬ��
## categorySize can be scaled by 'pvalue' or 'geneNum'


biocLite("topGO")
library(topGO)
tiff(file="down_MF_topGO.tif",res=400,units='in',width=10,height=10)
#��tiff��ͼ�豸
plotGOgraph(down_ego_MF,firstSigNodes =5)
#topGO�����޻�ͼDAG��
dev.off()
#�ر��豸������ͼƬ��


#====================both up and down====
dif_ego_CC<-enrichGO(gene=EG_dif,universe= EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
#�����������ccϸ������Ͻ���GO����������
dif_ego_MF<- enrichGO(gene=EG_dif,universe= EGU,OrgDb= org.Mm.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
#�����������MF���ӹ����Ͻ���GO����������
dif_ego_BP<- enrichGO(gene=EG_dif,universe= EGU,OrgDb= org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
#�����������BP��������Ͻ���GO����������

write.table(dif_ego_CC,file="dif_ego_CC.txt",sep="\t",quote=F,row.names=F)
#�����������CC�ϵ�GO����������������д洢��
write.table(dif_ego_MF,file="dif_ego_MF.txt",sep="\t",quote=F,row.names=F)
#���������MF�ϵ�GO����������������д洢��
write.table(dif_ego_BP,file="dif_ego_BP.txt",sep="\t",quote=F,row.names=F)
#�����������BP�ϵ�GO����������������д洢��


tiff(file="dif_MF_barplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
barplot(dif_ego_MF, showCategory=10)
#font.size=26
#GO��������ͼ
dev.off()
#�ر��豸������ͼƬ��

tiff(file="dif_MF_dotplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
dotplot(dif_ego_MF,showCategory=10,font.size =32)
#GO����ɢ��ͼ
dev.off()
#�ر��豸������ͼƬ��

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ��pvalue��, ��p.adjust�� or ��qvalue��.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"




tiff(file="dif_MF_enrichMap.tif",res=100,units='in',width=20,height=20)
#��tiff��ͼ�豸������ͼƬ��
enrichMap(dif_ego_MF)
#������GO term֮���ϵͼ
dev.off()
#�ر��豸������ͼƬ��



tiff(file="dif_MF_cnetplot2.tif",res=150,units='in',width=20,height=15)
#��tiff��ͼ�豸������ͼƬ��
cnetplot(dif_ego_MF,showCategory=25 ,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#GO������ϵ������ͼ��
dev.off()
#�ر��豸������ͼƬ��
## categorySize can be scaled by 'pvalue' or 'geneNum'


biocLite("topGO")
library(topGO)
tiff(file="dif_CC_topGO.tif",res=400,units='in',width=10,height=10)
#��tiff��ͼ�豸
plotGOgraph(dif_ego_CC,firstSigNodes =5)
#topGO�����޻�ͼDAG��
dev.off()
#�ر��豸������ͼƬ��


#=====dif, BP==
  
tiff(file="dif_BP_barplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
barplot(dif_ego_BP, showCategory=10)
#GO��������ͼ
dev.off()
#�ر��豸������ͼƬ��

tiff(file="dif_BP_dotplot.tif",res=100,units='in',width=10,height=10)
#��tiff��ͼ�豸������ͼƬ��
dotplot(dif_ego_BP,showCategory=10)
#GO����ɢ��ͼ
dev.off()
#�ر��豸������ͼƬ��

#https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
dotplot(do, x="count", showCategory=20, colorBy="qvalue")
#We can set the x-axis to use gene count and dot color by one of ��pvalue��, ��p.adjust�� or ��qvalue��.
#x- "GeneRatio"(default) or "Count"
#color- "p/adjust" or "pvalue" or "qvalue"




tiff(file="dif_BP_enrichMap.tif",res=100,units='in',width=20,height=20)
#��tiff��ͼ�豸������ͼƬ��
enrichMap(dif_ego_BP)
#������GO term֮���ϵͼ
dev.off()
#�ر��豸������ͼƬ��



tiff(file="dif_BP_cnetplot2.tif",res=150,units='in',width=20,height=15)
#��tiff��ͼ�豸������ͼƬ��
cnetplot(dif_ego_BP,showCategory =15,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#GO������ϵ������ͼ��
dev.off()
#�ر��豸������ͼƬ��
## categorySize can be scaled by 'pvalue' or 'geneNum'

biocLite("topGO")
library(topGO)
tiff(file="dif_BP_topGO.tif",res=400,units='in',width=10,height=10)
#��tiff��ͼ�豸
plotGOgraph(dif_ego_BP,firstSigNodes =5)
#topGO�����޻�ͼDAG��
dev.off()
#�ر��豸������ͼƬ





