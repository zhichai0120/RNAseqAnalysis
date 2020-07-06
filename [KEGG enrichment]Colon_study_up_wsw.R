source("http://Bioconductor.org/biocLite.R")
#����bioconductor��������վ

setwd("c:/R_work/deseq2/")
#���ù���Ŀ¼

#biocLite("clusterProfiler")
#install.packages
#biocLite("stringi")
library(clusterProfiler)

allgene4<-read.csv("20181208allgene3.csv",header=T,row.names=1)
#allgene<-read.csv("allgene2.csv",header=T,row.names=1)
#��ȡ���л�����������ֵ��allgene
gene_up<-read.csv("20181208gene_up.csv",header=T,row.names=1)
#��ȡ���������������ֵ��gene_dif


EGU<-allgene4$entrez
#ѡȡ���л����EntrezID����ֵ��EGU
mode(EGU)
EGU<-as.character(EGU)
#��EGUת��Ϊ�ַ�������
mode(EGU)
EG_up<-gene_up$entrez
#ѡȡ��������EntrezID����ֵ��EG_dif
EG_up<-as.character(EG_up)
#��EG_difת��Ϊ�ַ�������

genelist<-allgene4$log2FoldChange
#��ȡ���л���ı仯����ֵ����ֵ��genelist
names(genelist)<-allgene4$entrez
#��EntrezID�Ա仯���ʽ�������
genelist<-sort(genelist,decreasing=T)
#�����ݴӴ�С��������

up_genelist<-gene_up$log2FoldChange
#��ȡ�������ı仯����ֵ����ֵ��dif_genelist
names(up_genelist)<-gene_up$entrez
#��EntrezID�Բ������ı仯���ʽ�������

#up_kk <- enrichKEGG(gene=EG_up,universe= EGU,organism= 'mmu',pvalueCutoff = 0.05)
up_kk <- enrichKEGG(gene=EG_up,organism= 'mmu',pvalueCutoff = 0.05)
#up_ego_CC<-enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#dif_kk <- enrichKEGG(gene=EG_dif,organism= 'hsa',pvalueCutoff = 0.05)
#other organisms ?enrichKEGG   --> supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'

#�����������KEGG���и���������
write.table(up_kk,file="up_kk.txt",sep="\t",quote=F,row.names=F)
#�������������KEGG���������Ľ�������д洢��
barplot(up_kk,showCategory=10)
#KEGG��������ͼ
dotplot(up_kk,showCategory=10)
#KEGG����ɢ��ͼ
enrichMap(up_kk)
#����ͨ·֮��Ĺ�ϵͼ
cnetplot(up_kk,showCategory =10,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#����ͨ·�����֮��Ĺ�ϵͼ



browseKEGG(up_kk, 'mmu05150')
#��ͨ·ͼ����ҳ��




#biocLite("stringi")
#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

#biocLite("pathview")
library(pathview)
mmu05150 <- pathview(gene.data=up_genelist,pathway.id = "mmu05150",species= "mmu",limit=list(gene=max(abs(genelist)), cpd=1))
#�����ӻ���ͨ·ͼ