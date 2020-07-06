source("http://Bioconductor.org/biocLite.R")
#进入bioconductor第三方网站

setwd("c:/R_work/deseq2/")
#设置工作目录

#biocLite("clusterProfiler")
#install.packages
#biocLite("stringi")
library(clusterProfiler)

allgene4<-read.csv("20181208allgene3.csv",header=T,row.names=1)
#allgene<-read.csv("allgene2.csv",header=T,row.names=1)
#读取所有基因检验表，赋值给allgene
gene_up<-read.csv("20181208gene_up.csv",header=T,row.names=1)
#读取差异基因检验表，赋值给gene_dif


EGU<-allgene4$entrez
#选取所有基因的EntrezID，赋值给EGU
mode(EGU)
EGU<-as.character(EGU)
#将EGU转化为字符串类型
mode(EGU)
EG_up<-gene_up$entrez
#选取差异基因的EntrezID，赋值给EG_dif
EG_up<-as.character(EG_up)
#将EG_dif转化为字符串类型

genelist<-allgene4$log2FoldChange
#提取所有基因的变化倍率值，赋值给genelist
names(genelist)<-allgene4$entrez
#用EntrezID对变化倍率进行命名
genelist<-sort(genelist,decreasing=T)
#将数据从大到小进行排列

up_genelist<-gene_up$log2FoldChange
#提取差异基因的变化倍率值，赋值给dif_genelist
names(up_genelist)<-gene_up$entrez
#用EntrezID对差异基因的变化倍率进行命名

#up_kk <- enrichKEGG(gene=EG_up,universe= EGU,organism= 'mmu',pvalueCutoff = 0.05)
up_kk <- enrichKEGG(gene=EG_up,organism= 'mmu',pvalueCutoff = 0.05)
#up_ego_CC<-enrichGO(gene=EG_up,universe= EGU,OrgDb= org.Mm.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)
#dif_kk <- enrichKEGG(gene=EG_dif,organism= 'hsa',pvalueCutoff = 0.05)
#other organisms ?enrichKEGG   --> supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'

#将差异基因在KEGG进行富集分析。
write.table(up_kk,file="up_kk.txt",sep="\t",quote=F,row.names=F)
#将差异基因在在KEGG富集分析的结果表进行存储。
barplot(up_kk,showCategory=10)
#KEGG富集条形图
dotplot(up_kk,showCategory=10)
#KEGG富集散点图
enrichMap(up_kk)
#富集通路之间的关系图
cnetplot(up_kk,showCategory =10,categorySize="pvalue", foldChange=genelist,vertex.label.cex=1)
#富集通路与基因之间的关系图



browseKEGG(up_kk, 'mmu05150')
#打开通路图的网页。




#biocLite("stringi")
#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

#biocLite("pathview")
library(pathview)
mmu05150 <- pathview(gene.data=up_genelist,pathway.id = "mmu05150",species= "mmu",limit=list(gene=max(abs(genelist)), cpd=1))
#做可视化的通路图