setwd("c:/R_work/deseq2/");
rm(list=ls(all=TRUE));
##install DESeq2
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")


library("DESeq2")

# loading data
dat0 <- read.table(file="data_chai/counts_R04.txt", header = T, sep = "\t")
dat1 <- dat0[,2:13]
#dat1 <- dat1[,-c(4,5,6,10,11,12)]

row.names(dat1) <- dat0[,1]

# removing low counts
hist(log(rowSums(dat1)),breaks=50)
flt1 <- which(rowSums(dat1<10)>8)
flt2 <- which(rowSums(dat1)<200)
flt <- unique(c(flt1,flt2)); length(flt)
dat1 <- dat1[-flt,]
hist(log(rowSums(dat1)),breaks=50)


# cluster & EDA
dat1log <- log(dat1+1)
d <- dist(t(dat1log),method = "euclidean")
hc <-hclust(d, method = "complete", members = NULL)
plot(hc)

lattice::levelplot(cor(dat1log))

# cluster & EDA using rank
dat1.rk <- apply(-dat1,2,rank,ties.method="random")
d <- dist(t(dat1.rk),method = "euclidean")
hc <-hclust(d, method = "complete", members = NULL)
plot(hc)

lattice::levelplot(cor(dat1log,method = "spearman"))





### construct dds

desi <- cbind(c(rep("VAD",6),rep("VAS",6)),
              rep(c(rep("N_inf",3),rep("inf",3)),2),
              rep(c("r1","r2","r3"),4))



row.names(desi) <- colnames(dat1)
colnames(desi) <- c("VA","is.inf","reps")
desi <- as.data.frame(desi)





## with interaction

dds2 <- DESeqDataSetFromMatrix(countData = dat1,colData = desi,design = ~VA+is.inf+VA:is.inf)

#set "VAD" as the reference level of vitamin A effect (vA)
#dds2$VA <- relevel(dds2$VA,ref="VAD")
#set "N_inf" as the reference level of infection effect (is.inf)
dds2$is.inf <- relevel(dds2$is.inf,ref="N_inf")
##check if the first one is the reference level
levels(dds2$is.inf)
levels(dds2$VA)

#####->###
dds2 <- DESeq(dds2)
resultsNames(dds2)



#the Infection effect for VAD (the main effect)(VAD is the reference level of VA effect, genotype I) 
res_infectioneffectforVAD <- results(dds2,cooksCutoff=FALSE, contrast=c("is.inf","inf","N_inf"),alpha=0.05)
#when cooksCutoff=FALSE, the automatic filtering for outlier is turned off, so there won't be outliers (p=adjp="NA")
summary(res_infectioneffectforVAD)
##infectioneffectforVAD is a DataFrame object
mcols(res_infectioneffectforVAD, use.names=TRUE)

sig.infectioneffectforVAD <- res_infectioneffectforVAD$padj<0.05
sum(sig.infectioneffectforVAD,na.rm = T)


####print normalized count, relatively high expression count, for all genes
write.csv(counts(dds2, normalized=TRUE), file = "go_enrichment/20180914_R04LSI_highcount_normalizedbyDESeq2.csv",quote = F)
## print an excel of all p and adj P values of all genes
write.csv(res_infectioneffectforVAD, file = "go_enrichment/20180914_R04LSI_results_infectioneffect_forVAD.csv",quote = F, row.names = rownames(dat1))


##MA plot, colored red means adj p less than 0.1
plotMA(res_infectioneffectforVAD,ylim=c(-2,2))






##the Infection effect for VAS
##this is, by definition, the main effect *plus* the interaction term (VAS is not the reference level of VA effect, VAS~genotype II) 
##interaction term: (the extra infection effect in VAS compared to VAD)
res_infectioneffectforVAS <- results(dds2,cooksCutoff=FALSE,list(c("is.inf_inf_vs_N_inf","VAVAS.is.infinf")))
summary(res_infectioneffectforVAS)
##infectioneffectforVAS is a DataFrame object
mcols(res_infectioneffectforVAS, use.names=TRUE)

sig.infectioneffectforVAS <- res_infectioneffectforVAS$padj<0.05
sum(sig.infectioneffectforVAS,na.rm = T)


#rownames(dat1)[which(sig.infectioneffectforVAS)]
#write.table(x = rownames(dat1)[which(sig.infectioneffectforVAS)], file = "go_enrichment/colon_infectioneffect_forVAS.txt",quote = F, row.names = F, col.names = F)

## print an excel of all count of significant genes
#write.table(counts(dds2[which(sig.infectioneffectforVAS)], normalized=TRUE), file = "go_enrichment/colon_expr_infectioneffect_forVAS_normalized.txt",quote = F, row.names = rownames(dat1)[which(sig.infectioneffectforVAS)])
# print an excel of all p and adj P values of significant genes
#write.table(res_infectioneffectforVAS[which(sig.infectioneffectforVAS),], file = "go_enrichment/colon_expr2_infectioneffect_forVAS.txt",quote = F, row.names = rownames(dat1)[which(sig.infectioneffectforVAS)])

####print normalized count, relatively high expression count, for all genes
#write.csv(counts(dds2, normalized=TRUE), file = "go_enrichment/20180914_R04LSI_highcount_normalizedbyDESeq2.csv",quote = F)
## print an excel of all p and adj P values of all genes
write.csv(res_infectioneffectforVAS, file = "go_enrichment/20180914_R04LSI_results_infectioneffect_forVAS.csv",quote = F, row.names = rownames(dat1))






##the interaction term, answering is the infection effect *different* across VA status?
res_interterm <- results(dds2,cooksCutoff=FALSE,names="VAVAS.is.infinf")
summary(res_interterm) 
##res_interterm is a DataFrame object
mcols(res_interterm, use.names=TRUE)

sig.interterm <- res_interterm$padj<0.05
sum(sig.interterm,na.rm = T)

#rownames(dat1)[which(sig.interterm)]
#write.table(x = rownames(dat1)[which(sig.interterm)], file = "go_enrichment/LSI_interterm.txt",quote = F, row.names = F, col.names = F)

## print an excel of all count of significant genes
#write.table(counts(dds2[which(sig.interterm)], normalized=TRUE), file = "go_enrichment/colon_expr_interterm_normalized.txt",quote = F, row.names = rownames(dat1)[which(sig.interterm)])
# print an excel of all p and adj P values of significant genes
#write.table(res_interterm[which(sig.interterm),], file = "go_enrichment/colon_expr2_interterm.txt",quote = F, row.names = rownames(dat1)[which(sig.interterm)])

## print an excel of all p and adj P values of all genes
write.csv(res_interterm, file = "go_enrichment/20180914_R04LSI_results_Interactioneffect.csv",quote = F, row.names = rownames(dat1))






###In general, the results for a comparison of any two levels of a variable can be extracted using the contrast argument to results. The user should specify three values: the name of the variable, 
###the name of the level for the numerator, and the name of the level for the denominator. Here we extract results for the log2 of the fold change of one cell line over another:
###  results(dds, contrast = c("cell", "N061011", "N61311"))
res_VA <- results(dds2, cooksCutoff=FALSE,name="VA_VAS_vs_VAD")
#check if there's any outliers
summary(res_VA)

##res_VA is a DataFrame object
mcols(res_VA, use.names=TRUE)

##subset
#res_VA_Sig <- subset(res_VA,padj<0.05)
#dim(res_VA_Sig)

#res_VA_Sig_big <- subset(res_VA,log2FoldChange>1)
#dim(res_VA_Sig_big)

sig.VA <- res_VA$padj<0.05
sum(sig.VA,na.rm = T)


#rownames(dat1)[which(sig.VA)]
#write.table(x = rownames(dat1)[which(sig.VA)], file = "go_enrichment/colon_gene_VA_woinfection LSI.txt",quote = F, row.names = F, col.names = F)

####print normalized count, relatively high expression count, for all genes
write.csv(counts(dds2, normalized=TRUE), file = "go_enrichment/20180914_R04LSI_highcount_normalizedbyDESeq2.csv",quote = F)
## print an excel of all p and adj P values of all genes
write.csv(res_VA, file = "go_enrichment/20180914_LSI_results_VAeffectWithoutinfection.csv",quote = F, row.names = rownames(dat1))






#volcano plot [实验万事屋-R语言轻松入门-筛选差异基因（下）]
#for VA effect without infection
allgene<- res_VA

#给allgene增加一列type，并给数据赋值为NA
allgene$type<-NA

for(i in 1:nrow(allgene))if(allgene[i,2]>1&allgene[i,6]<0.05){allgene[i,7]<-"up"}else if(-allgene[i,2]>1&allgene[i,6]<0.05){allgene[i,7]<-"down"}else{allgene[i,7]<-"normal"}
#通过循环语句，将基因类型分为三类（log2FoldChange>1，padj<0.05的为上调基因up；log2FoldChange<-1,adj.padj<0.05的为下调基因down，其他为normal）
sum(allgene$type== "up")
sum(allgene$type== "down")

write.csv(allgene,file="20180930allgene.csv",quote=F)
#将新的所有基因的检测报告，输出为csv格式。

gene_up<-allgene[allgene$type=="up",]
#筛选上调差异基因，并存储于变量gene_up中

gene_down<-allgene[allgene$type=="down",]
#筛选下调差异基因，并存储于变量gene_down中

gene_dif<-allgene[allgene$type!="normal",]
#筛选非差异基因，并存储于变量gene_dif中

write.csv(gene_up,file="20180930gene_up.csv",quote=F)
#将筛选的上调基因检测报告，保存为csv格式
write.csv(gene_down,file="20180930gene_down.csv",quote=F)
#将筛选的下调基因检测报告，保存为csv格式
write.csv(gene_dif,file="20180930gene_dif.csv",quote=F)
#将筛选的差异基因检测报告，保存为csv格式

######heatmap of all DE genes (only need DE gene names & expression levels
# don't cluster on column)
#source("http://Bioconductor.org/biocLite.R")
#biocLite("pheatmap")
#下载pheatmap包，用于画热图
####install the heatmap tool

library(pheatmap)
#加载pheatmap包

gene_dif_exp<-counts(dds2, normalized=TRUE)[rownames(gene_dif),]
#提取差异基因的基因表达值

#rownames(gene_dif_exp)<-gene_dif$symbols
#将差异基因表达矩阵的行名转化为gene symbol ID名

write.csv(gene_dif_exp,file="gene_dif_exp.csv",quote=F)
#将筛选的差异基因表达，保存为csv格式

tiff(file="gene_dif_pheatmap.tif",res=300,units='in',width=30,height=30)
#打开tiff绘图设备

pheatmap(gene_dif_exp,cluster_cols = FALSE,cluster_rows = TRUE,
         color=colorRampPalette(c("blue","white","red"))(100),fontsize_row=32,scale="row",border_color=NA)
#画差异基因热图

dev.off()
#保存，关闭绘图设备



###draw a heatmap for the interaction effect (already normalized by DESeq2)

# loading data
dat10 <- read.table(file="data_chai/interaction FC high.txt", header = T, sep = "\t")
dat11 <- dat10[,2:17]
row.names(dat11) <- dat10[,1]

####install the heatmap tool
#source("http://Bioconductor.org/biocLite.R")
#biocLite("pheatmap")

library(pheatmap)
#select genes in the color/module
gene_dif_exp <- dat11
#write a table
write.csv(gene_dif_exp,file="gene_dif_exp.csv",quote=F)
#turn on a graph tool called "tiff"
tiff(file="gene_dif_pheatmap.tif",res=300,units='in',width=30,height=30)
#draw the heatmap, the file will be saved, not shown in R because of the resolution
pheatmap(gene_dif_exp,color=colorRampPalette(c("green","black","red"))(100),fontsize_row=4,scale="row",border_color=NA)
#save and turn off the tiff tool
dev.off()









##Volcano plot

#install.packages("ggplot2")
##下载ggplot2包

library(ggplot2)
#加载ggplot2包

threshold<-factor(allgene$type)
#将type转化为因子向量，并存储于变量threshold中，用来区分火山图差异基因颜色鑹?
#here allgene does not include those whose EntrezID is missing


# load data
allgene2 <- read.csv(file="20180930allgene.csv", header=T,row.names=1)



# Install ggrepel package if needed
#install.packages("devtools")
#devtools::install_github("slowkow/ggrepel")
library("ggrepel") #Avoid overlapping labels

#install.packages("dplyr")
library(dplyr)


tiff(file="20180930volcano.tif",res=100,units='in',width=10,height=10)
#打开tiff格式的绘图设备，并命名为volcano;res为像素，unit为单位；details see ?tiff

ggplot(results2,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+
xlab("log2FC")+ylab("-log10adj.P")+
  ggtitle("Volcano Plot")+theme(plot.title=element_text(hjust=0.5))+
  geom_point()+geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+
  theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "black", "red"))+
  expand_limits(x=c(-8,8))
##original code from wanshiwu
#用ggplot??火山图,Criteria for DE, |Foldchange|=2, so draw lines x=Log2FoldChange=-1 and 1; Criteria for DE,adj p<0.05, so draw line y= -log10(0.05)
#aes对坐标进行设置；”+“ layer
dev.off()
#close drawing device 保存关闭绘图设备


#ggplot(allgene2,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+
#  xlab("log2FC")+ylab("-log10adj.P")+ggtitle("Volcano Plot")+
#  theme(plot.title=element_text(hjust=0.5))+geom_point()+
#  geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+
#  theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "black", "red"))+
#  expand_limits(x=c(-8,8))+
#  geom_text(aes(label=rownames(allgene2)))
#prototype- with text label,  but a lot of overlaps, all genes have labels

#ggplot(results,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+
#  xlab("log2FC")+ylab("-log10adj.P")+ggtitle("Volcano Plot")+
#  theme(plot.title=element_text(hjust=0.5))+geom_point()+
#  geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+
  theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "transparent", "red"))+
  expand_limits(x=c(-8,8))+
  geom_text_repel(aes(label=rownames(allgene2)))
#with text label & label won't interfere with dots, hide un-DE genes


#p=ggplot(allgene2,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+
#  xlab("log2FC")+ylab("-log10adj.P")+ggtitle("Volcano Plot")+
#  theme(plot.title=element_text(hjust=0.5))+geom_point()+
#  geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+
#  theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "black", "red"))+
#  expand_limits(x=c(-8,8))
#p2=p+geom_text_repel(aes(label=rownames(allgene2$type!="normal")))
#p2
# try to only label DE genes (original, but with mistake)

  
tiff(file="20180930volcano.tif",res=100,units='in',width=10,height=10)

results <-allgene2
inputs <- cbind(gene=rownames(results), results) #convert the rownames to a column
results2 = mutate(inputs, sig=ifelse(inputs$type!="normal", "DE gene", "Non-DE"))

p=ggplot(results2,aes(x=log2FoldChange,y=-log10(padj)),colour=threshold)+
  xlab("log2FC")+ylab("-log10adj.P")+ggtitle("Volcano Plot")+
  theme(plot.title=element_text(hjust=0.5))+geom_point(aes(col=type))+
  geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+
  theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "black", "red"))+
  expand_limits(x=c(-8,8))
p
p2=p+geom_text_repel(data=filter(results2, type!="normal"), aes(label=gene))
p2
dev.off()
# try to only label DE genes, without interfering with dots or other labels


tiff(file="20180930volcano.tif",res=100,units='in',width=10,height=10)

results <-allgene2
inputs <- cbind(gene=rownames(results), results) #convert the rownames to a column
results2 = mutate(inputs, sig=ifelse(inputs$type!="normal", "DE gene", "Non-DE"))

p=ggplot(results2,aes(x=log2FoldChange,y=-log10(padj)),colour=threshold)+
  xlab("log2FC")+ylab("-log10adj.P")+ggtitle("Volcano Plot")+
  theme(plot.title=element_text(hjust=0.5))+geom_point(aes(col=type))+
  geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+
  theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "transparent", "red"))+
  expand_limits(x=c(-8,8))
p
p2=p+geom_text_repel(data=filter(results2, type!="normal"), aes(label=gene))
p2
dev.off()
# try to only label DE genes, without interfering with dots or other labels
# transparent non-DE genes



tiff(file="20180930volcano.tif",res=100,units='in',width=10,height=10)
####aim: make size of point proportional to Mean expression strength (baseMean)
ggplot(allgene2,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+xlab("log2FC")+ylab("-log10adj.P-Value")+ggtitle("Volcano Plot")+theme(plot.title=element_text(hjust=0.5))+
geom_point(size=allgene2$baseMean*0.0005)+geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+theme(panel.grid =element_blank())+scale_colour_manual(values=c("blue", "transparent", "red"))+expand_limits(x=c(-8,8))
#set dots of non-DE genes as "transparent", increase size of dots to exaggerate DE genes
dev.off()



#ggplot(allgene2,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+xlab("log2FC")+ylab("-log10adj.P-Value")+ggtitle("Volcano Plot")+theme(plot.title=element_text(hjust=0.5))+
#geom_point(size=(allgene3$baseMean*0.00001)^2)+geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+theme(panel.grid =element_blank())+scale_colour_manual(values=c("green", "black", "red"))+expand_limits(x=c(-8,8))

#ggplot(allgene2,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+xlab("log2FC")+ylab("-log10adj.P-Value")+ggtitle("Volcano Plot")+theme(plot.title=element_text(hjust=0.5))+
#geom_point(size=allgene3$baseMean*0.00005)+geom_vline(xintercept = c(-1,1))+geom_hline(yintercept = -log10(0.05))+theme(panel.grid =element_blank())+scale_colour_manual(values=c("green", "black", "red"))+expand_limits(x=c(-8,8))
#set size=allgene3$baseMean*0.00005, dots of non-DE genes won't overlap with dots representing DE genes 








##online tutorial
#https://www.biostars.org/p/203876/

# Install ggrepel package if needed
#install.packages("devtools")
#devtools::install_github("slowkow/ggrepel")
library("ggrepel") #Avoid overlapping labels
#install.packages("dplyr")
library(dplyr)

tiff(file="20180930volcano.tif",res=100,units='in',width=10,height=10)

results <-allgene2
inputs <- cbind(gene=rownames(results), results) #convert the rownames to a column
#####
#results2 = mutate(inputs, sig=ifelse(inputs$padj<0.05, "FDR<0.05", "Not Sig"))
results2 = mutate(inputs, sig=ifelse(inputs$type!="normal", "DE gene", "Non-DE"))
#p = ggplot(results2, aes(log2FoldChange, -log10(pvalue))) +
#  geom_point(aes(col=sig)) +
#  scale_color_manual(values=c("red", "black"))
p = ggplot(results2, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p
#p2=p+geom_text_repel(data=filter(results2, padj<0.05), aes(label=gene))
#p2
p3=p+geom_text_repel(data=filter(results2, type!="normal"), aes(label=gene))
p3
dev.off()






###add a column of EntrezID (will be used for GO enrichment, using clusterProfiler)
###3.多变量差异分析分析和不同基因ID之间转换

# 加载基因注释软件包
("https://bioconductor.org/biocLite.R") 
#biocLite(" AnnotationDbi")
library("AnnotationDbi")

###source("https://bioconductor.org/biocLite.R") 
####biocLite("org.Hs.eg.db")
####library("org.Hs.eg.db")
##columns(org.Hs.eg.db)

#biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
#Genome wide annotation for Mouse, primarily based on mapping using Entrez Gene identifiers.

columns(org.Mm.eg.db)



allgene$entrez <- mapIds(org.Mm.eg.db,
                         keys=row.names(res_VA),
                         keytype="SYMBOL",
                         column="ENTREZID",
                         multiVals="first")
# 以res_VA结果的row names作为key
# 设置keytype为SYMBOL
# column="ENTREZID"表示提取EntrezID信息
# multiVals="first"表示如果key对应的EntrezID有多个，提取第一个,因为通常第一个EntrezId用得较为广泛
##will show: 'select()' returned 1:many mapping between keys and columns
dim(allgene)

allgene3<- na.omit(allgene)
#将entrez一列显示NA的gene行全部移除
dim(allgene3)

write.csv(allgene3,file="20181004allgene3.csv",quote=F)
#将新的所有基因的检测报告，输出为csv格式。

gene_up<-allgene3[allgene3$type=="up",]
#筛选上调差异基因，并存储于变量gene_up中

gene_down<-allgene3[allgene3$type=="down",]
#筛选下调差异基因，并存储于变量gene_down中

gene_dif<-allgene[allgene$type!="normal",]
#筛选非差异基因，并存储于变量gene_dif中

write.csv(gene_up,file="20181004gene_up.csv",quote=F)
#将筛选的上调基因检测报告，保存为csv格式
write.csv(gene_down,file="20181004gene_down.csv",quote=F)
#将筛选的下调基因检测报告，保存为csv格式
write.csv(gene_dif,file="20181004gene_dif.csv",quote=F)
#将筛选的差异基因检测报告，保存为csv格式







#VA effect under infection status
res_VAeffect_withinfection <- results(dds2,list(c("VA_VAS_vs_VAD","VAVAS.is.infinf")))
sig.VAeffect_withinfection <- res_VAeffect_withinfection$padj<0.05
sum(sig.VAeffect_withinfection,na.rm = T)

save(file = "res_chai/de_res_R04_20180916.Rdata",res_interterm,res_infectioneffectforVAD,res_infectioneffectforVAS,res_VA,res_VAeffect_withinfection)
sessionInfo()










# ~~~ Using a grouping variable ~~~
###insert at 
#####->###

# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.

dds2$group <- factor(paste0(dds2$VA, dds2$is.inf))
design(dds2) <- ~ group
dds2 <- DESeq(dds2)
resultsNames(dds2)

# the infection effect for VAS
res_01 <- results(dds2, contrast=c("group", "VASinf", "VASN_inf"))
sig.res_01 <- res_01$padj<0.05
sum(sig.res_01,na.rm = T)


# the VA effect without infection
res_01 <- results(dds2, contrast=c("group", "VASN_inf", "VADN_inf"))
sig.res_01 <- res_01$padj<0.05
sum(sig.res_01,na.rm = T)


# the VA effect under infection
res_VAeffect_underinfection <- results(dds2, contrast=c("group", "VASinf", "VADinf"))
sig.VAeffect_underinfection <- res_VAeffect_underinfection$padj<0.05
sum(sig.VAeffect_underinfection,na.rm = T)

save(file = "res_chai/de_res_R04_VAeffectunderinfection_June2018.Rdata",res_VAeffect_underinfection)



