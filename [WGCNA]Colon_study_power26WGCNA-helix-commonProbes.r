setwd("c:/R_work/20190319 R11 power 26 WGCNA-helix-commonProbes/");
rm(list=ls(all=TRUE));

install.packages("WGCNA")


# 01. Data input and cleaning

library(WGCNA)
#配色包
library(RColorBrewer)
options(stringsAsFactors=FALSE)
#allowWGCNAThreads()

#femData <- read.csv("LiverFemale3600.csv",header=T,sep=",",check.names=F)
#head(femData[,1:13])


dat0 <- read.table(file="20181206_R11colon_highcount_normalizedbyDESeq2.txt", header = T, sep = "\t")
#dat1 <- dat0[,2:19]
#dat1 <- dat1[,-c(4,5,6,10,11,12)]
dat1 <- dat0[,2:17]
row.names(dat1) <- dat0[,1]

head(dat1[,1:12])


datExpr0=t(dat1)
dim(datExpr0)


gsg <- goodSamplesGenes(datExpr0)
## Flagging genes and samples with too many missing values...（expression level NA/zero in more thant 50% samples）
##  ..step 1
gsg$allOK
## [1] TRUE

####if "FALSE"
if(!gsg$allOK){
if(sum(!gsg$goodGenes)>0){
printFlush(paste("Removing genes:",paste(names(datExpr0)[!gsg$goodGenes],collapse=",")))
}
if (sum(!gsg$goodSamples)>0){
printFlush(paste("Removing samples:",paste(rownames(datExpr0)[!gsg$goodSamples],collapse=",")))
}
datExpr0 <- datExpr0[gsg$goodSamples,gsg$goodGenes]
}




sampleTree <- hclust(dist(datExpr0),method="average")
pdf("01.samples_cluster_tree.pdf",width=25,height=8)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="")
#cutHeight <- 15
abline(h=cutHeight,col="red")
dev.off()


datExpr <- datExpr0

traitData = read.csv("20190312-ClinicalTrait_R11.csv");
dim(traitData)
names(traitData)
head(traitData)
##       Mice Number Mouse_ID          Strain sex        DOB parents Western_Diet
## 1 1 F2_290    290    306-4 BxH ApoE-/-, F2   2 2002-03-22  229232   2002-05-14
## 2 2 F2_291    291    307-1 BxH ApoE-/-, F2   2 2002-03-22     232   2002-05-14
## 3 3 F2_292    292    307-2 BxH ApoE-/-, F2   1 2002-03-22     232   2002-05-14
## 4 4 F2_293    293    307-3 BxH ApoE-/-, F2   1 2002-03-22     232   2002-05-14
## 5 5 F2_294    294    307-4 BxH ApoE-/-, F2   1 2002-03-22     232   2002-05-14
## 6 6 F2_295    295    308-1 BxH ApoE-/-, F2   1 2002-03-22     232   2002-05-14
##     Sac_Date weight_g length_cm ab_fat other_fat total_fat comments
## 1 2002-09-11     36.9       9.9   2.53      2.26      4.79     <NA>
## 2 2002-09-11     48.5      10.7   2.90      2.97      5.87     <NA>
## 3 2002-09-11     45.7      10.4   1.04      2.31      3.35     <NA>
## 4 2002-09-11     50.3      10.9   0.91      1.89      2.80     <NA>
## 5 2002-09-11     44.8       9.8   1.22      2.47      3.69     <NA>
## 6 2002-09-11     39.2      10.2   3.06      2.49      5.55     <NA>

rownames(datTraits) = traitData[traitRows, 1];
MiceSamples = rownames(datExpr);
traitRows = match(MiceSamples, traitData$Mice);
datTraits = traitData[traitRows, -c(1,2)];

collectGarbage();




sampleTree2 <- hclust(dist(datExpr),method="average")
traitColors <- numbers2colors(datTraits,colors=blueWhiteRed(50))
pdf("01.cluster_trait.pdf",width=25,height=10)
plotDendroAndColors(sampleTree2,traitColors,
				groupLabels=names(datTraits),
				main="Sample dendrogram and trait heatmap")
dev.off()

###################################################################################################
# 02. Choosing the soft-thresholding power

powers <- c(c(1:10),seq(from=12,to=30,by=2))
#指定网络类型（邻接adjacency矩阵类型）
networkType <- "signed"
sft <- pickSoftThreshold(datExpr,
						powerVector=powers,
						networkType="signed",
						verbose=5)
# networkType = "unsigned", adjacency = |cor|^power;
#signed可以区分极大正相关和极大负相关
# networkType = "signed", adjacency = ((1+cor)/2)^power;
# networkType = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise;
# networkType = "distance", adjacency = (1-(dist/max(dist))^2)^power.

## pickSoftThreshold: will use block size 3060.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 3060 of 3060
##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k. max.k.
## 1      1   0.0945  4.49          0.840 1560.000  1560.000 1720.0
## 2      2   0.0631 -1.89          0.849  854.000   855.000 1050.0
## 3      3   0.1880 -2.18          0.911  494.000   495.000  679.0
## 4      4   0.3050 -2.17          0.930  300.000   301.000  465.0
## 5      5   0.3030 -1.74          0.942  191.000   190.000  331.0
## 6      6   0.3160 -1.51          0.949  126.000   123.000  244.0
## 7      7   0.3400 -1.40          0.948   86.500    82.900  185.0
## 8      8   0.3570 -1.35          0.944   61.100    56.800  144.0
## 9      9   0.3850 -1.29          0.944   44.300    40.300  114.0
## 10    10   0.4190 -1.25          0.948   32.900    29.500   91.0
## 11    12   0.5360 -1.26          0.979   19.200    16.600   61.1
## 12    14   0.6360 -1.70          0.961   12.000     9.840   45.9
## 13    16   0.6730 -2.05          0.957    7.830     6.120   35.3
## 14    18   0.7230 -2.09          0.968    5.330     3.910   27.6
## 15    20   0.7540 -2.06          0.959    3.760     2.590   22.0
## 16    22   0.8790 -1.81          0.999    2.730     1.750   17.7
## 17    24   0.8970 -1.88          0.982    2.030     1.220   15.8
## 18    26   0.9180 -1.92          0.969    1.550     0.856   14.4
## 19    28   0.9140 -1.96          0.939    1.200     0.615   13.3
## 20    30   0.9180 -1.93          0.917    0.951     0.442   12.4

pdf(file="02.soft_threshold.pdf",width=8,height=4.5)
par(mfrow=c(1,2))
cex1 <- 0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft threshold (power)",
	ylab="Scale free topology model fit, signed R^2",
	type="n",
	main=paste("Scale independence"))
	text(sft$fitIndices[,1],
	-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1],sft$fitIndices[,5],
xlab="Soft threshold (power)",ylab="Mean connectivity",
type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")
dev.off()

###################################################################################################
# 03. One-step network construction and module detection

softPower <- 26
networkType <- "signed"
net <- blockwiseModules(datExpr,power=softPower,
					networkType=networkType,numericLabels=TRUE,
					mergeCutHeight=0.25,minModuleSize=30,
					maxBlockSize=30000,saveTOMs=TRUE,
					saveTOMFileBase="WGCNA_TOM",verbose=5)
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
pdf(file="03.auto_modules_color.pdf",width=12,height=4)
plotDendroAndColors(net$dendrograms[[1]],
				moduleColors[net$blockGenes[[1]]],"Module colors",
				dendroLabels=FALSE,addGuide=TRUE)
dev.off()

save(datExpr,sft,softPower,networkType,net,moduleColors,
	file="WGCNA_2.RData")





##Yafei output##
mods <- moduleEigengenes(datExpr, colors = moduleColors)

table(mods$validColors)
save(file = "WGCNA_mod_R11_15341gene_power26_signed_20190319.Rdata",mods)
##adjacency,TOM

## print as tables##
for (i in unique(mods$validColors)){
  write.table(file=paste("mod_",i,".txt",sep=""),colnames(datExpr)[which(mods$validColors == i)],quote = F,col.names = F,row.names = F)
}










###################################################################################################
# 04. Step-by-step network construction and module detection

softPower <- 22
networkType <- "signed"
adjacency <- adjacency(datExpr,power=softPower,type=networkType)
adjacency[1:5,1:5]
##              MMT00000044  MMT00000046  MMT00000051  MMT00000080  MMT00000102
## MMT00000044 1.000000e+00 1.017339e-07 1.695517e-06 3.969534e-06 1.011784e-06
## MMT00000046 1.017339e-07 1.000000e+00 2.801347e-15 5.707620e-08 1.315773e-15
## MMT00000051 1.695517e-06 2.801347e-15 1.000000e+00 6.993674e-08 1.313098e-06
## MMT00000080 3.969534e-06 5.707620e-08 6.993674e-08 1.000000e+00 3.334227e-08
## MMT00000102 1.011784e-06 1.315773e-15 1.313098e-06 3.334227e-08 1.000000e+00

TOM <- TOMsimilarity(adjacency,TOMType="signed")
dissTOM <- 1-TOM
dissTOM[1:5,1:5]
##           [,1]      [,2]      [,3]      [,4]      [,5]
## [1,] 0.0000000 0.9999996 0.9999810 0.9999815 0.9999964
## [2,] 0.9999996 0.0000000 1.0000000 0.9999928 1.0000000
## [3,] 0.9999810 1.0000000 0.0000000 0.9999909 0.9997470
## [4,] 0.9999815 0.9999928 0.9999909 0.0000000 0.9999983
## [5,] 0.9999964 1.0000000 0.9997470 0.9999983 0.0000000

geneTree <- hclust(as.dist(dissTOM),method="average")
pdf(file="04.man_genes_cluster_tree.pdf",width=12,height=4)
plot(geneTree,xlab="",
	sub="",main="Cluster dendrogram on TOM-based dissimilarity",
	labels=FALSE)
dev.off()

dynamicMods <- cutreeDynamic(dendro=geneTree,distM=dissTOM,
							deepSplit=2,cutHeight=0.995,minClusterSize=30)
dynamicColors <- labels2colors(dynamicMods)
pdf(file="04.man_modules_color.pdf",width=12,height=4)
plotDendroAndColors(geneTree,
				dynamicColors,"Module colors",
				dendroLabels=FALSE,addGuide=TRUE)
dev.off()

MEDissThres <- 0.25
merge <- mergeCloseModules(datExpr,dynamicColors,cutHeight=MEDissThres)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

pdf(file="04.comparison_modules_color.pdf",width=12,height=5)
plotDendroAndColors(geneTree,
				cbind(moduleColors,dynamicColors,mergedColors),
				c("Auto construction","Dynamic tree cut", "Merged dynamic"),
				dendroLabels=FALSE,addGuide=TRUE)
dev.off()











###################################################################################################
# 05. Module eigengenes Calculation

MEs0 <- moduleEigengenes(datExpr[,net$goodGenes],
						moduleColors[net$blockGenes[[1]]])$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames(datExpr[,net$goodGenes])
text <- cbind(rownames(MEs),MEs)
colnames(text)[1] <- "samples"
write.table(text,file="05.module_eigengenes.xls",
		quote=F,sep="\t",row.names=F)

names(MEs) <- substring(names(MEs),3)
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss),method="average")
pdf(file="05.modules_cluster_tree.pdf",width=7,height=5)
plot(METree,main="Clustering of module eigengenes",xlab="",sub="")
dev.off()

moduleCor <- corAndPvalue(MEs,use="p")
rowLabels <- paste("ME",names(MEs),sep="")
textMatrix <- paste(signif(moduleCor$cor,2),
				"\n(",signif(moduleCor$p,1),")",sep="")
dim(textMatrix) <- dim(moduleCor$cor)
pdf(file="05.modules_relationships.pdf",12,8)
par(mar=c(10,10,1,2))
labeledHeatmap(Matrix=moduleCor$cor,
			textMatrix=textMatrix,
			xLabels=rowLabels,yLabels=rowLabels,
			xSymbols=names(MEs),ySymbols=names(MEs),
			colorLabels=TRUE,colors=blueWhiteRed(50),
			setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

text <- paste("cor=",round(moduleCor$cor,4),
				";p-value=",round(moduleCor$p,4),sep="")
dim(text) <- dim(moduleCor$cor)
rownames(text) <- rowLabels
colnames(text) <- rowLabels
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="05.modules_relationships.xls",
		quote=F,sep="\t",row.names=F)

dir.create("05.expression_ME")
for(i in 1:(ncol(MEs)-1)) {
which.module <- labels2colors(i)
dir <- "05.expression_ME/"
pdf(file=paste(dir,"05.expression_ME_",which.module,".pdf",sep=""),
	25,10)
ME <- MEs[,which.module]
ME <- t(as.matrix(MEs[,which.module]))
colnames(ME) <- rownames(datExpr[,net$goodGenes])
layout(matrix(c(1,2)),heights=c(1.5,3))
par(mar=c(0.3,9,3,5))
plotMat(t(scale(datExpr[,net$goodGenes][,moduleColors[net$blockGenes[[1]]]==which.module])),
	nrgcols=30,rlabels=F,rcols=which.module,
	main=paste(which.module),cex.main=1)
###colors=blueWhiteRed(50),
###color=colorRampPalette(c("blue","white","red"))
par(mar=c(5,4,0,1))
barplot(ME,col=which.module,main="",cex.names=1,cex.axis=1,
ylab="module eigengene",las=3)
dev.off()
}


###################################################################################################
# 06. Relating modules to external information (samples/traits)

sample_cor <- cor(t(datExpr[,net$goodGenes]),t(datExpr[,net$goodGenes]),
				use='pairwise.complete.obs')
moduleSampleCor <- cor(MEs,sample_cor,use="p")
nSamples <- nrow(datExpr[,net$goodGenes])
moduleSamplePvalue <- corPvalueStudent(moduleSampleCor,nSamples)
textMatrix <- paste(signif(moduleSampleCor,2),
				"\n(",signif(moduleSamplePvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleSampleCor)
rowLabels <- paste("ME",names(MEs),sep="")
pdf(file="06.modules_samples_relationships.pdf",22,6)
par(mar=c(5,12,1,1))
labeledHeatmap(Matrix=moduleSampleCor,
			xLabels=colnames(sample_cor),
			yLabels=rowLabels,ySymbols=names(MEs),
			colorLabels=TRUE,colors=blueWhiteRed(50),
			setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

text <- paste("cor=",round(moduleSampleCor,4),
				";p-value=",round(moduleSamplePvalue,4),sep="")
dim(text) <- dim(moduleSampleCor)
rownames(text) <- rownames(moduleSampleCor)
colnames(text) <- colnames(moduleSampleCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="06.modules_samples_relationships.xls",
		quote=F,sep="\t",row.names=F)



#reload data for visualization
load(file="WGCNA_2.RData")
#load datTrait
traitData = read.csv("20190312-ClinicalTrait_R11.csv");
dim(traitData)
names(traitData)
head(traitData)
rownames(datTraits) = traitData[traitRows, 1];
MiceSamples = rownames(datExpr);
traitRows = match(MiceSamples, traitData$Mice);
datTraits = traitData[traitRows, -c(1,2)];

collectGarbage();




moduleTraitCor <- cor(MEs,datTraits,use="p")
nSamples <- nrow(datExpr[,net$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste(signif(moduleTraitCor,2),
				"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- paste("ME",names(MEs),sep="")
pdf(file="06.modules_traits_relationships_bigFont.pdf",15,11.56)
par(mar=c(12,11,1,2))
labeledHeatmap(Matrix=moduleTraitCor,
			textMatrix=textMatrix,
			xLabels=colnames(datTraits),
			yLabels=rowLabels,ySymbols=names(MEs),
			colorLabels=TRUE,colors=blueWhiteRed(50),
			setStdMargins=FALSE,xLabelsAngle=45,zlim=c(-1,1),
			cex.lab.x = 2, cex.lab.y = 2)
dev.off()

text <- paste("cor=",round(moduleTraitCor,4),
				";p-value=",round(moduleTraitPvalue,4),sep="")
dim(text) <- dim(moduleTraitCor)
rownames(text) <- rownames(moduleTraitCor)
colnames(text) <- colnames(moduleTraitCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="06.modules_traits_relationships.xls",
		quote=F,sep="\t",row.names=F)

###################################################################################################
# 07. Relating nodes to external information (samples/traits)

modNames <- names(MEs)
geneModuleMembership <- cor(datExpr[,net$goodGenes],MEs,use="p")
nSamples <- nrow(datExpr[,net$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership,nSamples)
colnames(geneModuleMembership) <- paste("MM",modNames,sep="")
colnames(MMPvalue) <- paste("p.MM",modNames,sep="")

text <- paste("cor=",round(geneModuleMembership,4),
				";p-value=",round(MMPvalue,4),sep="")
dim(text) <- dim(geneModuleMembership)
rownames(text) <- rownames(geneModuleMembership)
colnames(text) <- colnames(geneModuleMembership)
text <- cbind(rownames(text),moduleColors[net$blockGenes[[1]]],text)
colnames(text)[1] <- "modules"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_module_membership.xls",
		quote=F,sep="\t",row.names=F)








modNames <- names(MEs)
VA_status <- as.data.frame(datTraits$VA_status)
names(VA_status) <- "VA_status"
geneTraitSignificance <- as.data.frame(cor(datExpr[,net$goodGenes],VA_status,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(VA_status),sep="")
names(GSPvalue) <- paste("p.GS.",names(VA_status),sep="")
head(geneTraitSignificance)
##               GS.trigly
## MMT00000044  0.01253338
## MMT00000046  0.13950675
## MMT00000051 -0.01999059
## MMT00000080 -0.05179917
## MMT00000102 -0.06059991
## MMT00000149 -0.23484857
head(GSPvalue)
##             p.GS.trigly
## MMT00000044 0.885713080
## MMT00000046 0.107911902
## MMT00000051 0.818664712
## MMT00000080 0.552244343
## MMT00000102 0.486703013
## MMT00000149 0.006306398

dir.create("07.MM_vs_VA_status")
for(i in 1:(ncol(MEs)-1)) {
which.module <- labels2colors(i)
column <- match(which.module,modNames)
moduleGenes <- moduleColors[net$blockGenes[[1]]]==which.module
dir <- "07.MM_vs_VA_status/"
pdf(file=paste(dir,"07.",which.module,"_MM_vs_VA_status.pdf",sep=""),6,6)
verboseScatterplot(geneModuleMembership[moduleGenes,column],
				geneTraitSignificance[moduleGenes,1],
				xlab=paste("Module membership (MM) in",which.module,"module"),
				ylab="Gene significance for VA_status",
				main=paste("Module membership vs. gene significance\n"),
				col=which.module)
dev.off()
}

text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),moduleColors[net$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance_VA_status.xls",
		quote=F,sep="\t",row.names=F)

save(datExpr,sft,softPower,networkType,net,moduleColors,MEs,
	moduleSampleCor,moduleSamplePvalue,
	moduleTraitCor,moduleTraitPvalue,
	geneModuleMembership,MMPvalue,
	file="WGCNA_2.RData")



#####################################################3
############### #####################7.Infection_status
modNames <- names(MEs)
Infection_status <- as.data.frame(datTraits$Infection_status)
names(Infection_status) <- "Infection_status"
geneTraitSignificance <- as.data.frame(cor(datExpr[,net$goodGenes],Infection_status,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(Infection_status),sep="")
names(GSPvalue) <- paste("p.GS.",names(Infection_status),sep="")
head(geneTraitSignificance)
##               GS.Infection_status
## MMT00000044  0.01253338

head(GSPvalue)
##             p.GS.Infection_status
## MMT00000044 0.885713080


dir.create("07.MM_vs_Infection_status")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColors[net$blockGenes[[1]]]==which.module
  dir <- "07.MM_vs_Infection_status/"
  pdf(file=paste(dir,"07.",which.module,"_MM_vs_Infection_status.pdf",sep=""),6,6)
  verboseScatterplot(geneModuleMembership[moduleGenes,column],
                     geneTraitSignificance[moduleGenes,1],
                     xlab=paste("Module membership (MM) in",which.module,"module"),
                     ylab="Gene significance for Infection_status",
                     main=paste("Module membership vs. gene significance\n"),
                     col=which.module)
  dev.off()
}
text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),moduleColors[net$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-Infection_status.xls",
            quote=F,sep="\t",row.names=F)

##########################################
#########################################3

###################################################################################################
# 08. Exporting network

load(file="WGCNA_TOM-block.1.RData")
load(file="WGCNA_2.RData")
#save(datExpr,sft,softPower,networkType,net,moduleColors,MEs,
#     moduleSampleCor,moduleSamplePvalue,
#     moduleTraitCor,moduleTraitPvalue,
#     geneModuleMembership,MMPvalue,
#     file="WGCNA_2.RData")
ATOM <- as.matrix(TOM)
TOM1 <- ATOM[1:round((nrow(ATOM)/2)),1:round((nrow(ATOM)/2))]
TOM2 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),1:round((nrow(ATOM)/2))]
TOM3 <- ATOM[1:round((nrow(ATOM)/2)),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
TOM4 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),(round(nrow(ATOM)/2)+1):nrow(ATOM)]

dir.create("08.module_result__")
setwd("08.module_result")
for(i in 1:(ncol(MEs)-1)) {
module <- labels2colors(i)
inModule <- moduleColors[net$blockGenes[[1]]]==module
genename <- colnames(datExpr[,net$goodGenes])
modGenes <- genename[inModule]
modTOM1 <- TOM1[inModule[1:round((nrow(ATOM)/2))],
				inModule[1:round((nrow(ATOM)/2))]]
modTOM2 <- TOM2[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)],
				inModule[1:round((nrow(ATOM)/2))]]
modTOM3 <- TOM3[inModule[1:round((nrow(ATOM)/2))],
				inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
modTOM4 <- TOM4[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)],
				inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
modTOM <- rbind(cbind(modTOM1,modTOM3),cbind(modTOM2,modTOM4))

IMConn <- softConnectivity(datExpr[,net$goodGenes][,modGenes])


####weight=TOM value
cyt1 <- exportNetworkToCytoscape(modTOM,
		edgeFile=paste("CytoscapeInput-edges-",paste(module,collapse="-"),".txt",sep=""),
		nodeFile=paste("CytoscapeInput-nodes-",paste(module,collapse="-"),".txt",sep=""),
		weighted=TRUE,threshold=0.02,
		nodeNames=modGenes,altNodeNames=modGenes,
		nodeAttr=moduleColors[net$blockGenes[[1]]][inModule])

####sum of the adjacency to the other nodes within the same module->connectivity of a node-> rank(top 5%-15%) -> Hub edge (key regulatory role in the module)
out <- cbind(modGenes,IMConn)
###connectivity=sum of adjacency with all module members in the same module
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(module,"-module-gene.txt",sep=""),
			sep="\t",quote=F,row.names=F)


###hub edge extraction
nTop <- 0.05*length(modGenes)
top <- (rank(-IMConn) <= nTop)
out <- cbind(modGenes[top],IMConn[top])
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(module,"-5%hubgene.txt",sep=""),
			sep="\t",quote=F,row.names=F)

nTop <- 0.1*length(modGenes)
top <- (rank(-IMConn) <= nTop)
out <- cbind(modGenes[top],IMConn[top])
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(module,"-10%hubgene.txt",sep=""),
			sep="\t",quote=F,row.names=F)

nTop <- 25
top <- (rank(-IMConn) <= nTop)
out <- cbind(modGenes[top],IMConn[top])
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(module,"-top25_hubgene.txt",sep=""),
			sep="\t",quote=F,row.names=F)
}







genename <- colnames(datExpr[,net$goodGenes])
modGenes <- genename
IMConn <- softConnectivity(datExpr[,net$goodGenes][,modGenes])
cyt1 <- exportNetworkToCytoscape(ATOM,
		edgeFile=paste("CytoscapeInput-edges-overall.txt",sep=""),
		nodeFile=paste("CytoscapeInput-nodes-overall.txt",sep=""),
		weighted=TRUE,threshold=0.02,
		nodeNames=modGenes,altNodeNames=modGenes,
		nodeAttr=moduleColors[net$blockGenes[[1]]])
out <- cbind(modGenes,IMConn)
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste("overall_gene_connectivity.txt",sep=""),
			sep="\t",quote=F,row.names=F)





##intramodularConnectivity.=connectivity of nodes to other nodes within the same module.
IntraMConn <- intramodularConnectivity.fromExpr(datExpr[,net$goodGenes],
											moduleColors[net$blockGenes[[1]]],
											networkType=networkType,
											power=softPower)
merged_IntraMConn <- cbind(colnames(datExpr[,net$goodGenes]),
		moduleColors[net$blockGenes[[1]]],IntraMConn)
colnames(merged_IntraMConn)[1:2] <- c("gene","module")
head(merged_IntraMConn)
##         gene module     kTotal  kWithin      kOut    kDiff
## 1 MMT00000044   grey 0.03598708       NA        NA       NA
## 2 MMT00000046 yellow 4.33851798 3.473299 0.8652190 2.608080
## 3 MMT00000051  brown 3.03831719 2.138751 0.8995664 1.239184
## 4 MMT00000080  green 3.21982738 2.839456 0.3803716 2.459084
## 5 MMT00000102   grey 1.03882912       NA        NA       NA
## 6 MMT00000149   blue 6.41916108 5.849061 0.5701003 5.278960
write.table(merged_IntraMConn,file="intraModular_connectivity.xls",
		quote=F,sep="\t",row.names=FALSE)

