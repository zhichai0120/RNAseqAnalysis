#dir.create("helix-WGCNA-lesson5")
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller")
rm(list=ls(all=TRUE));
#fileUrl <- "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/metaAnalysisFiles.zip"
#download.file(fileUrl,destfile="metaAnalysisFiles.zip")
#unzip("metaAnalysisFiles.zip")
#list.files()
## [1] "__MACOSX" "metaAnalysisFiles" "metaAnalysisFiles.zip"

#setwd("metaAnalysisFiles")
#load("metaAnalysisData.RData")
#ls()
## [1] "datExprA1" "datExprA2" "datExprB1" "datExprB2" "genesA"    "genesI"   
## [7] "probesA"   "probesI"

###################################################################################################

library(WGCNA)
library(reshape2)


dat0 <- read.table(file="20180914_R04LSI_highcount_normalizedbyDESeq2.txt", header = T, sep = "\t")
#dat1 <- dat0[,2:19]
#dat1 <- dat1[,-c(4,5,6,10,11,12)]
dat1 <- dat0[,2:13]
row.names(dat1) <- dat0[,1]
datExprA1 <-dat1

dat3 <- read.table(file="20181206_R11colon_highcount_normalizedbyDESeq2.txt", header = T, sep = "\t")
#dat1 <- dat0[,2:19]
#dat1 <- dat1[,-c(4,5,6,10,11,12)]
dat4 <- dat3[,2:17]
row.names(dat4) <- dat3[,1]
datExprA2 <-dat4

commonProbesA <- intersect(rownames(datExprA1),rownames(datExprA2))
datExprA1p <- datExprA1[commonProbesA,]
datExprA2p <- datExprA2[commonProbesA,]

dim(datExprA1p)
#[1] 14059    12
dim(datExprA2p)
#[1] 14059    16

softPower1 <- 7
softPower2<- 26
rankExprA1 <- rank(rowMeans(datExprA1p))
rankExprA2 <- rank(rowMeans(datExprA2p))
random5000 <- sample(row.names(datExprA1p),5000)
rankConnA1 <- rank(softConnectivity(t(datExprA1p[random5000,]),
				type="signed",power=softPower1))
rankConnA2 <- rank(softConnectivity(t(datExprA2p[random5000,]),
				type="signed",power=softPower2))

pdf("01.generalNetworkProperties.pdf",height=5,width=10)
par(mfrow=c(1,2))
verboseScatterplot(rankExprA1,rankExprA2,xlab="Ranked Expression (SI)",
				ylab="Ranked Expression (Colon)")
verboseScatterplot(rankConnA1,rankConnA2,xlab="Ranked Connectivity (SI)",
				ylab="Ranked Connectivity (Colon)")
dev.off()

###################################################################################################

networkType <- "signed"
netA1 <- blockwiseModules(t(datExprA1p),power=softPower1,
				networkType=networkType,numericLabels=TRUE,
				mergeCutHeight=0.25,minModuleSize=30,
				maxBlockSize=30000,saveTOMs=TRUE,
				saveTOMFileBase="WGCNA_TOM_A1",verbose=5)
moduleLabelsA1 <- netA1$colors
moduleColorsA1 <- labels2colors(moduleLabelsA1)

##lnames = load(file = "20190218-R04-16075gene_power6-input cleaned as Tutorial 2.1-FemaleLiver-02-networkConstruction-BlockWise.RData");
##lnames
##moduleLabelsA1 <- BWmoduleLabels
##moduleColorsA1 <- BWmoduleColors

netA2 <- blockwiseModules(t(datExprA2p),power=softPower2,
                networkType=networkType,numericLabels=TRUE,
                mergeCutHeight=0.25,minModuleSize=30,
                maxBlockSize=30000,saveTOMs=TRUE,
                saveTOMFileBase="WGCNA_TOM_A2",verbose=5)
moduleLabelsA2 <- netA2$colors
moduleColorsA2 <- labels2colors(moduleLabelsA2)

#lnames = load(file = "20190220-R11-16075gene_power6-MaleLiver-02-networkConstruction-BlockWise.RData");
#lnames
#moduleLabelsA2 <- BWmoduleLabels
#moduleColorsA2 <- BWmoduleColors
#lnames = load(file = "WGCNA_mod_R11_16075gene_power6_20190220.RData");
#lnames

pdf(file="02.comparison_modules_color.pdf",width=12,height=5)
plotDendroAndColors(netA1$dendrograms[[1]],
				cbind(moduleColorsA1[netA1$blockGenes[[1]]],
					moduleColorsA2[netA2$blockGenes[[1]]]),
				c("Modules SI","Modules Colon"),
				dendroLabels=FALSE,addGuide=TRUE)
dev.off()

###################################################################################################

multiExpr <- list(A1=list(data=t(datExprA1p)),A2=list(data=t(datExprA2p)))
multiColor <- list(A1=moduleColorsA1[netA1$blockGenes[[1]]])
mp <- modulePreservation(multiExpr,multiColor,
				referenceNetworks=1,verbose=5,
				networkType="signed",
				nPermutations=50)
stats <- mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)]
##              moduleSize Zsummary.pres
## brown              1000      54.22666
## blue               1000      41.47798
## black               638      36.43791
## turquoise          1000      36.06894
## yellow             1000      33.92899
stats[order(-stats[,2])]

###################################################################################################

MEs0_A1 <- moduleEigengenes(t(datExprA1p)[,netA1$goodGenes],
						moduleColorsA1[netA1$blockGenes[[1]]])$eigengenes
MEs_A1 <- orderMEs(MEs0_A1)
geneModuleMembership1 <- signedKME(t(datExprA1p)[,netA1$goodGenes],MEs_A1)
name <- colnames(MEs_A1)
colnames(geneModuleMembership1) <- paste("PC",name,".cor",sep="")
geneModuleMembership1 <-
geneModuleMembership1[,order(colnames(geneModuleMembership1))]
MMPvalue1 <- corPvalueStudent(as.matrix(geneModuleMembership1),
			dim(datExprA1p)[[2]])
colnames(MMPvalue1) <- paste("PC",name[order(name)],".pval",sep="")
Gene1 <- rownames(datExprA1p)[netA1$goodGenes]
Gene1 <- as.data.frame(Gene1)
kMEtable1 <- cbind(Gene1,Gene1,moduleColorsA1[netA1$blockGenes[[1]]])
for (i in 1:nrow(t(MEs_A1))){
kMEtable1 <- cbind(kMEtable1,geneModuleMembership1[,i],MMPvalue1[,i])
}
colnames(kMEtable1) <- c("PSID","Gene","Module",
		sort(c(colnames(geneModuleMembership1),colnames(MMPvalue1))))
kMEtable1 <- kMEtable1[kMEtable1$Module!="grey",]
write.csv(na.omit(kMEtable1),"03.kMEtable_A1.csv",row.names=FALSE)
#write.csv(kMEtable1,"03.no omit-kMEtable_A1.csv",row.names=FALSE)

MEs0_A2 <- moduleEigengenes(t(datExprA2p)[,netA2$goodGenes],
						moduleColorsA2[netA2$blockGenes[[1]]])$eigengenes
MEs_A2 <- orderMEs(MEs0_A2)
geneModuleMembership2 <- signedKME(t(datExprA2p)[,netA2$goodGenes],MEs_A2)
name <- colnames(MEs_A2)
colnames(geneModuleMembership2) <- paste("PC",name,".cor",sep="")
geneModuleMembership2 <-
geneModuleMembership2[,order(colnames(geneModuleMembership2))]
MMPvalue2 <- corPvalueStudent(as.matrix(geneModuleMembership2),
			dim(datExprA2p)[[2]])
colnames(MMPvalue2) <- paste("PC",name[order(name)],".pval",sep="")
Gene2 <- rownames(datExprA2p)[netA2$goodGenes]
Gene2 <- as.data.frame(Gene2)
kMEtable2 <- cbind(Gene2,Gene2,moduleColorsA2[netA2$blockGenes[[1]]])
for (i in 1:nrow(t(MEs_A2))){
kMEtable2 <- cbind(kMEtable2,geneModuleMembership2[,i],MMPvalue2[,i])
}
colnames(kMEtable2) <- c("PSID","Gene","Module",
		sort(c(colnames(geneModuleMembership2),colnames(MMPvalue2))))
kMEtable2 <- kMEtable2[kMEtable2$Module!="grey",]
write.csv(na.omit(kMEtable2),"03.kMEtable_A2.csv",row.names=FALSE)   
#write.csv(kMEtable2,"03.no omit-kMEtable_A2.csv",row.names=FALSE)





#######===========20190427 to re-make heatmap with bigger font size for x and y labels
#reload data for visualization
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190427 big font")
rm(list=ls(all=TRUE));
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190320_individual_R04")
load(file="WGCNA_2_individualR04.RData")
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190320_individual_R11")
load(file="WGCNA_2_individualR11.RData")
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190427 big font")
library(WGCNA)
library(reshape2)

datExprA1p <- t(datExprA1q)
datExprA2p <- t(datExprA2q)
##########===================


enrichments_A1_A2 <- userListEnrichment(rownames(datExprA1p)[netA1$goodGenes],
				moduleColorsA1[netA1$blockGenes[[1]]],
				"03.kMEtable_A2.csv",
				"","03.significant_enrichment_A1_A2.csv")
enrichments_A1_A2$pValues[,2] <- sub("__","",enrichments_A1_A2$pValues[,2])
head(enrichments_A1_A2$pValues)
##   InputCategories UserDefinedCategories Type NumOverlap      Pvalues
## 1           black                 black User        148 1.771052e-53
## 2           black                  blue User        260 6.720260e-61
## 3           black                 brown User          2 1.000000e+00
## 4           black                  cyan User          2 9.997414e-01
## 5           black                 green User          3 1.000000e+00
## 6           black           greenyellow User          1 1.000000e+00
##   CorrectedPvalues
## 1     6.021578e-51
## 2     2.284888e-58
## 3     1.000000e+00
## 4     1.000000e+00
## 5     1.000000e+00
## 6     1.000000e+00
write.csv(enrichments_A1_A2$pValues, file = "20190509_enrichment_A1_A2.csv",quote = F)
ovGenes=enrichments_A1_A2$ovGenes
save(ovGenes,
     file="20190509_enrichment_overlapGenes.RData")

load("20190319_enrichment_overlapGenes.Rdata")

### e.g. 
enrichments_A1_A2$ovGenes$'blue -- magenta'
##OR
ovGenes$'blue -- blue'
ovGenes$'black -- greenyellow'
ovGenes$'red -- greenyellow'

###all significant overlaps
ovGenes$'black -- brown'
ovGenes$'black -- greenyellow'

ovGenes$'blue -- magenta'
ovGenes$'blue -- turquoise'
write.table(ovGenes$'blue -- turquoise',file="ovGenes'blue -- turquoise'.txt",sep="\t",quote=F,row.names=F)

ovGenes$'brown -- red'
ovGenes$'brown -- yellow'

ovGenes$'green -- green'
ovGenes$'green -- tan'
ovGenes$'green -- turquoise'

ovGenes$'red -- brown'
ovGenes$'red -- greenyellow'

write.table(ovGenes$'turquoise -- blue',file="ovGenes'turquoise -- blue'.txt",sep="\t",quote=F,row.names=F)
write.table(ovGenes$'turquoise -- brown',file="ovGenes'turquoise -- brown'.txt",sep="\t",quote=F,row.names=F)
write.table(ovGenes$'turquoise -- yellow',file="ovGenes'turquoise -- yellow'.txt",sep="\t",quote=F,row.names=F)


ovGenes$'yellow -- green'
ovGenes$'yellow -- pink'
write.table(ovGenes$'yellow -- turquoise',file="ovGenes'yellow -- turquoise'.txt",sep="\t",quote=F,row.names=F)



N_A1_A2 <- dcast(enrichments_A1_A2$pValues,
			InputCategories~UserDefinedCategories,
			value.var="NumOverlap")
rownames(N_A1_A2) <- N_A1_A2[,1]
N_A1_A2 <- N_A1_A2[,-1]
N_A1_A2[1:5,1:5]
##       black blue brown cyan green
## black   148  260     2    2     3
## blue     21    9  1460    3   114
## brown     9    7   633   15   953
## cyan      6   48     0   74     5
## green     8    9   297   18   140
N_A1_A2
#t(N_A1_A2)


P_A1_A2 <- dcast(enrichments_A1_A2$pValues,
			InputCategories~UserDefinedCategories,
			value.var="Pvalues")
rownames(P_A1_A2) <- P_A1_A2[,1]
P_A1_A2 <- P_A1_A2[,-1]

CP_A1_A2 <- dcast(enrichments_A1_A2$pValues,
			InputCategories~UserDefinedCategories,
			value.var="CorrectedPvalues")
rownames(CP_A1_A2) <- CP_A1_A2[,1]
CP_A1_A2 <- CP_A1_A2[,-1]
LCP_A1_A2 <- -log10(CP_A1_A2)
LCP_A1_A2[LCP_A1_A2>=50] <- 50

textMatrix <- paste(as.matrix(N_A1_A2),
				"\n(",as.matrix(signif(CP_A1_A2,1)),")",sep="")

#pdf(file="03.t2-internetwork_modules_relationships_bigFont.pdf",15,10)
par(mar=c(8,8,1,2))
labeledHeatmap(Matrix=(LCP_A1_A2),
			textMatrix=textMatrix,
			xLabels=colnames(LCP_A1_A2),yLabels=rownames(LCP_A1_A2),
			xSymbols=paste("ME",colnames(LCP_A1_A2),":____",sep=""),
			ySymbols=paste("ME",rownames(LCP_A1_A2),":____",sep=""),
			colorLabels=TRUE,colors=blueWhiteRed(100)[51:100],
			setStdMargins=FALSE,xLabelsAngle=45,zlim=c(0,50),
			cex.lab.x = 2, cex.lab.y = 2)
dev.off()


pdf(file="03.t2-internetwork_modules_relationships_bigFont_SI rownames changed.pdf",15,10)
par(mar=c(15,18,1,2))
labeledHeatmap(Matrix=(LCP_A1_A2),
               textMatrix=textMatrix,
               xLabels=c("mpColon(Brown)  351","mpColon(Turquoise)  4485","mpColon(Black2)  1387","mpColon(Yellow)  512","mpColon(Black1)  94","mpColon(Magenta)  186","mpColon(Pink)  227","mpColon(Purple)  113","mpColon(Red)  409","mpColon(Tan)  56","mpColon(Blue)  4972","mpColon(Green)  682"),
               yLabels=c("mpSI(Pink)  303","mpSI(Blue)  3610","mpSI(Brown)  1700","mpSI(Yellow)  430","mpSI(Black)  179","mpSI(Red)  338","mpSI(Turquoise)  5079","mpSI(Green)  1122"),
               xSymbols=paste("ME",colnames(LCP_A1_A2),":____",sep=""),
               ySymbols=paste("ME",rownames(LCP_A1_A2),":____",sep=""),
               colorLabels=TRUE,colors=blueWhiteRed(100)[51:100],
               setStdMargins=FALSE,xLabelsAngle=45,zlim=c(0,50),
               cex.lab.x = 2, cex.lab.y = 2)
dev.off()







dir.create("04.kMEtable1_vs_kMEtable2")
topGenesKME <- NULL
topNames <- NULL
for(i in 1:ncol(geneModuleMembership1)){
for(j in 1:ncol(geneModuleMembership2)){
M1 <- sub("PCME","",names(geneModuleMembership1)[i])
M1 <- sub(".cor","",M1)
M2 <- sub("PCME","",names(geneModuleMembership2)[j])
M2 <- sub(".cor","",M2)
inMod1 <- moduleColorsA1[netA1$blockGenes[[1]]]==M1
inMod2 <- moduleColorsA2[netA2$blockGenes[[1]]]==M2
inModp <- intersect(rownames(geneModuleMembership1)[inMod1],
rownames(geneModuleMembership2)[inMod2])

if(length(inModp)>10){
fn <- paste("04.kMEtable1_",M1,"_vs_kMEtable2_",M2,".pdf",sep="")
pdf(paste("04.kMEtable1_vs_kMEtable2/",fn,sep=""),height=8,width=8)
verboseScatterplot(geneModuleMembership1[inModp,i],
geneModuleMembership2[inModp,j],
main=paste("A1_",M1,"_vs_A2_",M2,sep=""),
xlab=paste("kME in A1 ",M1,sep=""),
ylab=paste("kME in A2 ",M2,sep=""))
dev.off()

kMErank1 <- rank(-geneModuleMembership1[inModp,i])
kMErank2 <- rank(-geneModuleMembership2[inModp,j])
maxKMErank <- rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
topGenesKME <- cbind(topGenesKME,inModp[maxKMErank<=10])
comp <- paste("A1_",M1,"_vs_A2_",M2,sep="")
topNames <- cbind(topNames,comp)
}
}
}
colnames(topGenesKME) <- topNames
write.table(topGenesKME,
	"04.shared_top10_hub_internetwork_modules.xls",
	quote=F,row.names=F,sep="\t")
























setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190320_individual_R04")

gsg <- goodSamplesGenes(t(datExprA1p))
## Flagging genes and samples with too many missing values...（expression level NA/zero in more thant 50% samples）
##  ..step 1
gsg$allOK
## [1] TRUE

sampleTree <- hclust(dist(t(datExprA1p)),method="average")
pdf("01.samples_cluster_tree.pdf",width=25,height=8)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="")
#cutHeight <- 15
abline(h=cutHeight,col="red")
dev.off()


datExprA1q <- t(datExprA1p)

traitData = read.csv("ClinicalTraits_R04.csv");
dim(traitData)
names(traitData)
head(traitData)

rownames(datTraits) = traitData[traitRows, 1];
#Error in `[.data.frame`(traitData, traitRows, 1) : 
#  object 'traitRows' not found
MiceSamples = rownames(datExprA1q);
traitRows = match(MiceSamples, traitData$Mice);
datTraits = traitData[traitRows, -c(1,2)];

collectGarbage();


sampleTree2 <- hclust(dist(datExprA1q),method="average")
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
sft <- pickSoftThreshold(datExprA1q,
                         powerVector=powers,
                         networkType="signed",
                         verbose=5)
# networkType = "unsigned", adjacency = |cor|^power;
#signed可以区分极大正相关和极大负相关
# networkType = "signed", adjacency = ((1+cor)/2)^power;
# networkType = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise;
# networkType = "distance", adjacency = (1-(dist/max(dist))^2)^power.

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
# 03. One-step network construction and module detection--after

#moduleLabelsA1 <- netA1$colors
#moduleColorsA1 <- labels2colors(moduleLabelsA1)
pdf(file="03.auto_modules_color.pdf",width=12,height=4)
plotDendroAndColors(netA1$dendrograms[[1]],
                    moduleColorsA1[netA1$blockGenes[[1]]],"Module colors",
                    dendroLabels=FALSE,addGuide=TRUE)
dev.off()

save(datExprA1q,sft,softPower1,networkType,netA1,moduleColorsA1,moduleLabelsA1,
     file="WGCNA_2_individualR04.RData")


##Yafei output##
mods <- moduleEigengenes(datExprA1q, colors = moduleColorsA1)

table(mods$validColors)
save(file = "WGCNA_mod_R04_individual_14059genes_power7_signed_20190320.Rdata",mods)
##adjacency,TOM

## print as tables##
for (i in unique(mods$validColors)){
  write.table(file=paste("mod_",i,".txt",sep=""),colnames(datExprA1q)[which(mods$validColors == i)],quote = F,col.names = F,row.names = F)
}


###################################################################################################
# 05. Module eigengenes Calculation

MEs0 <- moduleEigengenes(datExprA1q[,netA1$goodGenes],
                         moduleColorsA1[netA1$blockGenes[[1]]])$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames(datExprA1q[,netA1$goodGenes])
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
  colnames(ME) <- rownames(datExprA1q[,netA1$goodGenes])
  layout(matrix(c(1,2)),heights=c(1.5,3))
  par(mar=c(0.3,9,3,5))
  plotMat(t(scale(datExprA1q[,netA1$goodGenes][,moduleColorsA1[netA1$blockGenes[[1]]]==which.module])),
          nrgcols=30,rlabels=F,rcols=which.module,
          main=paste(which.module),cex.main=1)
  ###colors=blueWhiteRed(50),
  ###color=colorRampPalette(c("blue","white","red"))
  par(mar=c(5,4,0,1))
  barplot(ME,col=which.module,main="",cex.names=1,cex.axis=1,
          ylab="module eigengene",las=3)
  dev.off()
}
#####Jeremy Miller tutorial
modulesA1=mods$validColors
PCs1A = moduleEigengenes(datExprA1q, colors = modulesA1)
#wrong PCs1A = moduleEigengenes(datExprA1q, colors=moduleColorsA1)
#Error in moduleEigengenes(datExprA1q, colors = modulesA1) : 
#object 'modulesA1' not found
ME_1A = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a")
MDS_1A = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesA1))
#Error in moduleEigengenes(datExprA1q, colors = modulesA1) : 
#object 'modulesA1' not found
save.image("tutorial_R04individual.RData")
pdf("ModuleEigengeneVisualizations_R04individual.pdf",height=6,width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1, main="MDS plot", cex=2, pch=19)

ordergenes = netA1$dendrograms[[1]]$order
#Error: object 'geneTreeA1' not found
datExprA1g=datExprA1p
plotMat(scale(log(datExprA1g[ordergenes,])) , rlabels= modulesA1[ordergenes], clabels=
          colnames(datExprA1g), rcols=moduleColorsA1[ordergenes])
#wrong plotMat(scale(log(datExprA1g[ordergenes,])) , rlabels= moduleColorsA1[ordergenes], clabels=
 #          colnames(datExprA1g), rcols=moduleColorsA1[ordergenes])
#Error in plot.mat(scale(log(datExprA1g[ordergenes, ])), rlabels = modulesA1[ordergenes],  : 
#                    could not find function "plot.mat"
for (which.module in names(table(moduleColorsA1))){
  ME = ME_1A[, paste("ME",which.module, sep="")]
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
};
dev.off();



###################################################################################################
# 06. Relating modules to external information (samples/traits)

sample_cor <- cor(t(datExprA1q[,netA1$goodGenes]),t(datExprA1q[,netA1$goodGenes]),
                  use='pairwise.complete.obs')
moduleSampleCor <- cor(MEs,sample_cor,use="p")
nSamples <- nrow(datExprA1q[,netA1$goodGenes])
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

moduleTraitCor <- cor(MEs,datTraits,use="p")
nSamples <- nrow(datExprA1q[,netA1$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste(signif(moduleTraitCor,2),
                    "\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- paste("ME",names(MEs),sep="")
pdf(file="06.modules_traits_relationships.pdf",15,8)
par(mar=c(9,8,1,2))
labeledHeatmap(Matrix=moduleTraitCor,
               textMatrix=textMatrix,
               xLabels=colnames(datTraits),
               yLabels=rowLabels,ySymbols=names(MEs),
               colorLabels=TRUE,colors=blueWhiteRed(50),
               setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
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
geneModuleMembership <- cor(datExprA1q[,netA1$goodGenes],MEs,use="p")
nSamples <- nrow(datExprA1q[,netA1$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership,nSamples)
colnames(geneModuleMembership) <- paste("MM",modNames,sep="")
colnames(MMPvalue) <- paste("p.MM",modNames,sep="")

text <- paste("cor=",round(geneModuleMembership,4),
              ";p-value=",round(MMPvalue,4),sep="")
dim(text) <- dim(geneModuleMembership)
rownames(text) <- rownames(geneModuleMembership)
colnames(text) <- colnames(geneModuleMembership)
text <- cbind(rownames(text),moduleColorsA1[netA1$blockGenes[[1]]],text)
colnames(text)[1] <- "modules"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_module_membership.xls",
            quote=F,sep="\t",row.names=F)







modNames <- names(MEs)
VA_status <- as.data.frame(datTraits$VA_status)
names(VA_status) <- "VA_status"
geneTraitSignificance <- as.data.frame(cor(datExprA1q[,netA1$goodGenes],VA_status,use="p"))
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
  moduleGenes <- moduleColorsA1[netA1$blockGenes[[1]]]==which.module
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
text <- cbind(rownames(text),moduleColorsA1[netA1$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-VAstatus.xls",
            quote=F,sep="\t",row.names=F)

save(datExprA1q,sft,softPower1,networkType,netA1,moduleColorsA1,MEs,
     moduleSampleCor,moduleSamplePvalue,
     moduleTraitCor,moduleTraitPvalue,
     geneModuleMembership,MMPvalue,
     file="WGCNA_2_individualR04.RData")



#####################################################3
############### #####################7.BW_change_percent
modNames <- names(MEs)
BW_change_percent <- as.data.frame(datTraits$BW_change_percent)
names(BW_change_percent) <- "BW_change_percent"
geneTraitSignificance <- as.data.frame(cor(datExprA1q[,netA1$goodGenes],BW_change_percent,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(BW_change_percent),sep="")
names(GSPvalue) <- paste("p.GS.",names(BW_change_percent),sep="")
head(geneTraitSignificance)
##               GS.BW_change_percent
## MMT00000044  0.01253338

head(GSPvalue)
##             p.GS.BW_change_percent
## MMT00000044 0.885713080


dir.create("07.MM_vs_BW_change_percent")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColorsA1[netA1$blockGenes[[1]]]==which.module
  dir <- "07.MM_vs_BW_change_percent/"
  pdf(file=paste(dir,"07.",which.module,"_MM_vs_BW_change_percent.pdf",sep=""),6,6)
  verboseScatterplot(geneModuleMembership[moduleGenes,column],
                     geneTraitSignificance[moduleGenes,1],
                     xlab=paste("Module membership (MM) in",which.module,"module"),
                     ylab="Gene significance for BW_change_percent",
                     main=paste("Module membership vs. gene significance\n"),
                     col=which.module)
  dev.off()
}

text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),moduleColorsA1[netA1$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-BW_change_percent.xls",
            quote=F,sep="\t",row.names=F)

##########################################
#########################################3



#####################################################3
############### #####################7.BW
modNames <- names(MEs)
BW <- as.data.frame(datTraits$BW)
names(BW) <- "BW"
geneTraitSignificance <- as.data.frame(cor(datExprA1q[,netA1$goodGenes],BW,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(BW),sep="")
names(GSPvalue) <- paste("p.GS.",names(BW),sep="")
head(geneTraitSignificance)
##               GS.BW_change_percent
## MMT00000044  0.01253338

head(GSPvalue)
##             p.GS.BW_change_percent
## MMT00000044 0.885713080


dir.create("07.MM_vs_BW")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColorsA1[netA1$blockGenes[[1]]]==which.module
  dir <- "07.MM_vs_BW/"
  pdf(file=paste(dir,"07.",which.module,"_MM_vs_BW.pdf",sep=""),6,6)
  verboseScatterplot(geneModuleMembership[moduleGenes,column],
                     geneTraitSignificance[moduleGenes,1],
                     xlab=paste("Module membership (MM) in",which.module,"module"),
                     ylab="Gene significance for BW",
                     main=paste("Module membership vs. gene significance\n"),
                     col=which.module)
  dev.off()
}

text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),moduleColorsA1[netA1$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-BW.xls",
            quote=F,sep="\t",row.names=F)

##########################################
#########################################3


#####################################################3
############### #####################7.Gender
modNames <- names(MEs)
Gender <- as.data.frame(datTraits$Gender)
names(Gender) <- "Gender"
geneTraitSignificance <- as.data.frame(cor(datExprA1q[,netA1$goodGenes],Gender,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(Gender),sep="")
names(GSPvalue) <- paste("p.GS.",names(Gender),sep="")
head(geneTraitSignificance)
##               GS.BW_change_percent
## MMT00000044  0.01253338

head(GSPvalue)
##             p.GS.BW_change_percent
## MMT00000044 0.885713080


dir.create("07.MM_vs_Gender")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColorsA1[netA1$blockGenes[[1]]]==which.module
  dir <- "07.MM_vs_Gender/"
  pdf(file=paste(dir,"07.",which.module,"_MM_vs_Gender.pdf",sep=""),6,6)
  verboseScatterplot(geneModuleMembership[moduleGenes,column],
                     geneTraitSignificance[moduleGenes,1],
                     xlab=paste("Module membership (MM) in",which.module,"module"),
                     ylab="Gene significance for Gender",
                     main=paste("Module membership vs. gene significance\n"),
                     col=which.module)
  dev.off()
}

text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),moduleColorsA1[netA1$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-Gender.xls",
            quote=F,sep="\t",row.names=F)

##########################################
#########################################3

###################################################################################################
# 08. Exporting network
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller")
load(file="WGCNA_TOM_A1-block.1.RData")
ATOM <- as.matrix(TOM)
TOM1 <- ATOM[1:round((nrow(ATOM)/2)),1:round((nrow(ATOM)/2))]
TOM2 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),1:round((nrow(ATOM)/2))]
TOM3 <- ATOM[1:round((nrow(ATOM)/2)),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
TOM4 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),(round(nrow(ATOM)/2)+1):nrow(ATOM)]

setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190320_individual_R04")
load(file="WGCNA_2_individualR04.RData")

dir.create("08.module_result")
setwd("08.module_result")
for(i in 1:(ncol(MEs)-1)) {
  module <- labels2colors(i)
  inModule <- moduleColorsA1[netA1$blockGenes[[1]]]==module
  genename <- colnames(datExprA1q[,netA1$goodGenes])
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
  IMConn <- softConnectivity(datExprA1q[,netA1$goodGenes][,modGenes])
  ####weight=TOM value
  cyt1 <- exportNetworkToCytoscape(modTOM,
                                   edgeFile=paste("CytoscapeInput-edges-",paste(module,collapse="-"),".txt",sep=""),
                                   nodeFile=paste("CytoscapeInput-nodes-",paste(module,collapse="-"),".txt",sep=""),
                                   weighted=TRUE,threshold=0.02,
                                   nodeNames=modGenes,altNodeNames=modGenes,
                                   nodeAttr=moduleColorsA1[netA1$blockGenes[[1]]][inModule])
  
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







genename <- colnames(datExprA1q[,netA1$goodGenes])
modGenes <- genename
IMConn <- softConnectivity(datExprA1q[,netA1$goodGenes][,modGenes])
cyt1 <- exportNetworkToCytoscape(modTOM,
                                 edgeFile=paste("CytoscapeInput-edges-overall.txt",sep=""),
                                 nodeFile=paste("CytoscapeInput-nodes-overall.txt",sep=""),
                                 weighted=TRUE,threshold=0.02,
                                 nodeNames=modGenes,altNodeNames=modGenes,
                                 nodeAttr=moduleColorsA1[netA1$blockGenes[[1]]])
#> cyt1 <- exportNetworkToCytoscape(ATOM,
#Error: cannot allocate vector of size 754.0 Mb
out <- cbind(modGenes,IMConn)
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste("overall_gene_connectivity.txt",sep=""),
            sep="\t",quote=F,row.names=F)





##intramodularConnectivity.=connectivity of nodes to other nodes within the same module.
IntraMConn <- intramodularConnectivity.fromExpr(datExprA1q[,netA1$goodGenes],
                                                moduleColorsA1[netA1$blockGenes[[1]]],
                                                networkType=networkType,
                                                power=softPower1)
merged_IntraMConn <- cbind(colnames(datExprA1q[,netA1$goodGenes]),
                           moduleColorsA1[netA1$blockGenes[[1]]],IntraMConn)
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



























#################individual modules-  A2, R11

#networkType <- "signed"
#netA2 <- blockwiseModules(t(datExprA2p),power=softPower2,
#                          networkType=networkType,numericLabels=TRUE,
#                          mergeCutHeight=0.25,minModuleSize=30,
#                          maxBlockSize=30000,saveTOMs=TRUE,
#                          saveTOMFileBase="WGCNA_TOM_A2",verbose=5)
#moduleLabelsA2 <- netA2$colors
#moduleColorsA2 <- labels2colors(moduleLabelsA2)


setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190320_individual_R11")

gsg <- goodSamplesGenes(t(datExprA2p))
## Flagging genes and samples with too many missing values...（expression level NA/zero in more thant 50% samples）
##  ..step 1
gsg$allOK
## [1] TRUE

sampleTree <- hclust(dist(t(datExprA2p)),method="average")
pdf("01.samples_cluster_tree.pdf",width=25,height=8)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="")
#cutHeight <- 15
abline(h=cutHeight,col="red")
dev.off()


datExprA2q <- t(datExprA2p)

traitData = read.csv("20190312-ClinicalTrait_R11.csv");
dim(traitData)
names(traitData)
head(traitData)

rownames(datTraits) = traitData[traitRows, 1];
#Error in `[.data.frame`(traitData, traitRows, 1) : 
#  object 'traitRows' not found
MiceSamples = rownames(datExprA2q);
traitRows = match(MiceSamples, traitData$Mice);
datTraits = traitData[traitRows, -c(1,2)];

collectGarbage();


sampleTree2 <- hclust(dist(datExprA2q),method="average")
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
sft <- pickSoftThreshold(datExprA2q,
                         powerVector=powers,
                         networkType="signed",
                         verbose=5)
# networkType = "unsigned", adjacency = |cor|^power;
#signed可以区分极大正相关和极大负相关
# networkType = "signed", adjacency = ((1+cor)/2)^power;
# networkType = "signed hybrid", adjacency = cor^power if cor>0 and 0 otherwise;
# networkType = "distance", adjacency = (1-(dist/max(dist))^2)^power.

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
# 03. One-step network construction and module detection--after

#moduleLabelsA1 <- netA1$colors
#moduleColorsA1 <- labels2colors(moduleLabelsA1)
pdf(file="03.auto_modules_color.pdf",width=12,height=4)
plotDendroAndColors(netA2$dendrograms[[1]],
                    moduleColorsA2[netA2$blockGenes[[1]]],"Module colors",
                    dendroLabels=FALSE,addGuide=TRUE)
dev.off()

save(datExprA2q,sft,softPower2,networkType,netA2,moduleColorsA2,moduleLabelsA2,
     file="WGCNA_2_individualR11.RData")


##Yafei output##
mods <- moduleEigengenes(datExprA2q, colors = moduleColorsA2)

table(mods$validColors)
save(file = "WGCNA_mod_R11_individual_14059genes_power26_signed_20190321.Rdata",mods)
##adjacency,TOM

## print as tables##
for (i in unique(mods$validColors)){
  write.table(file=paste("mod_",i,".txt",sep=""),colnames(datExprA2q)[which(mods$validColors == i)],quote = F,col.names = F,row.names = F)
}


###################################################################################################
# 05. Module eigengenes Calculation

MEs0 <- moduleEigengenes(datExprA2q[,netA2$goodGenes],
                         moduleColorsA2[netA2$blockGenes[[1]]])$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames(datExprA2q[,netA2$goodGenes])
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
  colnames(ME) <- rownames(datExprA2q[,netA2$goodGenes])
  layout(matrix(c(1,2)),heights=c(1.5,3))
  par(mar=c(0.3,9,3,5))
  plotMat(t(scale(datExprA2q[,netA2$goodGenes][,moduleColorsA2[netA2$blockGenes[[1]]]==which.module])),
          nrgcols=30,rlabels=F,rcols=which.module,
          main=paste(which.module),cex.main=1)
  ###colors=blueWhiteRed(50),
  ###color=colorRampPalette(c("blue","white","red"))
  par(mar=c(5,4,0,1))
  barplot(ME,col=which.module,main="",cex.names=1,cex.axis=1,
          ylab="module eigengene",las=3)
  dev.off()
}
#####Jeremy Miller tutorial
modulesA1=mods$validColors
PCs1A = moduleEigengenes(datExprA2q, colors = modulesA1)
#wrong PCs1A = moduleEigengenes(datExprA1q, colors=moduleColorsA1)
#Error in moduleEigengenes(datExprA1q, colors = modulesA1) : 
#object 'modulesA1' not found
ME_1A = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a")
MDS_1A = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesA1))
#Error in moduleEigengenes(datExprA1q, colors = modulesA1) : 
#object 'modulesA1' not found
save.image("tutorial_R11individual.RData")
pdf("ModuleEigengeneVisualizations_R11individual.pdf",height=6,width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1, main="MDS plot", cex=2, pch=19)

ordergenes = netA2$dendrograms[[1]]$order
#Error: object 'geneTreeA1' not found
datExprA2g=datExprA2p
plotMat(scale(log(datExprA2g[ordergenes,])) , rlabels= modulesA1[ordergenes], clabels=
          colnames(datExprA2g), rcols=moduleColorsA2[ordergenes])
#wrong plotMat(scale(log(datExprA1g[ordergenes,])) , rlabels= moduleColorsA1[ordergenes], clabels=
#          colnames(datExprA1g), rcols=moduleColorsA1[ordergenes])
#Error in plot.mat(scale(log(datExprA1g[ordergenes, ])), rlabels = modulesA1[ordergenes],  : 
#                    could not find function "plot.mat"
for (which.module in names(table(moduleColorsA2))){
  ME = ME_1A[, paste("ME",which.module, sep="")]
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
};
dev.off();



###################################################################################################
# 06. Relating modules to external information (samples/traits)

sample_cor <- cor(t(datExprA2q[,netA2$goodGenes]),t(datExprA2q[,netA2$goodGenes]),
                  use='pairwise.complete.obs')
moduleSampleCor <- cor(MEs,sample_cor,use="p")
nSamples <- nrow(datExprA2q[,netA2$goodGenes])
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

moduleTraitCor <- cor(MEs,datTraits,use="p")
nSamples <- nrow(datExprA2q[,netA2$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste(signif(moduleTraitCor,2),
                    "\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- paste("ME",names(MEs),sep="")
pdf(file="06.modules_traits_relationships.pdf",15,8)
par(mar=c(9,8,1,2))
labeledHeatmap(Matrix=moduleTraitCor,
               textMatrix=textMatrix,
               xLabels=colnames(datTraits),
               yLabels=rowLabels,ySymbols=names(MEs),
               colorLabels=TRUE,colors=blueWhiteRed(50),
               setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
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
geneModuleMembership <- cor(datExprA2q[,netA2$goodGenes],MEs,use="p")
nSamples <- nrow(datExprA2q[,netA2$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership,nSamples)
colnames(geneModuleMembership) <- paste("MM",modNames,sep="")
colnames(MMPvalue) <- paste("p.MM",modNames,sep="")

text <- paste("cor=",round(geneModuleMembership,4),
              ";p-value=",round(MMPvalue,4),sep="")
dim(text) <- dim(geneModuleMembership)
rownames(text) <- rownames(geneModuleMembership)
colnames(text) <- colnames(geneModuleMembership)
text <- cbind(rownames(text),moduleColorsA2[netA2$blockGenes[[1]]],text)
colnames(text)[1] <- "modules"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_module_membership.xls",
            quote=F,sep="\t",row.names=F)







modNames <- names(MEs)
VA_status <- as.data.frame(datTraits$VA_status)
names(VA_status) <- "VA_status"
geneTraitSignificance <- as.data.frame(cor(datExprA2q[,netA2$goodGenes],VA_status,use="p"))
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
  moduleGenes <- moduleColorsA2[netA2$blockGenes[[1]]]==which.module
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
text <- cbind(rownames(text),moduleColorsA2[netA2$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-VAstatus.xls",
            quote=F,sep="\t",row.names=F)

save(datExprA2q,sft,softPower2,networkType,netA2,moduleColorsA2,moduleLabelsA2,MEs,
     moduleSampleCor,moduleSamplePvalue,
     moduleTraitCor,moduleTraitPvalue,
     geneModuleMembership,MMPvalue,
     file="WGCNA_2_individualR11.RData")



#####################################################3
############### #####################7.Infection_status
modNames <- names(MEs)
Infection_status <- as.data.frame(datTraits$Infection_status)
names(Infection_status) <- "Infection_status"
geneTraitSignificance <- as.data.frame(cor(datExprA2q[,netA2$goodGenes],Infection_status,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(Infection_status),sep="")
names(GSPvalue) <- paste("p.GS.",names(Infection_status),sep="")
head(geneTraitSignificance)
##               GS.BW_change_percent
## MMT00000044  0.01253338

head(GSPvalue)
##             p.GS.BW_change_percent
## MMT00000044 0.885713080


dir.create("07.MM_vs_Infection_status")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColorsA2[netA2$blockGenes[[1]]]==which.module
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
text <- cbind(rownames(text),moduleColorsA2[netA2$blockGenes[[1]]],text)
#text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
colnames(text)[2] <- "modulesColors"
write.table(text,file="07.genes_trait_significance-Infection_status.xls",
            quote=F,sep="\t",row.names=F)

##########################################
#########################################3









###################################################################################################
# 08. Exporting network
setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller")
load(file="WGCNA_TOM_A2-block.1.RData")
ATOM <- as.matrix(TOM)
TOM1 <- ATOM[1:round((nrow(ATOM)/2)),1:round((nrow(ATOM)/2))]
TOM2 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),1:round((nrow(ATOM)/2))]
TOM3 <- ATOM[1:round((nrow(ATOM)/2)),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
TOM4 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),(round(nrow(ATOM)/2)+1):nrow(ATOM)]

setwd("C:/R_work/20190320 module preservation softpower 7 n 26-helix lesson5-Jeremy Miller/20190320_individual_R11")
load(file="WGCNA_2_individualR11.RData")

dir.create("08.module_result")
setwd("08.module_result")
for(i in 1:(ncol(MEs)-1)) {
  module <- labels2colors(i)
  inModule <- moduleColorsA2[netA2$blockGenes[[1]]]==module
  genename <- colnames(datExprA2q[,netA2$goodGenes])
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
  IMConn <- softConnectivity(datExprA2q[,netA2$goodGenes][,modGenes])
  ####weight=TOM value
  cyt1 <- exportNetworkToCytoscape(modTOM,
                                   edgeFile=paste("CytoscapeInput-edges-",paste(module,collapse="-"),".txt",sep=""),
                                   nodeFile=paste("CytoscapeInput-nodes-",paste(module,collapse="-"),".txt",sep=""),
                                   weighted=TRUE,threshold=0.02,
                                   nodeNames=modGenes,altNodeNames=modGenes,
                                   nodeAttr=moduleColorsA2[netA2$blockGenes[[1]]][inModule])
  
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







genename <- colnames(datExprA2q[,netA2$goodGenes])
modGenes <- genename
IMConn <- softConnectivity(datExprA2q[,netA2$goodGenes][,modGenes])
cyt1 <- exportNetworkToCytoscape(modTOM,
                                 edgeFile=paste("CytoscapeInput-edges-overall.txt",sep=""),
                                 nodeFile=paste("CytoscapeInput-nodes-overall.txt",sep=""),
                                 weighted=TRUE,threshold=0.02,
                                 nodeNames=modGenes,altNodeNames=modGenes,
                                 nodeAttr=moduleColorsA2[netA2$blockGenes[[1]]])
#> cyt1 <- exportNetworkToCytoscape(ATOM,
#Error: cannot allocate vector of size 754.0 Mb
out <- cbind(modGenes,IMConn)
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste("overall_gene_connectivity.txt",sep=""),
            sep="\t",quote=F,row.names=F)





##intramodularConnectivity.=connectivity of nodes to other nodes within the same module.
IntraMConn <- intramodularConnectivity.fromExpr(datExprA2q[,netA2$goodGenes],
                                                moduleColorsA2[netA2$blockGenes[[1]]],
                                                networkType=networkType,
                                                power=softPower2)
merged_IntraMConn <- cbind(colnames(datExprA2q[,netA2$goodGenes]),
                           moduleColorsA2[netA2$blockGenes[[1]]],IntraMConn)
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
