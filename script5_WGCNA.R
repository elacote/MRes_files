#BiocManager::install("tximportData")
library(tximportData)
library(tximport)

#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") 
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

#BiocManager::install("GenomicFeatures")
library("GenomicFeatures")

#BiocManager::install("DESeq2")
library("DESeq2")

#BiocManager::install("ggplot2")
library("ggplot2")

#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

#BiocManager::install("AnnotationDbi")
library("AnnotationDbi")

#BiocManager::install("clusterProfiler")
library("clusterProfiler")
library(enrichplot)

#install.packages("readxl")
library(readxl)

#BiocManager::install("WGCNA")
library("WGCNA")

options(stringsAsFactors = FALSE)

#BiocManager::install("flashClust")
library("flashClust")

#WGCNA on NPCs using dataset 2 ASD samples from iPSC- adjusted_counts_NPC_griesi variable from DE analysis
#counts normalised with Combat_seq

#matrix with information on the experiment
coldata2 <- readRDS("~/Desktop/MRes/20220412/coldata2.rds")
NPC_coldata <-subset(coldata2, subset = cell_type =="NPC")
NPC_coldata$cell_type <-NULL
NPC_coldata$condition <-c(rep(1,16), rep(0,13)) #ASD: 1; CTL:0
NPC_coldata$group <-NULL
datTraits_NPC_griesi <- NPC_coldata

#construct a DESeqDataSet NPC
dds_NPC <- DESeqDataSetFromMatrix(countData = adjusted_counts_NPC_griesi,
                                  colData = datTraits_NPC_griesi,
                                  design = ~1)

vsd_NPC <-vst(dds_NPC, blind=FALSE)
mat <- assay(vsd_NPC)

#check data
mean_expr <-apply(mat,2, mean, na.rm =T)
number_missing <-apply(is.na(data.frame(mat)),2,sum)

sizeGrWindow(12, 8)
barplot(mean_expr,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:29))
number_missing

NumberMissingByGene =apply( is.na(data.frame(mat)),1, sum)
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(mat),1,var, na.rm=T))

no.presentdatExpr=as.vector(apply(!is.na(as.matrix(mat)),1, sum) )
table(no.presentdatExpr)

# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)
mat=mat[KeepGenes, ]

matNPC_griesi <-t(mat)

#same analysis for dataset 1_ IFN_ iPSC
dir <- "~/Desktop/MRes/20210402/IFN_from_kallisto/"
samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)
#list of files to be used to get abundance data for all the samples
files <- file.path(dir,  samples$V1, "abundance.tsv")
names(files) <- paste0("sample", 1:18)
all(file.exists(files))
#create a df with a correspondance between transcripts and genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
#select NPCs only
txi_NPC_IFN <- lapply(txi, function(x) if(is.matrix(x)) return(x[,1:6]) else return(x))

#import matrix with metadata
IFN_coldata <- readRDS("~/Desktop/MRes/20210214/coldata1.rds")
NPC_coldata_IFN <-subset(IFN_coldata, subset = cell_type =="NPC")
NPC_coldata_IFN$D30 <-NULL
NPC_coldata_IFN$cell_type <-NULL
NPC_coldata_IFN$cell_line <- c(rep(c(1,2,3),2))
NPC_coldata_IFN$treatment <- c(rep(1,3), rep(0,3))
NPC_coldata_IFN$D18 <- NULL
datTraits_NPC_IFN <- NPC_coldata_IFN

colnames(txi_NPC_IFN$counts) <-rownames(datTraits_NPC_IFN)

dds <- DESeqDataSetFromTximport(txi_NPC_IFN,
                                colData = datTraits_NPC_IFN,
                                design = ~ 1)

#prefiltering: only keep genes with at least 10 reads aligned
#Gene removal with less than 10 reads on all the conditions
keep <- rowSums(counts(dds)) >= 10
table(keep)
#FALSE  TRUE 
# 10081 20565    
dds <- dds[keep,]
vsd <-vst(dds, blind=FALSE)
mat_NPCs_D <- assay(vsd)

#check data
mean_expr <-apply(mat_NPCs_D,2, mean, na.rm =T)
number_missing <-apply(is.na(data.frame(mat_NPCs_D)),2,sum)

sizeGrWindow(12, 8)
barplot(mean_expr,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:6))
number_missing

NumberMissingByGene =apply( is.na(data.frame(mat_NPCs_D)),1, sum)
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(mat_NPCs_D),1,var, na.rm=T))

no.presentdatExpr=as.vector(apply(!is.na(as.matrix(mat_NPCs_D)),1, sum) )
table(no.presentdatExpr)

# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)
mat_NPCs_D=mat_NPCs_D[KeepGenes, ]

matNPC_IFN <-t(mat_NPCs_D)

#ensure the variables are comparable
#limit the analysis to genes that are expressed in both datasets
commonGenes <-intersect(rownames(mat), rownames(mat_NPCs_D))
length(commonGenes)  #17746
mat2NPC_griesi <-mat[commonGenes,]
mat2NPC_deepak<-mat_NPCs_D[commonGenes,]

#find appropriate soft power IFN dataset
powers <- c(c(1:20), seq(from=22, to=30, by=2))
sft <-pickSoftThreshold(t(mat2NPC_deepak), powerVector = powers, verbose =5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 <- 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.9, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <-22

pdf("soft_power_IFN.pdf",height=8,width=12)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.9, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#find appropriate soft power ASD dataset
powers <- c(c(1:20), seq(from=22, to=30, by=2))
sft <-pickSoftThreshold(t(mat2NPC_griesi), powerVector = powers, verbose =5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 <- 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.9, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <-13

pdf("soft_power_ASD.pdf",height=8,width=12)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.8, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#correlate measures of average gene expression and overall connectivity between 2 datasets
#the higher the correaltion of these properties, the better chance of finding similarities
#random sample of 5000 genes in the 10000 genes with highest variance in griesi's dataset
variancedatExpr=as.vector(apply(as.matrix(mat2NPC_griesi),1,var, na.rm=T))
variancedatExpr_order <-as.data.frame(variancedatExpr)
variancedatExpr_order$gene <-rownames(mat2NPC_griesi)
variancedatExpr_order <-variancedatExpr_order[order(variancedatExpr_order$variancedatExpr, decreasing = TRUE),]
variancedatExpr_order_10000 <-variancedatExpr_order[1:10000,]
random5000<- sample(variancedatExpr_order_10000$variancedatExpr,5000)

rankExpr_griesi<- rank(rowMeans(mat2NPC_griesi))
rankExpr_deepak<- rank(rowMeans(mat2NPC_deepak))
rankConn_griesi<- rank(softConnectivity(t(mat2NPC_griesi[random5000,]),type="signed",power=13)) 
rankConn_deepak<- rank(softConnectivity(t(mat2NPC_deepak[random5000,]),type="signed",power=22))

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
pdf("generalNetworkProperties.pdf", height=10, width=9)
par(mfrow=c(2,2))
verboseScatterplot(rankExpr_griesi,rankExpr_deepak, xlab="Ranked Expression (ASD)",ylab="Ranked Expression (IFN)") 
verboseScatterplot(rankConn_griesi,rankConn_deepak, xlab="Ranked Connectivity (ASD)", ylab="Ranked Connectivity (IFN)")
dev.off()

#run WGCNA on the datasets
adjacency_deepak <- adjacency(t(mat2NPC_deepak),power=22,type="signed") 
diag(adjacency_deepak)<-0
dissTOM_deepak <- 1-TOMsimilarity(adjacency_deepak, TOMType="signed") 
geneTree_deepak <- flashClust(as.dist(dissTOM_deepak), method="average")

adjacency_griesi <- adjacency(t(mat2NPC_griesi),power=13,type="signed") 
diag(adjacency_griesi)=0
TOM_griesi_NPC <-TOMsimilarity(adjacency_griesi, TOMType="signed")
dissTOM_griesi <- 1-TOMsimilarity(adjacency_griesi, TOMType="signed") 
geneTree_griesi <- flashClust(as.dist(dissTOM_griesi), method="average")

pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree_griesi,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (ASD)", labels=FALSE,hang=0.04);
plot(geneTree_deepak,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (IFN)", labels=FALSE,hang=0.04);
dev.off()
#"good" data when a lot of distict branches

#determine modules based on griesi's dataset
mColorh <-NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTree_griesi, pamStage=FALSE, cutHeight = 0.99, minClusterSize = (150), deepSplit = ds, distM = dissTOM_griesi)
  mColorh=cbind(mColorh,labels2colors(tree$labels)) }
pdf("Module_choices2.pdf", height=10,width=25)
plotDendroAndColors(geneTree_griesi, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE)
dev.off()
modules_griesi1 = mColorh[,1] #we choose ds=0

#calculate the principle components for visualisation and quantitative assessemnt
#first PC= module eigengene (ME) = single value that represents the highest percent of variance for all genes in a module (1 par module et par sample)
#most genes in the module "do the same thing" as ME
PCs1_deepak <- moduleEigengenes(t(mat2NPC_deepak), colors=modules_griesi1)

ME_1_deepak <-PCs1_deepak$eigengenes #table sample/module
distPC1_deepak <-1-abs(cor(ME_1_deepak,use="p")) #matrice module-module
distPC1_deepak <-ifelse(is.na(distPC1_deepak), 0, distPC1_deepak)
pcTree1_deepak <-hclust(as.dist(distPC1_deepak),method="a")
MDS_1_deepak <-cmdscale(as.dist(distPC1_deepak),2)
colors_deepak1 <-names(table(modules_griesi1))

PCs1_griesi = moduleEigengenes(t(mat2NPC_griesi), colors=modules_griesi1)
ME_1_griesi = PCs1_griesi$eigengenes
colors_griesi1 <-names(table(modules_griesi1))

pdf("ModuleEigengeneVisualizations.pdf",height=6,width=6) 
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)

plot(pcTree1_deepak, xlab="",ylab="",main="",sub="")
plot(MDS_1_deepak, col= colors_deepak1, main="MDS plot", cex=2, pch=19)
ordergenes = geneTree_deepak$order
for (which.module in names(table(modules_griesi1))){
  ME = ME_1_deepak[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, ylab="eigengene expression",xlab=" sample") };
dev.off();

#associate with traits
#treatment:1 ; no treatment: 0
moduleTraitCor = cor(ME_1_deepak, datTraits_NPC_IFN , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 6)

sizeGrWindow(12,6)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_NPC_IFN),
               yLabels = names(ME_1_deepak),
               ySymbols = names(ME_1_deepak),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#tan module associated with treatment

#for griesi's dataset
ME_1_griesi = orderMEs(ME_1_griesi)
moduleTraitCor = cor(ME_1_griesi, datTraits_NPC_griesi, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 29)

sizeGrWindow(12,6)
#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_NPC_griesi),
               yLabels = names(ME_1_griesi),
               ySymbols = names(ME_1_griesi),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#assess how well modules from ASD dataset are preserved in IFN dataset
pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTree_deepak, modules_griesi1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (IFN-gamma)") 
plotDendroAndColors(geneTree_griesi, modules_griesi1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (ASD)") 

dev.off()

#how well a model in a study is preserved in another one (outputs a single Z-score summary)
multiExpr <- list(griesi=list(data=t(mat2NPC_griesi)), deepak=list(data=t(mat2NPC_deepak)))
multiColor <- list(griesi = modules_griesi1) 
mp<-modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                       nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000) 
stats <- mp$preservation$Z$ref.griesi$inColumnsAlsoPresentIn.deepak 
stats_NPC <-stats[order(-stats[,2]),c(1:2)]
View(stats_NPC)

#5<Z<10: moderate preservation- Z>10: high preservation
#grey: uncharacterized genes; gold: random genes

#module membership=kME
#value used to measure correlations between each gene and each ME
# => genes not initially assigned to a module can be included in between-network comparisons
#get kME values and associated p-values
geneModuleMembership1 <- signedKME(t(mat2NPC_deepak), ME_1_deepak) 
colnames(geneModuleMembership1)<-paste("PC",colors_griesi1,".cor",sep="")

MMPvalue1<- corPvalueStudent(as.matrix(geneModuleMembership1),dim(mat2NPC_deepak)[[2]]) 
colnames(MMPvalue1)<- paste("PC",colors_griesi1,".pval",sep="")

Gene_NPC_griesi <- rownames(mat2NPC_griesi) 
kMEtable1 <- cbind(Gene_NPC_griesi,Gene_NPC_griesi,modules_griesi1) 
for (i in 1:length(colors_griesi1))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"kMEtable1deepakNPCintersection.csv",row.names=FALSE)

#kME values and p-values for griesi
geneModuleMembership2 <- signedKME(t(mat2NPC_griesi), ME_1_griesi ) 
colnames(geneModuleMembership2)<-paste("PC",colors_griesi1,".cor",sep="")

MMPvalue2<- corPvalueStudent(as.matrix(geneModuleMembership2),dim(mat2NPC_griesi)[[2]])
colnames(MMPvalue2)<- paste("PC",colors_griesi1,".pval",sep="")

kMEtable2 <- cbind(Gene_NPC_griesi,Gene_NPC_griesi,modules_griesi1) 
for (i in 1:length(colors_griesi1))
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i]) 
colnames(kMEtable2)<-colnames(kMEtable1)
write.csv(kMEtable2,"kMEtable2griesiNPCintersection.csv",row.names=FALSE)


#which genes are hubs in both networks (genes with extremely high kME values in both netwroks)
topGenesKME <- NULL
for (c in 1:length(colors_griesi1)){
  kMErank1 = rank(-geneModuleMembership1[,c])
  kMErank2 = rank(-geneModuleMembership2[,c])
  maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max)) 
  topGenesKME = cbind(topGenesKME,Gene_NPC_griesi[maxKMErank<=20])
} 
colnames(topGenesKME) <- colors_griesi1 
View(topGenesKME)

topGeneKME <-as.data.frame(topGenesKME)

module_brown <-data.frame(
  gene_ID <- topGenesKME[,3],
  SYMBOL <-topGenesKME[,3]
)
View(module_brown)
colnames(module_brown) <-c('gene_ID', 'SYMBOL')

# Add gene symbol column
module_brown$SYMBOL = mapIds(org.Hs.eg.db,
                             keys=gene_ID, 
                             column="SYMBOL",
                             keytype="ENTREZID",
                             multiVals="first")

#gene significance= correlation between the gene and the trait
condition <- as.data.frame(datTraits_NPC_griesi$condition)
names(condition) <- "condition" #CTL or ASD
# names (colors) of the modules
modNames <- substring(names(ME_1_griesi), 3)
names(geneModuleMembership2) <- paste("MM", modNames, sep="");
names(MMPvalue2) <- paste("p.MM", modNames, sep="");
geneTraitSignificance <- as.data.frame(cor(t(mat2NPC_griesi), condition, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 29));
names(geneTraitSignificance) <- paste("GS.", names(condition), sep="");
names(GSPvalue) <- paste("p.GS.", names(condition), sep="");


#intramodular analysis
#genes with high GS and high MM
module ="brown"
#module ="red"
#module ="tan"
column <-match(module, modNames)
moduleGenes <- modules_griesi1 ==module

sizeGrWindow(7,7)
par(mfrow =c(1,1))
verboseScatterplot(abs(geneModuleMembership2[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab ="Gene significance for condition",
                   main = paste("module membership vs gene significance\n"),
                   cex.main =1.2, cex.lab = 1.2, cex.axis =1.2, col = module)


geneInfo0_NPC_griesi <- data.frame(ENTREZID = rownames(mat2NPC_griesi),
                                  moduleColor = modules_griesi1,
                                  geneTraitSignificance,
                                  GSPvalue)
View(geneInfo0_NPC_griesi)
geneInfo0_NPC_griesi$SYMBOL <- mapIds(org.Hs.eg.db,
                                     keys=rownames(geneInfo0_NPC_griesi), 
                                     column="SYMBOL",
                                     keytype="ENTREZID",
                                     multiVals="first")

modOrder <- order(-abs(cor(ME_1_griesi, condition, use = "p")))
for (mod in 1:ncol(geneModuleMembership2))
{
  oldNames <- names(geneInfo0_NPC_griesi)
  geneInfo0_NPC_griesi<- data.frame(geneInfo0_NPC_griesi, geneModuleMembership2[, modOrder[mod]], 
                                   MMPvalue2[, modOrder[mod]]);
  names(geneInfo0_NPC_griesi) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                  paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder <- order(geneInfo0_NPC_griesi$moduleColor, -abs(geneInfo0_NPC_griesi$GS.condition))
geneInfo_NPC_griesi <- geneInfo0_NPC_griesi[geneOrder, ]

View(geneInfo_NPC_griesi)

brown <- subset(geneInfo_NPC_griesi, moduleColor == "brown")
View(brown )

brown_orderMMbrown <- brown[order(-abs(brown$MM.brown)),]
View(brown_orderMMbrown)
#save brown module NPC griesi WGCNA
write.csv(brown_orderMMbrown ,"NPC_WGCNA_MEbrown_griesi.csv",row.names=F)
#30 hub genes
hub_genes_30_MEbrown_NPC <-head(brown_orderMMbrown$ENTREZID, 30)
hub_genes_30_MEbrown_NPC_symbol <-head(brown_orderMMbrown$SYMBOL, 30)
writeLines(hub_genes_30_MEbrown_NPC_symbol, "MEbrownNPC_hub30.txt")


gene_studied <-"HLA"
line_gene <-grep(gene_studied,brown_orderMMbrown$SYMBOL)
line_gene

#GO enrichment of genes from brown module signif associated with ASD
# enrichment in term of BP of the set of  genes from interesting module
#NPC_tan <- tan3$ENTREZID
NPC_brown <- brown$ENTREZID
#red1 <-red[order(abs(red$MM.red), decreasing= TRUE),]

#save red module NPC griesi WGCNA
#write.csv(red1,"NPC_WGCNA_MEred_griesi.csv",row.names=F)
#save tan module NPC griesi WGCNA
#write.csv(tan3,"NPC_WGCNA_MEtan_griesi.csv",row.names=F)

#NPC_red <-red1$ENTREZID
#NPC_tan <-tan3$ENTREZID

ego2_NPC <- clusterProfiler::enrichGO(gene          = NPC_brown,
                                      OrgDb         = OrgDb,
                                      ont           = "BP",
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.01, 
                                      readable      = TRUE)

x2 <- pairwise_termsim(ego2_NPC)
emapplot(x2, showCategory = 15, node_label = "category")

writeLines(NPC_brown, "NPC_brown_genes.txt")
#writeLines(NPC_tan, "NPC_tan_genes.txt")
#writeLines(NPC_red, "NPC_red_genes.txt")

#overlap avec DE list IFN
MEbrown_NPC_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = NPC_brown, columns ="SYMBOL", keytype = "ENTREZID")

go.obj <-newGeneOverlap(MEbrown_NPC_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 20068 ) #64, 0.014
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap_MEbrown_DE_neurons_IFN <-getIntersection(go.obj)
overlap_MEbrown_DE_neurons_IFN
#[1] "RIBC1"      "BRINP2"     "RNF213-AS1" "NRN1"       "OPCML"      "C2CD4C"     "ECEL1"      "SIX3"       "SP8"        "ZIC3"       "LINC01833" 
#[12] "SYT13"      "LAMB2"      "DGLUCY"     "SDK2"       "FAM111A-DT" "DACH1"      "SLITRK1"    "RASGRP1"    "LRATD2"     "SPINT2"     "CMPK2"     
#[23] "RAB3B"      "TMEM255B"   "MYO5B"      "FN1"        "JUP"        "SAMD9"      "ZNF503"     "EPHA7"      "PRECSIT"    "LOC643201"  "HS3ST2"    
#[34] "KCNN1"      "UNC5D"      "CA11"       "AXIN2"      "SERPING1"   "DDX60L"     "CDH13"      "SLFN11"     "HLA-L"      "MIR10393"   "SP140L"    
#[45] "SLC12A7"    "REC8"       "MYOF"       "TIMP2"      "EFNB1"      "IFI6"       "LTBP4"      "RNASET2"    "PRMT2"      "LAMA1"      "IFI30"     
#[56] "RELN"       "RAB27B"     "PTGER3"     "TMEM255A"   "TNR"        "EXOC3L1"    "PRPH"       "IFITM1"     "IFI35" 

go.obj <-newGeneOverlap(MEbrown_NPC_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )
go.obj <-testGeneOverlap(go.obj) #113, p-val =0.063
go.obj
overlap_MEbrown_DE_NPC_IFN <-getIntersection(go.obj)
overlap_MEbrown_DE_NPC_IFN
#"DCT"          "PCP4"         "LOC105376095" "CXCL14"       "EFEMP1"       "GCNT1"        "RSPO1"        "CDH20"        "LHX2"        
#[10] "NRN1"         "CHL1"         "FEZF2"        "FGF8"         "SPAG17"       "FBXO44"       "C1S"          "FBXO2"        "ELN"         
#[19] "BLACAT1"      "DGLUCY"       "TNFSF10"      "ICAM5"        "DLX6-AS1"     "FAM111A-DT"   "OTX2-AS1"     "DNAAF3"       "PHYKPL"      
#[28] "LRATD2"       "CMPK2"        "ARX"          "TMEM255B"     "MYO5B"        "ATP1A2"       "FN1"          "KLHDC8B"      "MEG3"        
#[37] "MIR9-1HG"     "SAMD9"        "EPHA7"        "TEC"          "PLEKHM1P1"    "SYT17"        "PRECSIT"      "LOC643201"    "ALPL"        
#[46] "NRXN3"        "LINC02609"    "CFAP47"       "UNC5D"        "CA11"         "SMIM14"       "BMP7"         "NLGN3"        "SERPING1"    
#[55] "DDX60L"       "KCNB1"        "IGSF1"        "CHRNB1"       "LEF1"         "PLEKHF1"      "TSPAN11"      "SLFN11"       "CA12"        
#[64] "HLA-L"        "MIR10393"     "BCAM"         "IQSEC3"       "ARRDC4"       "SCRG1"        "LRP10"        "HSPA12A"      "PAIP2B"      
#[73] "SP140L"       "TLR3"         "CAPS2"        "SLC12A7"      "REC8"         "NTNG2"        "YPEL5"        "CTNNBIP1"     "EFNB1"       
#[82] "PGGHG"        "TKTL1"        "IFI6"         "ZBBX"         "CHRNA4"       "CCDC74B"      "FOXF1"        "PRMT2"        "LAMA1"       
#[91] "IFI30"        "LMO1"         "ARL6IP5"      "SYT7"         "C6orf141"     "KIF19"        "SUSD6"        "BMP4"         "MPP7"        
#[100] "AGTPBP1"      "TMEM255A"     "C2"           "B3GNTL1"      "PRPH"         "CKMT1B"       "IFITM1"       "BRI3"         "DUSP10"      
#[109] "GRIPAP1"      "RPS6KA5"      "IFI35"        "NSMF"         "CETN2" 

go.obj <-newGeneOverlap(list_SFARI, overlap_MEbrown_DE_neurons_IFN, 21581 )
go.obj <-testGeneOverlap(go.obj)
go.obj
getIntersection(go.obj) #"CDH13" "CMPK2" "LAMA1" "RELN" 

#overlap IFN MEred NPC griesi WGCNA
#MEred_NPC_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = NPC_red, columns ="SYMBOL", keytype = "ENTREZID")

#go.obj <-newGeneOverlap(MEred_NPC_genes_symbol$SYMBOL, DE_list_NPC_IFN $SYMBOL, 20068 )
#go.obj <-newGeneOverlap(MEred_NPC_genes_symbol$SYMBOL, DE_list_neurons$SYMBOL, 20068 ) #7, 0.99
#go.obj <-testGeneOverlap(go.obj) #41, p-val =1
#go.obj



#hub genes for the modules
hub_NPC_griesi <-chooseTopHubInEachModule(t(mat2NPC_griesi),  modules_griesi1, omitColors = "grey", power =13 )
hub_NPC_griesi <-as.data.frame(hub_NPC_griesi)
colnames(hub_NPC_griesi)
hub_NPC_griesi$SYMBOL = mapIds(org.Hs.eg.db,
                               keys=hub_NPC_griesi$hub_NPC_griesi, 
                               column="SYMBOL",
                               keytype="ENTREZID",
                               multiVals="first")
View(hub_NPC_griesi)
hub_NPC_griesi

#enrichment cell markers
list_cell_markers <- read_excel("~/Desktop/MRes/20220629/cell_markers.xlsx",  col_names =TRUE)
#convert to df
df_list_cell_markers  <-as.data.frame(list_cell_markers )
View(df_list_cell_markers)
df_list_cell_markers <-df_list_cell_markers[-1,1:3]
colnames(df_list_cell_markers) <- c("species", "gene_Symbol", "category" )
df_list_cell_markers <-subset(df_list_cell_markers, species !="Mm")
df_list_cell_markers <-df_list_cell_markers[,2:3]
df_list_cell_markers$ENTREZID = mapIds(org.Hs.eg.db,
                                       keys=df_list_cell_markers$gene_Symbol, 
                                       column="ENTREZID",
                                       keytype="SYMBOL",
                                       multiVals="first")


df_list_cell_markers<- df_list_cell_markers[,c(3,2,1)]
write.csv(df_list_cell_markers,"cells_markers.csv",row.names=FALSE)

#import spreadsheet with sfari genes
Sfari_genes <- read.csv("~/Desktop/MRes/20220611/SFARI.csv", header =TRUE, stringsAsFactors = FALSE)
Sfari_genes <-Sfari_genes[,c(2,6,7)]
Sfari_genes$ENTREZID = mapIds(org.Hs.eg.db,
                              keys=Sfari_genes$gene.symbol, 
                              column="ENTREZID",
                              keytype="SYMBOL",
                              multiVals="first")
Sfari_genes<- Sfari_genes[,c(4,2,3)]
write.csv(Sfari_genes,"Sfari_genes_all.csv",row.names=FALSE)

Sfari_high_score <- subset(Sfari_genes, gene.score %in% c(1,2))


write.csv(Sfari_high_score,"Sfari_genes.csv",row.names=FALSE)

enrichments_NPC_griesi <- userListEnrichment(Gene_NPC_griesi,
                                             modules_griesi1,
                                             c("cells_markers.csv","Sfari_genes_all.csv"),
                                             c("cells_markers", "ASD_genes"), 
                                             "enrichment_NPC.csv",
                                             useBrainLists = F)

enrichments_NPC_griesi$sigOverlaps
#InputCategories                           UserDefinedCategories Type CorrectedPvalues
#1            blue                     Interneurons__cells_markers User     6.737614e-06
#2            blue                          Neurons__cells_markers User     1.034648e-05
#3       turquoise                      Fibroblasts__cells_markers User     4.228787e-05
#4           brown                  Ependymal cells__cells_markers User     4.561479e-04
#5            blue                 Immature neurons__cells_markers User     2.494331e-03
#6            blue                 Purkinje neurons__cells_markers User     2.823592e-03
#7           black Rare Single Gene Mutation, Syndromic__ASD_genes User     1.275739e-02
#8            blue                      Neuroblasts__cells_markers User     1.758488e-02
#9            pink              Gamma delta T cells__cells_markers User     3.682147e-02



#neurons from normalised counts RUVSeq #variable normalized_counts2 with outlier removed
neurons_coldata <-subset(coldata2, subset = cell_type !="NPC")
neurons_coldata2 <-neurons_coldata[-9,]
neurons_coldata2$cell_type <-NULL
neurons_coldata2$condition <-c(rep(1,8), rep(0,10)) #ASD: 1; CTL:0
neurons_coldata2$group <-NULL
datTraits_neurons_griesi <- neurons_coldata2

txi_neurons_IFN <- lapply(txi, function(x) if(is.matrix(x)) return(x[,7:18]) else return(x))
neuron_coldata_IFN <-subset(IFN_coldata, subset = cell_type =="neuron")
neuron_coldata_IFN $cell_type <-NULL
neuron_coldata_IFN $cell_line <- c(rep(c(1,2,3),4))
neuron_coldata_IFN $treatment <- c(rep(2,3), rep(1,6), rep(0,3)) #): no treatment; 1: single treatment (early or late); 2: repeated treatment
neuron_coldata_IFN $D18 <- NULL
neuron_coldata_IFN $D30 <- NULL
datTraits_neurons_IFN <- neuron_coldata_IFN 

#ASD
dds <- DESeqDataSetFromMatrix(countData = normalized_counts2,
                              colData = datTraits_neurons_griesi,
                              design = ~1)
keep <- rowSums(counts(dds)) >= 10
table(keep)
#FALSE  TRUE 
# 4987 25659   
dds <- dds[keep,]
vsd <-vst(dds, blind=FALSE)
mat_G_n <- assay(vsd)

#check data
mean_expr <-apply(mat_G_n,2, mean, na.rm =T)
number_missing <-apply(is.na(data.frame(mat_G_n)),2,sum)

sizeGrWindow(12, 8)
barplot(mean_expr,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:18))

number_missing

NumberMissingByGene <- apply( is.na(data.frame(mat_G_n)),1, sum)
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variancedatExpr<-as.vector(apply(as.matrix(mat_G_n),1,var, na.rm=T))
no.presentdatExpr<-as.vector(apply(!is.na(as.matrix(mat_G_n)),1, sum) )
table(no.presentdatExpr)

# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes <- variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)

mat_neurons_griesi <-t(mat_G_n)

#same analysis with IFN
colnames(txi_neurons_IFN$counts) <-rownames(datTraits_neurons_IFN)

dds <- DESeqDataSetFromTximport(txi_neurons_IFN,
                                colData = datTraits_neurons_IFN,
                                design = ~ 1)

keep <- rowSums(counts(dds)) >= 10
table(keep)
#FALSE  TRUE 
# 8129 22517    
dds <- dds[keep,]

vsd <-vst(dds, blind=FALSE)
mat_D_n <- assay(vsd)

#check data
mean_expr <-apply(mat_D_n,2, mean, na.rm =T)
number_missing <-apply(is.na(data.frame(mat_D_n)),2,sum)

sizeGrWindow(12, 8)
barplot(mean_expr,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:12))

number_missing

NumberMissingByGene <- apply( is.na(data.frame(mat_D_n)),1, sum)
summary(NumberMissingByGene)

# Calculate the variances of the probes and the number of present entries
variancedatExpr<- as.vector(apply(as.matrix(mat_D_n),1,var, na.rm=T))

no.presentdatExpr<- as.vector(apply(!is.na(as.matrix(mat_D_n)),1, sum) )
table(no.presentdatExpr)

# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes<- variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)

mat_neurons_deepak <-t(mat_D_n)

#limit the analysis to genes that are expressed in both datasets
commonGenes <-intersect(rownames(mat_G_n), rownames(mat_D_n ))
length(commonGenes)  #18266
mat2neurons_griesi <-mat_G_n[commonGenes,]
mat2neurons_deepak<-mat_D_n[commonGenes,]

#find appropriate soft power
powers <- c(c(1:20), seq(from=22, to=30, by=2))
sft <-pickSoftThreshold(t(mat2neurons_deepak), powerVector = powers, verbose =5)
sft <-pickSoftThreshold(t(mat2neurons_griesi), powerVector = powers, verbose =5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 <- 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.9, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <-7  #IFN
softPower <- 19 #ASD

#correlate measures of average gene expression and overall connectivity between 2 datasets
#the higher the correaltion of these properties, the better chance of finding similarities

#random sample of 5000 genes in the 10000 genes with highest variance in griesi's dataset
variancedatExpr=as.vector(apply(as.matrix(mat2neurons_griesi),1,var, na.rm=T))
variancedatExpr_order <-as.data.frame(variancedatExpr)
variancedatExpr_order$gene <-rownames(mat2neurons_griesi)
variancedatExpr_order <-variancedatExpr_order[order(variancedatExpr_order$variancedatExpr, decreasing = TRUE),]
variancedatExpr_order_10000 <-variancedatExpr_order[1:10000,]
random5000 <- sample(variancedatExpr_order_10000$variancedatExpr,5000)

rankExpr_griesi<- rank(rowMeans(mat2neurons_griesi))
rankExpr_deepak<- rank(rowMeans(mat2neurons_deepak))
rankConn_griesi<- rank(softConnectivity(t(mat2neurons_griesi[random5000,]),type="signed",power=19)) 
rankConn_deepak<- rank(softConnectivity(t(mat2neurons_deepak[random5000,]),type="signed",power=7))

sizeGrWindow(9, 5)
par(mfrow=c(1,2))
verboseScatterplot(rankExpr_griesi,rankExpr_deepak, xlab="Ranked Expression (ASD)",ylab="Ranked Expression (IFN)") 
verboseScatterplot(rankConn_griesi,rankConn_deepak, xlab="Ranked Connectivity (ASD)", ylab="Ranked Connectivity (IFN)")


#run WGCNA on the datasets
adjacency_deepak <- adjacency(t(mat2neurons_deepak),power=7,type="signed") 
diag(adjacency_deepak)<-0
dissTOM_deepak <- 1-TOMsimilarity(adjacency_deepak, TOMType="signed") 
geneTree_deepak <- flashClust(as.dist(dissTOM_deepak), method="average")

adjacency_griesi_neurons <- adjacency(t(mat2neurons_griesi),power=19,type="signed") 
diag(adjacency_griesi_neurons)<-0
TOM_griesi_neurons <-TOMsimilarity(adjacency_griesi_neurons, TOMType="signed")
dissTOM_griesi <- 1-TOM_griesi_neurons  
geneTree_griesi <- flashClust(as.dist(dissTOM_griesi), method="average")

pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree_griesi,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (griesi)", labels=FALSE,hang=0.04);
plot(geneTree_deepak,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (deepak)", labels=FALSE,hang=0.04);
dev.off()
#"good" data when a lot of distict branches

#determine modules based on griesi's dataset
#4 module splits we can choose from
mColorh <-NULL
for (ds in 0:3){
  tree <- cutreeHybrid(dendro = geneTree_griesi, pamStage=FALSE, cutHeight = 0.99, minClusterSize = (150), deepSplit = ds, distM = dissTOM_griesi)
  mColorh=cbind(mColorh,labels2colors(tree$labels)) }
pdf("Module_choices2.pdf", height=10,width=25)
plotDendroAndColors(geneTree_griesi, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE)
dev.off()
modules_griesi2 = mColorh[,3] #ds =2

#calculate the principle components for visualisation and quantitative assessemnt
#first PC= module eigengene (ME) = single value that represents the highest percent of variance for all genes in a module (1 par module et par sample)
#most genes in the module "do the same thing" as ME
PCs2_deepak <- moduleEigengenes(t(mat2neurons_deepak), colors=modules_griesi2)

ME_2_deepak <-PCs2_deepak$eigengenes #table sample/module
distPC2_deepak <-1-abs(cor(ME_2_deepak,use="p")) #matrice module-module
distPC2_deepak <-ifelse(is.na(distPC2_deepak), 0, distPC2_deepak)
pcTree2_deepak <-hclust(as.dist(distPC2_deepak),method="a")
MDS_2_deepak <-cmdscale(as.dist(distPC2_deepak),2)
colors_deepak2 <-names(table(modules_griesi2))

PCs2_griesi <- moduleEigengenes(t(mat2neurons_griesi), colors=modules_griesi2)
ME_2_griesi <- PCs2_griesi$eigengenes
colors_griesi2 <-names(table(modules_griesi2))

#associate with traits
MEs <- orderMEs(ME_2_griesi)
moduleTraitCor <- cor(MEs, datTraits_neurons_griesi, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 18)

sizeGrWindow(10,6)

#Displaying correlations and its p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 5, 5.5))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_neurons_griesi),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#correlate with IFN-g treatment
ME1s <- orderMEs(ME_2_deepak)
moduleTraitCor <- cor(ME1s, datTraits_neurons_IFN, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 12)

sizeGrWindow(10,6)

#Displaying correlations and its p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 5, 5))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_neurons_IFN),
               yLabels = names(ME1s),
               ySymbols = names(ME1s),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#how well a model in a study is preserved in another one (outputs a single Z-score summary)
# (This step will take ~10-30 minutes)
multiExpr <- list(griesi=list(data=t(mat2neurons_griesi)), deepak=list(data=t(mat2neurons_deepak)))
multiColor <- list(griesi = modules_griesi2) 
mp<-modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                       nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000) 
stats <- mp$preservation$Z$ref.griesi$inColumnsAlsoPresentIn.deepak 
stats_neurons <-stats[order(-stats[,2]),c(1:2)]
stats_neurons <-stats[order(-stats[,2]),]
View(stats_neurons)

#module membership=kME
#value used to measure correlations between each gene and each ME
# => genes not initially assigned to a module can be included in between-network comparisons
#get kME values and associated p-values
geneModuleMembership3 <- signedKME(t(mat2neurons_deepak), ME_2_deepak) 
colnames(geneModuleMembership3)<-paste("PC",colors_griesi2,".cor",sep="")

MMPvalue3<- corPvalueStudent(as.matrix(geneModuleMembership3),dim(mat2neurons_deepak)[[2]]) 
colnames(MMPvalue3)<- paste("PC",colors_griesi2,".pval",sep="")

Gene_neurons_griesi <- rownames(mat2neurons_griesi) 
kMEtable3 <- cbind(Gene_neurons_griesi,Gene_neurons_griesi,modules_griesi2) 
for (i in 1:length(colors_griesi2))
  kMEtable3 <- cbind(kMEtable3, geneModuleMembership3[,i], MMPvalue3[,i])
colnames(kMEtable3) <-c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership3), colnames(MMPvalue3))))
write.csv(kMEtable3,"kMEtable3deepakneuronsintersection.csv",row.names=FALSE)

#kME values and p-values for griesi
geneModuleMembership4 <- signedKME(t(mat2neurons_griesi), ME_2_griesi ) 
colnames(geneModuleMembership4) <-paste("PC",colors_griesi2,".cor",sep="")

MMPvalue4<- corPvalueStudent(as.matrix(geneModuleMembership4),dim(mat2neurons_griesi)[[2]])
colnames(MMPvalue4) <-paste("PC",colors_griesi2,".pval",sep="")

kMEtable4 <- cbind(Gene_neurons_griesi,Gene_neurons_griesi,modules_griesi2) 
for (i in 1:length(colors_griesi2))
  kMEtable4 <- cbind(kMEtable4, geneModuleMembership4[,i], MMPvalue4[,i]) 
colnames(kMEtable4) <- colnames(kMEtable3)
write.csv(kMEtable4,"kMEtable4griesineuronsintersection.csv",row.names=FALSE)


#which genes are hubs in both networks (genes with extremely high kME values in both netwroks)
topGenesKME <- NULL
for (c in 1:length(colors_griesi2)){
  kMErank3 <- rank(-geneModuleMembership3[,c])
  kMErank4 <- rank(-geneModuleMembership4[,c])
  maxKMErank <- rank(apply(cbind(kMErank3,kMErank4+.00001),1,max)) 
  topGenesKME <- cbind(topGenesKME,Gene_neurons_griesi[maxKMErank<=20])
} 
colnames(topGenesKME) <- colors_griesi2
View(topGenesKME)

topGeneKME <-as.data.frame(topGenesKME)

module_turquoise <-data.frame(
  gene_ID <- topGenesKME[,15],
  SYMBOL <-topGenesKME[,15]
)
View(module_turquoise)
colnames(module_turquoise) <-c('gene_ID', 'SYMBOL')

# Add gene symbol column
module_turquoise$SYMBOL <- mapIds(org.Hs.eg.db,
                                 keys=gene_ID, 
                                 column="SYMBOL",
                                 keytype="ENTREZID",
                                 multiVals="first")

#gene significance
condition <- as.data.frame(datTraits_neurons_griesi$condition)
names(condition) <- "condition"
# names (colors) of the modules
modNames <- substring(names(ME_2_griesi), 3)
names(geneModuleMembership4) <- paste("MM", modNames, sep="");
names(MMPvalue4) <- paste("p.MM", modNames, sep="");
geneTraitSignificance <- as.data.frame(cor(t(mat2neurons_griesi), condition, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 18));
names(geneTraitSignificance) <- paste("GS.", names(condition), sep="");
names(GSPvalue) <- paste("p.GS.", names(condition), sep="");

#intramodular analysis
#genes with high GS and high MM
module <- "greenyellow"
column <-match(module, modNames)
moduleGenes <- modules_griesi2 ==module

sizeGrWindow(7,7)
par(mfrow =c(1,1))
verboseScatterplot(abs(geneModuleMembership4[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab ="Gene significance for condition",
                   main = paste("module membership vs gene significance\n"),
                   cex.main =1.2, cex.lab = 1.2, cex.axis =1.2, col = module)


geneInfo0_neurons_griesi <- data.frame(ENTREZID = rownames(mat2neurons_griesi),
                                      moduleColor = modules_griesi2,
                                      geneTraitSignificance,
                                      GSPvalue)

geneInfo0_neurons_griesi$SYMBOL <- mapIds(org.Hs.eg.db,
                                         keys=rownames(geneInfo0_neurons_griesi), 
                                         column="SYMBOL",
                                         keytype="ENTREZID",
                                         multiVals="first")

modOrder <- order(-abs(cor(ME_2_griesi, condition, use = "p")))
for (mod in 1:ncol(geneModuleMembership4))
{
  oldNames <- names(geneInfo0_neurons_griesi)
  geneInfo0_neurons_griesi<- data.frame(geneInfo0_neurons_griesi, geneModuleMembership4[, modOrder[mod]], 
                                       MMPvalue4[, modOrder[mod]]);
  names(geneInfo0_neurons_griesi) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                      paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder <- order(geneInfo0_neurons_griesi$moduleColor, -abs(geneInfo0_neurons_griesi$GS.condition))
geneInfo_neurons_griesi <- geneInfo0_neurons_griesi[geneOrder, ]

white <- subset(geneInfo_neurons_griesi, moduleColor == "white")
white2  <- subset(white, p.GS.condition <0.05)
white2 <-white2[order(abs(white2$GS.condition), decreasing= TRUE),]
white3 <-white[order(abs(white$MM.white), decreasing= TRUE),]
write.csv(white3,"neurons_WGCNA_griesi_MEwhite.csv",row.names=FALSE)

white2$ENTREZID <-white2$SYMBOL
white2$SYMBOL <-NULL
colnames(white2)[1] <- "SYMBOL"

greenyellow <- subset(geneInfo_neurons_griesi, moduleColor == "greenyellow")
greenyellow2  <- subset(greenyellow, p.GS.condition <0.05)
greenyellow2 <-greenyellow2[order(abs(greenyellow2$GS.condition), decreasing= TRUE),]
greenyellow3 <-greenyellow[order(abs(greenyellow$MM.greenyellow), decreasing= TRUE),]
write.csv(greenyellow3,"neurons_WGCNA_griesi_MEgreenyellow.csv",row.names=FALSE)

greenyellow2$ENTREZID <-greenyellow2$SYMBOL
colnames(greenyellow2)[1] <- "SYMBOL"

salmon_neurons <- subset(geneInfo_neurons_griesi, moduleColor == "salmon")

salmon_neurons3 <-salmon_neurons[order(abs(salmon_neurons$MM.salmon), decreasing= TRUE),]
write.csv(salmon_neurons3,"neurons_WGCNA_griesi_MEsalmon.csv",row.names=FALSE)

turquoise_neurons <- subset(geneInfo_neurons_griesi, moduleColor == "turquoise")

turquoise_neurons2  <- subset(turquoise_neurons, p.GS.condition <0.05)
View(turquoise_neurons2)
turquoise_neurons2 <-turquoise_neurons2[order(abs(turquoise_neurons2$GS.condition), decreasing= TRUE),]
turquoise_neurons3 <-turquoise_neurons[order(abs(turquoise_neurons$MM.turquoise), decreasing= TRUE),]
write.csv(turquoise_neurons3,"neurons_WGCNA_griesi_MEturquoise.csv",row.names=FALSE)

darkred_neurons <- subset(geneInfo_neurons_griesi, moduleColor == "darkred")
darkred_neurons3 <-darkred_neurons[order(abs(darkred_neurons$MM.darkred), decreasing= TRUE),]
write.csv(darkred_neurons3,"neurons_WGCNA_griesi_darkred.csv",row.names=FALSE)

darkturquoise_neurons <- subset(geneInfo_neurons_griesi, moduleColor == "darkturquoise")
darkturquoise_neurons3 <-darkturquoise_neurons[order(abs(darkturquoise_neurons$MM.darkturquoise), decreasing= TRUE),]
write.csv(darkturquoise_neurons3,"neurons_WGCNA_griesi_MEdarkturquoise.csv",row.names=FALSE)


#GO enrichment of genes from brown module signif associated with ASD
# enrichment in term of BP of the set of  genes from interesting module
list_genes <- white3$ENTREZID
list_genes <- greenyellow3$ENTREZID
list_genes <- darkred_neurons3$ENTREZID
list_genes <- salmon_neurons3$ENTREZID
list_genes <- darkturquoise_neurons3$ENTREZID
list_genes <- turquoise_neurons3$ENTREZID

ego2_neurons <- clusterProfiler::enrichGO(gene          = list_genes,
                                          OrgDb         = OrgDb,
                                          ont           = "BP",
                                          pAdjustMethod = "fdr",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.01, 
                                          readable      = TRUE)

x2 <- pairwise_termsim(ego2_neurons)
emapplot(x2, showCategory = 12, node_label = "category")

writeLines(greenyellow3$ENTREZID, "neurons_greenyellow_genes.txt")
writeLines(turquoise_neurons3$ENTREZID, "neurons_turquoise_genes.txt")
writeLines(white3$ENTREZID, "neurons_white_genes.txt")
writeLines(salmon_neurons3$ENTREZID, "neurons_salmon_genes.txt")
writeLines(darkred_neurons3$ENTREZID, "neurons_darkred_genes.txt")
writeLines(darkturquoise_neurons3$ENTREZID, "neurons_darkturquoise_genes.txt")

#overlap DE IFN dataset1 avec MEgreenyellow  neurons griesi WGCNA
MEgreenyellow_neurons_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = greenyellow3$ENTREZID, columns ="SYMBOL", keytype = "ENTREZID")

go.obj <-newGeneOverlap(MEgreenyellow_neurons_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )  #22, p-val =0.92
go.obj <-testGeneOverlap(go.obj)
go.obj
go.obj <-newGeneOverlap(MEgreenyellow_neurons_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 21581 ) #13, 0.54
go.obj <-testGeneOverlap(go.obj)
go.obj

#overlap DE IFN dataset1 avec MEturquoise  neurons griesi WGCNA
MEturquoise_neurons_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = turquoise_neurons3$ENTREZID, columns ="SYMBOL", keytype = "ENTREZID")

go.obj <-newGeneOverlap(MEturquoise_neurons_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 21581 ) #94, 1.9e-10
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap_MEturquoise_DE_neurons_IFN <-getIntersection(go.obj)
overlap_MEturquoise_DE_neurons_IFN
#[1] "GPRIN3"    "RAB27B"    "LGI3"      "RASGRF1"   "PPP1R1A"   "VSNL1"     "SYT13"     "PEG10"     "RIMKLA"    "SLC12A5"   "LDB2"     
#[12] "RAP1GAP2"  "SYT16"     "RAB11FIP4" "RAB3C"     "RAB3B"     "NETO1"     "GABRG2"    "DISP2"     "SPHKAP"    "GREB1"     "ARG2"     
#[23] "PLCH1"     "LRSAM1"    "SV2C"      "GABRA3"    "C22orf42"  "PGM2L1"    "SCN9A"     "INSRR"     "GNAS"      "CASP7"     "SPINT2"   
#[34] "GPC6"      "SLC2A13"   "NEFM"      "LRRTM3"    "DNM1"      "DAB1"      "RNF115"    "CDH12"     "SLC8A1"    "CDKL2"     "SHFL"     
#[45] "ADGRL3"    "CA10"      "HS3ST4"    "SEMA5A"    "ADGRB1"    "HLA-H"     "CSMD2"     "DDC"       "SLC16A7"   "TMEM255B"  "FBXO25"   
#[56] "KCNQ2"     "SIK3"      "HOOK1"     "GJD2"      "NRSN1"     "LRFN2"     "CHRNB2"    "SYNPR"     "PPP2R5B"   "NEFL"      "RALYL"    
#[67] "GALNT14"   "GRIK3"     "DGLUCY"    "SORCS1"    "NECTIN1"   "PCSK1"     "KCNAB1"    "IFIH1"     "SYT4"      "WSCD2"     "SCD"      
#[78] "CHRM3"     "NRIR"      "XKR4"      "GSTK1"     "APOBEC3G"  "XPR1"      "IGFLR1"    "HMGCS1"    "ITGAV"     "TRHDE-AS1" "SULF2"    
#[89] "NDFIP1"    "KCTD8"     "LANCL3"    "ARHGEF11"  "EFNB1"     "SOX9"  

go.obj <-newGeneOverlap(list_SFARI, overlap_MEturquoise_DE_neurons_IFN, 21581 )
go.obj <-testGeneOverlap(go.obj)
go.obj #11, 1.8e-03
getIntersection(go.obj)
#"CHRM3"   "DDC"     "GALNT14" "GNAS"    "GPC6"    "GRIK3"   "KCNQ2"   "LRFN2"   "SCN9A"   "SEMA5A"  "SLC12A5"

go.obj <-newGeneOverlap(MEturquoise_neurons_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )  #138, p-val =2.9e-04
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap_MEturquoise_DE_NPC_IFN <-getIntersection(go.obj)
overlap_MEturquoise_DE_NPC_IFN
#[1] "RIMS1"     "FBXO41"    "CKMT1B"    "PSD"       "STX1B"     "RIMS4"     "IFNLR1"    "CGREF1"    "SH3GL3"    "PLEKHB1"   "MCF2L"     "TRIM67"   
#[13] "SLC1A1"    "RIC3"      "KCNC1"     "LDB2"      "SEPTIN3"   "RAB11FIP4" "ADAM11"    "CHRNA4"    "RTN1"      "SLC1A6"    "FBXO44"    "ACSL3"    
#[25] "TUBB4A"    "GPR88"     "GPRIN1"    "KIF21B"    "CPNE7"     "KCNQ3"     "TBC1D9"    "ARG2"      "TDRD7"     "YBX1"      "FBXO2"     "PLCH1"    
#[37] "GNG4"      "ACSL6"     "PDE1B"     "ZNF385B"   "ELMOD1"    "SV2C"      "GABRA3"    "CHD5"      "DNER"      "NLGN3"     "NSG1"      "CADM3"    
#[49] "SCN9A"     "PARP8"     "CASP7"     "SPTB"      "GABRB3"    "GPC6"      "INPP1"     "SLC6A11"   "GRM4"      "FAM171B"   "LRRTM3"    "DNM1"     
#[61] "DAB1"      "ZBTB8B"    "MRAS"      "RNF115"    "CDH12"     "CPEB3"     "SLC8A1"    "SHFL"      "SYT7"      "PPFIA4"    "ADGRL3"    "CFH"      
#[73] "HS3ST4"    "SEMA5A"    "PARP6"     "ADGRB1"    "HLA-H"     "SC5D"      "TMEM255B"  "SLC29A1"   "NRXN3"     "GJD2"      "CHRNB2"    "ATP1A2"   
#[85] "HTR2C"     "TCAF1"     "PPP2R5B"   "C1R"       "ABCC5"     "GRIK3"     "PAK1"      "DGLUCY"    "SORCS1"    "PJA2"      "ARX"       "MLC1"     
#[97] "MN1"       "LINC02609" "PCSK1"     "EFR3B"     "TXLNB"     "IFIH1"     "CTSO"      "SCD"       "HSD17B7"   "KCNMB2"    "JADE2"     "ANKS1B"   
#[109] "HMGCR"     "GSTK1"     "APOBEC3G"  "B3GALT2"   "VRK2"      "IGFLR1"    "HMGCS1"    "STARD4"    "ITGAV"     "CYP51A1"   "TENM2"     "SULF2"    
#[121] "GPR50"     "IQSEC3"    "MSMO1"     "CORO1A"    "NALF2"     "FAM120B"   "PRKCQ-AS1" "FDFT1"     "GRIP2"     "ASB13"     "LHX1"      "ACSBG1"   
#[133] "SLC16A10"  "EFNB1"     "PDGFRB"    "PGBD5"     "PDLIM3"    "SOX9" 

go.obj <-newGeneOverlap(list_SFARI, overlap_MEturquoise_DE_NPC_IFN, 20068 )
go.obj <-testGeneOverlap(go.obj)
go.obj #17, 1.4e-04
getIntersection(go.obj)
#"ANKS1B"  "ARX"     "CORO1A"  "DNER"    "GABRB3"  "GPC6"    "GRIK3"   "INPP1"   "KCNC1"   "KCNQ3"   "NLGN3"   "NRXN3"   "RIMS1"   "SCN9A"  
#[15] "SEMA5A"  "SLC1A1"  "ZNF385B"

#overlap DE IFN dataset1 avec MEsalmon  neurons griesi WGCNA
MEsalmon_neurons_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = salmon_neurons3$ENTREZID, columns ="SYMBOL", keytype = "ENTREZID")

go.obj <-newGeneOverlap(MEsalmon_neurons_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )  #26, p-val =0.63
go.obj <-newGeneOverlap(MEsalmon_neurons_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 21581 ) #8, 0.94
go.obj <-testGeneOverlap(go.obj)
go.obj

#overlap DE IFN dataset1 avec MEdarkred  neurons griesi WGCNA
MEreddark_neurons_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = darkred_neurons3$ENTREZID, columns ="SYMBOL", keytype = "ENTREZID")
go.obj <-newGeneOverlap(MEreddark_neurons_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )  #15, p-val =0.76
go.obj <-newGeneOverlap(MEreddark_neurons_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 21581 ) #9, 0.49
go.obj <-testGeneOverlap(go.obj)
go.obj

#overlap DE IFN dataset1 avec MEwhite  neurons griesi WGCNA
MEwhite_neurons_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = white3$ENTREZID, columns ="SYMBOL", keytype = "ENTREZID")

go.obj <-newGeneOverlap(MEwhite_neurons_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )  #19, p-val =0.13
go.obj <-newGeneOverlap(MEwhite_neurons_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 21581 ) #5, 0.84
go.obj <-testGeneOverlap(go.obj)
go.obj

#overlap DE IFN dataset1 avec MEdarkturquoise  neurons griesi WGCNA
MEdarkturquoise_neurons_genes_symbol <-AnnotationDbi::select(org.Hs.eg.db, keys = darkturquoise_neurons3$ENTREZID, columns ="SYMBOL", keytype = "ENTREZID")

go.obj <-newGeneOverlap(MEdarkturquoise_neurons_genes_symbol$SYMBOL, DE_list_NPC_IFN$SYMBOL, 20068 )  #2, p-val =0.12
go.obj <-newGeneOverlap(MEdarkturquoise_neurons_genes_symbol$SYMBOL, DE_genes_neurons_IFN, 21581 ) #12, 0.13
go.obj <-testGeneOverlap(go.obj)
go.obj

#enrichment cell markers and sfari
enrichments_neurons_griesi <- userListEnrichment(Gene_neurons_griesi,modules_griesi2,  
                                                 c("cells_markers.csv","Sfari_genes_all.csv") , 
                                                 c("cells_markers", "ASD_genes"), 
                                                 "enrichment_neurons_griesi.csv",
                                                 useBrainLists =TRUE)

enrichments_neurons_griesi$sigOverlaps
#   InputCategories                                     UserDefinedCategories Type CorrectedPvalues
#1       lightgreen                        Gamma delta T cells__cells_markers User     6.389408e-15
#2        turquoise                                    Neurons__cells_markers User     6.600643e-13
#3           yellow                                Fibroblasts__cells_markers User     8.243809e-11
#4        turquoise                               Interneurons__cells_markers User     1.077503e-07
#5        turquoise Rare Single Gene Mutation, Genetic Association__ASD_genes User     6.685922e-07
#6           yellow                        Smooth muscle cells__cells_markers User     4.065763e-05
#7           yellow                             Myofibroblasts__cells_markers User     1.107936e-04
#8        turquoise                     Retinal ganglion cells__cells_markers User     2.825005e-04
#9      greenyellow                            Ependymal cells__cells_markers User     7.905520e-04
#10          purple           Rare Single Gene Mutation, Syndromic__ASD_genes User     1.464394e-03
#11       turquoise                            Enteric neurons__cells_markers User     1.920899e-03
#12      lightgreen                           Epithelial cells__cells_markers User     3.604475e-03
#13       turquoise                           Purkinje neurons__cells_markers User     1.099194e-02
#14       turquoise                         Trigeminal neurons__cells_markers User     1.477159e-02
#15          yellow                            Mesangial cells__cells_markers User     2.065697e-02
#16            blue                                Basal cells__cells_markers User     3.452814e-02
#17            blue                        Smooth muscle cells__cells_markers User     3.772395e-02


###module preservation between griesi NPC and griesi neurons :
enrichments_neurons_vs_NPC_griesi <- userListEnrichment(Gene_neurons_griesi, 
                                                        modules_griesi2, 
                                                        "kMEtable2griesiNPCintersection.csv",
                                                        "modules_NPC_griesi",
                                                        "enrichment_neurons_vs_NPC_griesi.csv")



#WGCNA on IFN dataset
bwnet7 <-blockwiseModules(mat_neurons_deepak, maxBlockSize = 5000,
                          power =7, TOMType = "signed", minModuleSize = 150,
                          reassignThreshold = 0, mergeCutHeight = 0.35, deepSplit=2,
                          numericLabels = TRUE,
                          saveTOMs = FALSE,
                          verbose=5)


#how many modules and their size
table(bwnet7$colors)
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15      16   17       18   19    20
#grey tur blue brown ye  green red black pink mag  purp gry tan salm cyan   midNbl  lcyan grey60 lgre lyel  royal blue
#960 4204 4166 2686 2536 1288 1114  821  788  727  545  459  308  280  256  254     252  247    245   216  165 


#display dendrogram with color assignment
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors <- labels2colors(bwnet7$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet7$dendrograms[[1]], mergedColors[bwnet7$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels <- bwnet7$colors
moduleColors <- labels2colors(bwnet7$colors)
moduleColors
colors_N_deepak <-names(table(moduleColors))
table(moduleColors)
#
unique(moduleColors)
MEList <- moduleEigengenes(mat_neurons_deepak, colors = moduleColors)
MEs <- MEList$eigengenes
MEs <- orderMEs(MEs)
moduleTraitCor <- cor(MEs, datTraits_neurons_IFN, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 12)

sizeGrWindow(10,6)

#Displaying correlations and its p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 8.5, 5, 5))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_neurons_IFN),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Identifying most important genes for one determined characteristic inside of the cluster
#GS quantifies associations of individual genes with trait of interest: abs of the correlation between the gene and the trait
#MM for each module: correlation between module eigengene and gene expression profile
geneModuleMembership_neurons_deepak <- signedKME(mat_neurons_deepak, MEs)
colnames(geneModuleMembership_neurons_deepak) <- paste("PC", colors_N_deepak, ".cor", sep="")

MMPvalue_neurons_deepak <- corPvalueStudent(as.matrix(geneModuleMembership_neurons_deepak), dim(t(mat_neurons_deepak))[[2]])
colnames(MMPvalue_neurons_deepak) <-paste("PC", colors_N_deepak, ".pval", sep="")

Gene_neuronsD <- colnames(mat_neurons_deepak)
kMEtable_neurons_deepak <- cbind(Gene_neuronsD,Gene_neuronsD,moduleColors) 
for (i in 1:length(colors_N_deepak))
  kMEtable_neurons_deepak <- cbind(kMEtable_neurons_deepak, geneModuleMembership_neurons_deepak[,i], MMPvalue_neurons_deepak[,i])
colnames(kMEtable_neurons_deepak)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership_neurons_deepak), colnames(MMPvalue_neurons_deepak))))
write.csv(kMEtable_neurons_deepak,"kMEtable_neurons_deepak.csv",row.names=FALSE)

# Define variable treatment containing the status column of datTrait
treatment <- as.data.frame(datTraits_neurons_IFN$treatment)
names(treatment) = "treatment"
# names (colors) of the modules
modNames <- substring(names(MEs), 3)
names(geneModuleMembership_neurons_deepak) = paste("MM", modNames, sep="");
names(MMPvalue_neurons_deepak) <- paste("p.MM", modNames, sep="");
geneTraitSignificance <- as.data.frame(cor(mat_neurons_deepak, treatment, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 12));
names(geneTraitSignificance) <- paste("GS.", names(treatment), sep="");
names(GSPvalue) <- paste("p.GS.", names(treatment), sep="");

#intramodular analysis
#genes with high GS and high MM
module <- "green"
column <-match(module, modNames)
moduleGenes <- moduleColors ==module

sizeGrWindow(7,7)
par(mfrow =c(1,1))
verboseScatterplot(abs(geneModuleMembership_neurons_deepak[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab ="Gene significance for condition",
                   main = paste("module membership vs gene significance\n"),
                   cex.main =1.2, cex.lab = 1.2, cex.axis =1.2, col = module)


geneInfo0ND <- data.frame(ENTREZID = colnames(mat_neurons_deepak),
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)

geneInfo0ND$SYMBOL <- mapIds(org.Hs.eg.db,
                            keys=rownames(geneInfo0ND), 
                            column="SYMBOL",
                            keytype="ENTREZID",
                            multiVals="first")

modOrder <- order(-abs(cor(MEs, treatment, use = "p")))
for (mod in 1:ncol(geneModuleMembership_neurons_deepak))
{
  oldNames <- names(geneInfo0ND)
  geneInfo0ND <- data.frame(geneInfo0ND, geneModuleMembership_neurons_deepak[, modOrder[mod]], 
                           MMPvalue_neurons_deepak[, modOrder[mod]]);
  names(geneInfo0ND) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder <- order(geneInfo0ND$moduleColor, -abs(geneInfo0ND$GS.treatment))
geneInfoND <- geneInfo0ND[geneOrder, ]

green <- subset(geneInfoND, moduleColor == "green")
View(green)
green2 <- subset(green, p.GS.treatment <0.05) #significantly correlated
View(green2)


green2$ENTREZID <-green2$SYMBOL
green2$SYMBOL <-NULL
colnames(green2)[1] <- "SYMBOL"

green3 <-green[order(green$MM.green),]
View(green3)

write.csv(green3,"neurons_WGCNA_deepak_MEgreen.csv",row.names=FALSE)
writeLines(green3$ENTREZID, "neurons_green_deepak_genes.txt")

purple <- subset(geneInfoND, moduleColor == "purple")
purple2 <-purple[order(purple$MM.purple),]
View(purple2)

write.csv(purple2,"neurons_WGCNA_deepak_MEpurple.csv",row.names=FALSE)
writeLines(purple2$ENTREZID, "neurons_purple_deepak_genes.txt")

grey60 <- subset(geneInfoND, moduleColor == "grey60")
grey60_2 <-grey60[order(grey60$MM.grey60),]
View(grey60_2)

write.csv(grey60_2,"neurons_WGCNA_deepak_MEgrey60.csv",row.names=FALSE)
writeLines(grey60_2$ENTREZID, "neurons_grey60_deepak_genes.txt")

#GO enrichment of genes from brown module signif associated with ASD
# enrichment in term of BP of the set of  genes from interesting module
neurons_green <- green3$ENTREZID
neurons_list <- purple2$ENTREZID
neurons_list <- grey60_2$ENTREZID

ego2_neurons <- clusterProfiler::enrichGO(gene          = neurons_list,
                                          OrgDb         = OrgDb,
                                          ont           = "BP",
                                          pAdjustMethod = "fdr",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.01, 
                                          readable      = TRUE)

#no enriched term found...

#enrichment cell markers and sfari
enrichments_neurons_deepak <- userListEnrichment(Gene_neuronsD ,moduleColors,  
                                                 c("cells_markers.csv","Sfari_genes_all.csv") , 
                                                 c("cells_markers", "ASD_genes"), 
                                                 "enrichment_neurons_deepak.csv",
                                                 useBrainLists =FALSE)

enrichments_neurons <- userListEnrichment(Gene_neuronsD ,moduleColors,  
                                          "cells_markers.csv" , 
                                          "cells_markers", 
                                          "enrichment_neurons.csv",
                                          useImmunePathwayLists  =TRUE)

enrichments_neurons <- userListEnrichment(Gene_neuronsD ,moduleColors,  
                                          "cells_markers.csv" , 
                                          "cells_markers", 
                                          "enrichment_neurons.csv")



enrichments_neurons_deepak$sigOverlaps
enrichments_neurons$pValues
enrichments_neurons$ovGenes$"green -- cells_markers"

enrichments_neurons$ovGenes

##VISANT
df_genes_neurons_IFN <- data.frame(
  ENTREZID <-Gene_neuronsD
)
df_genes_neurons_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                     keys=ENTREZID, 
                                     column="SYMBOL",
                                     keytype="ENTREZID",
                                     multiVals="first")

colnames(df_genes_neurons_IFN) <-c("ENTREZID", "SYMBOL")


adjacency_neurons_deepak = adjacency(mat_neurons_deepak,power=7,type="signed") 
diag(adjacency_neurons_deepak)=0
TOM_deepak_neurons = TOMsimilarity(adjacency_neurons_deepak, TOMType="signed") 
module ="green"
genes <-Gene_neuronsD
inModule <-(moduleColors==module)
modGenes <-genes[inModule]
modTOM <-TOM_deepak_neurons[inModule, inModule]
dimnames(modTOM) <-list(modGenes,modGenes)
nTOP =20
IMConn <-softConnectivity(mat_neurons_deepak[,modGenes])
top <-(rank(-IMConn) <= nTOP)

vis_neurons_deepak20 <-exportNetworkToVisANT(
  modTOM[top,top],
  file =paste("VisantInput_neurons_deepak", module, "top20.txt", sep=""),
  weighted = TRUE, 
  threshold = 0,
  probeToGene = df_genes_neurons_IFN)


#overlap between DE genes Gandal and WGCNA IFN
#overlap with genes from MEgreen WGCNA 
go.obj <-newGeneOverlap(green3$SYMBOL,genes_DE_gandal , 21642 )
go.obj <-testGeneOverlap(go.obj)
go.obj
#listA size=1287
#listB size=1611
#Intersection size=123
#Overlapping p-value=2.3e-03
#Jaccard Index=0.0
overlap_Gandal_MEgreen <- getIntersection(go.obj)
overlap_Gandal_MEgreen
#"ATP6V1A"    "FNBP1L"     "SLC38A1"    "EPS15"      "SPTSSB"     "KIF17"      "RALGAPA1"   "SPIN3"      "GLCE"       "ZNF836"    
# "TRIM37"     "UBE2N"      "SBK1"       "CABP1"      "TUBG1"      "CTSZ"       "MIR600HG"   "KCNAB3"     "USP46"      "TOR1B"     
# "FNDC5"      "DCUN1D2"    "TACR3"      "KDM1B"      "R3HDM2"     "P4HA1"      "CD99L2"     "SEMA7A"     "SCRT1"      "SESN2"     
# "CAND1"      "PLCB3"      "RNF165"     "IRGQ"       "PABPC5-AS1" "AKAP12"     "TGFBR3"     "TRIM27"     "GMPS"       "C1S"       
#"PPP2R5A"    "C1QL1"      "SYNPO2"     "KCNK13"     "PSMB2"      "ATP10B"     "LAMB1"      "LGI3"       "STX1B"      "EID2"      
#"STXBP5-AS1" "NMI"        "AK4"        "PLEKHG2"    "PXDC1"      "USP33"      "RRAGD"      "RAP2B"      "SLC15A3"    "CMC2"      
#"LRCH2"      "AK3"        "CTH"        "KALRN"      "NFIL3"      "SMG1"       "IFI16"      "TNIP2"      "TRIM14"     "TRIM56"    
# "GADD45G"    "TRIM5"      "SP110"      "PSMA4"      "LAP3"       "CACHD1"     "IFIH1"      "HERC5"      "FERMT2"     "ATR"       
# "LZTS1"      "SCN1A"      "IRF7"       "DDX60L"     "TXLNA"      "APOBEC3F"   "BST2"       "PML"        "MOV10"      "USP18"     
#"SERPING1"   "HLA-B"      "DDX58"      "EIF2AK2"    "PARP12"     "CLMP"       "SPATS2L"    "PLSCR1"     "SCN9A"      "TRIM25"    
#"TRIM21"     "SBNO2"      "STAT3"      "OAS3"       "PARP9"      "PARP14"     "DNPEP"      "DTX3L"      "IFI6"       "BTN3A1"    
#"CSF1"       "GBP1"       "PARP10"     "ZC3HAV1"    "BTN3A3"     "TAP1"       "PSME2"      "APOL6"      "UBE2L6"     "PSMB8"     
#"HELZ2"      "IRF9"       "ISG15" 

##overlap post-mortem WGCNA IFN MEgreen
overlap_N <- NULL
overlap1_N <- NULL
overlap2_N <- NULL

for (name in green3$SYMBOL) {
  if (name %in% genes_DE_gandal){
    overlap1_N <- c(overlap1_N, name)
  }
}

for (name in overlap1_N) {
  if (name %in% list_DE_Gupta){
    overlap2_N <- c(overlap2_N, name)
  }
}

for (name in overlap2_N) {
  if (name %in% genes_DE_haney2){
    overlap_N <- c(overlap_N, name)
  }
}

overlap_N






