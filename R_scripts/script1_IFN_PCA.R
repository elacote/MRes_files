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

#BiocManager::install("PCAtools")
library("PCAtools")

#BiocManager::install("clusterProfiler")
library("clusterProfiler")

#BiocManager::install("ggnewscale")
library("ggnewscale")

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)


#################################
#directory with the files containing the abundance data
dir <- "~/Desktop/MRes/20210402/IFN_from_kallisto/"
list.files(dir)

#names of the folders containing the abundance data
samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)
samples

#list of files to be used to get abundance data for all the samples
files <- file.path(dir,  samples$V1, "abundance.tsv")
names(files)
names(files) <- paste0("sample", 1:18)
all(file.exists(files))

#create a df with a correspondance between transcripts and genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

#create an object with the abundance, counts and length, as well as the correspondance between transcripts and gene
#gene counts are estimated from the abundance data
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)


#PCA on the complete dataset
#import matrix with metadata
deepak_coldata <- readRDS("~/Desktop/MRes/20210214/coldata1.rds")
#obtain a matrix with only a treatment column instead of D18 and D30
deepak_coldata2 <- deepak_coldata
deepak_coldata2$treatment <-c("T", "T", "T", "U", "U", "U", "TT", "TT", "TT", "TU", "TU", "TU", "UT", "UT", "UT", "UU", "UU", "UU")
deepak_coldata2$D18 <-NULL
deepak_coldata2$D30 <-NULL

#convert the values as factors
deepak_coldata2$cell_line <- factor(deepak_coldata2$cell_line, levels = c("M1","M2", "M3"))
deepak_coldata2$treatment <- factor(deepak_coldata2$treatment, levels = c("U","T", "UU", "UT", "TU", "TT"))
deepak_coldata2$cell_type <- factor(deepak_coldata2$cell_type, levels = c("NPC","neuron"))

rownames(deepak_coldata2) <-colnames(txi$counts)

#construct a DESeqDataSet 
dds <- DESeqDataSetFromTximport(txi,
                                    colData = deepak_coldata2,
                                    design = ~  cell_line  + treatment)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
vsd <- vst(dds,blind=TRUE)

pcaData <- plotPCA(vsd, intgroup = c( "cell_line", "treatment"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment, shape = cell_line )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA on vst: neurons and NPCs with or without IFN treatment")

#select NPCs only
txi_NPC <- lapply(txi, function(x) if(is.matrix(x)) return(x[,1:6]) else return(x))

NPC_coldata <-subset(deepak_coldata, subset = cell_type =="NPC")
NPC_coldata$D30 <-NULL
NPC_coldata$cell_type <-NULL

#convert the data as factors
NPC_coldata$cell_line <- factor(NPC_coldata$cell_line, levels = c("M1","M2", "M3"))
NPC_coldata$D18 <- factor(NPC_coldata$D18, levels = c("untreated","treated"))

rownames(NPC_coldata) <-colnames(txi_NPC$counts)

#construct a DESeqDataSet 
dds_NPC1 <- DESeqDataSetFromTximport(txi_NPC,
                                   colData = NPC_coldata,
                                   design = ~  cell_line + D18)

#prefiltering: only keep genes with at least 10 reads aligned
#Gene removal with less than 10 reads on all the conditions
keep <- rowSums(counts(dds_NPC1)) >= 10
table(keep)
dds_NPC1 <- dds_NPC1[keep,]
dds_NPC1 <- estimateSizeFactors(dds_NPC1)

# Get log2 counts
vsd <- vst(dds_NPC1,blind=TRUE)

#PCA combining colours and shapes
pcaData <- plotPCA(vsd, intgroup = c( "cell_line", "D18"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = D18, shape = cell_line )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs")


#analysis with pcatools
#removing the lower 10% of variables based on variance
vsd2 <-assay(vsd)
p<-pca(vsd2, metadata=NPC_coldata, removeVar=0.1)

#replace ENTREZ ids with genes symbols
newnames <- mapIds(org.Hs.eg.db,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype = 'ENTREZID')
newnames <- ifelse(is.na(newnames) | duplicated(newnames),names(newnames), newnames)
rownames(p$loadings) <- newnames

biplot(p, x= "PC1", y="PC2")
biplot(p, x= "PC3", y="PC4")
             
#get the names of the genes contributing the most to PC1
genes_PC1 <- getLoadings(p, components="PC1")
ord3 <-genes_PC1[order(abs(genes_PC1$PC1), decreasing=TRUE), , drop=FALSE]
head(ord3, 20)

#heatmap for NPC on the 50 genes with highest variance whether or not they are DE
mat_NPC <- assay(vsd)

variancedatExpr=as.vector(apply(as.matrix(mat_NPC),1,var, na.rm=T))
variancedatExpr_order <-as.data.frame(variancedatExpr)
variancedatExpr_order$gene <-rownames(mat_NPC)
variancedatExpr_order <-variancedatExpr_order[order(variancedatExpr_order$variancedatExpr, decreasing = TRUE),]
variancedatExpr_order_50 <-variancedatExpr_order[1:50,]

variancedatExpr_order_50$SYMBOL = mapIds(org.Hs.eg.db,
                                         keys=variancedatExpr_order_50$gene, 
                                         column="SYMBOL",
                                         keytype="ENTREZID",
                                         multiVals="first")

mat <- counts(dds_NPC1, normalized=T)[variancedatExpr_order_50$gene,]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <-c("M1_T", "M2_T", "M3_T", "M1_U", "M2_U", "M3_U")

h_NPC <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), 
                 name ="Z-score", row_labels =variancedatExpr_order_50$SYMBOL)
h_NPC

#pca on neurons
#select neurons only
txi_neurons <- lapply(txi, function(x) if(is.matrix(x)) return(x[,7:18]) else return(x))

neuron_coldata <-subset(deepak_coldata, subset = cell_type =="neuron")
neuron_coldata$cell_type <-NULL


rownames(neuron_coldata) <- c("M1_repeated", "M2_repeated", "M3_repeated", 
                              "M1_early", "M2_early", "M3_early", 
                              "M1_late", "M2_late", "M3_late", "M1_untreated", "M2_untreated", "M3_untreated")

#convert the data as factors
neuron_coldata$cell_line <- factor(neuron_coldata$cell_line, levels = c("M1","M2", "M3"))
neuron_coldata$D18 <- factor(neuron_coldata$D18, levels = c("untreated", "treated"))
neuron_coldata$D30 <- factor(neuron_coldata$D30, levels = c("untreated", "treated"))

colnames(txi_neurons$counts) <-rownames(neuron_coldata)

#construct a DESeqDataSet 
dds_neuron1 <- DESeqDataSetFromTximport(txi_neurons,
                                       colData = neuron_coldata,
                                       design = ~  cell_line + D18 + D30)


#prefiltering: only keep genes with at least 10 reads aligned
keep <- rowSums(counts(dds_neuron1)) >= 10
table(keep)
dds_neuron1 <- dds_neuron1[keep,]
dds_neuron1 <- estimateSizeFactors(dds_neuron1)

# Get log2 counts
vsd <- vst(dds_neuron1,blind=TRUE)
#PCA combining colours and shapes
pcaData <- plotPCA(vsd, intgroup = c( "cell_line", "D30"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = D30, shape = cell_line )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons")


#analysis with pcatools
#removing the lower 10% of variables based on variance
vsd2 <-assay(vsd)
rownames(neuron_coldata) <- c("M1_repeated", "M2_repeated", "M3_repeated", 
                              "M1_early", "M2_early", "M3_early", 
                              "M1_late", "M2_late", "M3_late", "M1_untreated", "M2_untreated", "M3_untreated")

p<-pca(vsd2, metadata=neuron_coldata, removeVar=0.1)

#replace ENTREZ ids with genes symbols
newnames <- mapIds(org.Hs.eg.db,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype = 'ENTREZID')
newnames <- ifelse(is.na(newnames) | duplicated(newnames),names(newnames), newnames)
rownames(p$loadings) <- newnames

biplot(p, x="PC1", y="PC2")
biplot(p, x="PC3", y="PC4")
       
#get the names of the genes contributing the most to PC1
genes_PC1 <- getLoadings(p, components="PC1")
ord3 <-genes_PC1[order(abs(genes_PC1$PC1), decreasing=TRUE), , drop=FALSE]
head(ord3, 20)

#heatmap on neuron samples with the 50 genes with most variance 
txi_neurons <- lapply(txi, function(x) if(is.matrix(x)) return(x[,7:18]) else return(x))
deepak_coldata <- readRDS("~/Desktop/MRes/20210214/coldata1.rds")
neurons_coldata <-subset(deepak_coldata, subset = cell_type !="NPC")
neurons_coldata$cell_type <-NULL
neurons_coldata$treatment <- c(rep("treated", 9), rep("untreated", 3))
neurons_coldata$D18 <-NULL
neurons_coldata$D30 <-NULL

colnames(txi_neurons$counts) <- rownames(neurons_coldata)

dds_neurons2 <- DESeqDataSetFromTximport(txi_neurons,
                                        colData = neurons_coldata,
                                        design = ~ cell_line + treatment)

keep <- rowSums(counts(dds_neurons2)) >= 10
dds_neurons2 <- dds_neurons2[keep,]
dds_neurons2 <- estimateSizeFactors(dds_neurons2)

vsd <- vst(dds_neurons2,blind=TRUE)
mat_neurons <- assay(vsd)

variancedatExpr=as.vector(apply(as.matrix(mat_neurons),1,var, na.rm=T))
variancedatExpr_order <-as.data.frame(variancedatExpr)
variancedatExpr_order$gene <-rownames(mat_neurons)
variancedatExpr_order <-variancedatExpr_order[order(variancedatExpr_order$variancedatExpr, decreasing = TRUE),]
variancedatExpr_order_50 <-variancedatExpr_order[1:50,]

variancedatExpr_order_50$SYMBOL = mapIds(org.Hs.eg.db,
                                         keys=variancedatExpr_order_50$gene, 
                                         column="SYMBOL",
                                         keytype="ENTREZID",
                                         multiVals="first")

mat <- counts(dds_neurons2, normalized=T)[variancedatExpr_order_50$gene,]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <-c("M1_TT", "M2_TT", "M3_TT", "M1_TU", "M2_TU", "M3_TU", "M1_UT", "M2_UT", "M3_UT", "M1_UU", "M2_UU", "M3_UU")

h_neurons <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), 
                     name ="Z-score", row_labels =variancedatExpr_order_50$SYMBOL)
h_neurons
