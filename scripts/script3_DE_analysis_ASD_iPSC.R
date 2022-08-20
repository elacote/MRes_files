#BiocManager::install("tximportData")
library(tximportData)
library(tximport)

#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") 
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

#BiocManager::install("GenomicFeatures")
library("GenomicFeatures")

#BiocManager::install("PCAtools")
library("PCAtools")

#BiocManager::install("DESeq2")
library("DESeq2")

#BiocManager::install("ggplot2")
library("ggplot2")

#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

#BiocManager::install("AnnotationDbi")
library("AnnotationDbi")

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#BiocManager::install("sva")
library("sva")

#BiocManager::install("variancePartition")
library(variancePartition)

#BiocManager::install("EDASeq")
library(EDASeq)

#BiocManager::install("RUVSeq")
library("RUVSeq")

#BiocManager::install("clusterProfiler")
library("clusterProfiler")
library(enrichplot)

#install.packages("readxl")
library(readxl)

###########################
#directory with the files containing the abundance data
dir <- "~/Desktop/MRes/20220412/kal_results/"

#names of the folders containing the abundance data
samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)

#list of files to be used to get abundance data for all the samples
files <- file.path(dir,  samples$V1, "abundance.tsv")
names(files) <- samples$V1
all(file.exists(files))

#create a df with a correspondance between transcripts and genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

#create an object with the abundance, counts and length, as well as the correspondance between transcripts and gene
#gene counts are estimated from the abundance data
txi_griesi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

coldata2 <- readRDS("~/Desktop/MRes/20220412/coldata2.rds")
#add the info about the individual from whom each sample originates
individual <-c("F2613", "F2613", "F2613", "F2688", "F2688", "F2688", "F2735", "F2735", "F2735","F6136", "F7511", "F7511", "F7511", "F6281", "F6281", "F6281", "F5541", "F5541", "F5541", "F5618", "F5618", "F6119", "F6119","F7315", "F7315", "F7647", "F7647", "F7647","F7672", "F2613","F2613","F2735", "F2735", "F2735","F6281", "F6281", "F6281", "F7511", "F5541", "F5541", "F5618", "F5618", "F7647", "F7647","F7672","F7672","F7672","F7672")   
coldata2$individual <- individual

colnames(txi_griesi$counts) <-rownames(coldata2)

#analysis for NPCs
txi_griesi_NPC <- lapply(txi_griesi, function(x) if(is.matrix(x)) return(x[,1:29]) else return(x))
NPC_coldata <-subset(coldata2, subset = cell_type =="NPC")
NPC_coldata$cell_type <-NULL
NPC_coldata$condition <-factor(NPC_coldata$condition, levels = c("CTL","ASD"))
NPC_coldata$library_batch <-factor(NPC_coldata$library_batch, levels = c(1, 2, 3))
#add the info on the group, depending on the proportion of neurons in the samples, and found in the publication
NPC_coldata$group <-factor(NPC_coldata$group, levels = c("A","B", "C", "D"))

dds_NPC <- DESeqDataSetFromTximport(txi_griesi_NPC,
                                    colData = NPC_coldata,
                                    design = ~ group + library_batch + condition)
#filter genes with low expression
keep <- rowSums(counts(dds_NPC)) >= 10
table(keep)
#FALSE  TRUE 
#5179 25467
dds_NPC <- dds_NPC[keep,]
#boxplot(log10(assay(dds_NPC)))
dds_NPC <- estimateSizeFactors(dds_NPC)

# Get log2 counts
vsd_NPC <- vst(dds_NPC,blind=TRUE)
# Check distributions of samples using boxplots
boxplot(assay(vsd_NPC), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(assay(vsd_NPC)), col="blue")

#heatmap with the 50 genes with most variance
mat_NPC <- assay(vsd_NPC)
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

mat <- counts(dds_NPC, normalized=T)[variancedatExpr_order_50$gene,]
mat.z <- t(apply(mat, 1, scale))

colnames(mat.z) <-c("ASD_1_B1", "ASD_1_B2", "ASD_1_C2", 
                    "ASD_2_C1", "ASD_2_C2", "ASD_2_C2", 
                    "ASD_3_B1", "ASD_3_C2", "ASD_3_D2", 
                    "ASD_4_A1", "ASD_5_B1", "ASD_5_C2", "ASD_5_A2",
                    "ASD_6_C1", "ASD_6_C2", "ASD_6_C2",
                    "CTL_7_D2", "CTL_7_D2", "CTL_7_B1",
                    "CTL_8_C2", "CTL_8_B2", "CTL_9_B2", "CTL_9_B3",
                    "CTL_10_D2", "CTL_10_C1", "CTL_11_C2", "CTL_11_C2", "CTL_11_C1", "CTL_12_B1")

h_NPC_ASD <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), 
             name ="Z-score", row_labels =variancedatExpr_order_50$SYMBOL, row_names_gp = gpar(fontsize = 10))
h_NPC_ASD


#########
#PCA before any mormalisation
#condition and proportion of neurons in sample
pcaData <- plotPCA(vsd_NPC, intgroup = c( "group", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs- neuron proportion and status")

#condition and library batch
pcaData <- plotPCA(vsd_NPC, intgroup = c( "library_batch", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = library_batch )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs- library batch and status")

#individual and condition
pcaData <- plotPCA(vsd_NPC, intgroup = c( "individual", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = individual, shape =condition )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs- individuals")

#analysis with pcatools
#removing the lower 10% of variables based on variance
vsd2 <-assay(vsd_NPC)
NPC_coldata2 <-NPC_coldata
rownames(NPC_coldata2) <- c("ASD_F2613_B1", "ASD_F2613_B2", "ASD_F2613_C2", 
                           "ASD_F2688_C1", "ASD_F2688_C2", "ASD_F2688_C2_2", 
                           "ASD_F2735_B1", "ASD_F2735_C2", "ASD_F2735_D2", 
                           "ASD_F6136_A1", "ASD_F7511_B1", "ASD_F7511_C2", "ASD_F7511_A2",
                           "ASD_F6281_C1", "ASD_F6281_C2", "ASD_F6281_C2_2",
                           "CTL_F5541_D2", "CTL_F5541_D2_2", "CTL_F5541_B1",
                           "CTL_F5618_C2", "CTL_F5618_B2", "CTL_F6119_B2", "CTL_F6119_B3",
                           "CTL_F7315_D2", "CTL_F7315_C1", "CTL_F7647_C2", "CTL_F7647_C2_2", "CTL_F7647_C1", "CTL_F7672_B1")
colnames(vsd2) <-rownames(NPC_coldata2)

p<-pca(vsd2, metadata=NPC_coldata2, removeVar=0.1)

#replace ENTREZ ids with genes symbols
newnames <- mapIds(org.Hs.eg.db,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype = 'ENTREZID')
newnames <- ifelse(is.na(newnames) | duplicated(newnames),names(newnames), newnames)
rownames(p$loadings) <- newnames

biplot(p, x="PC3", y="PC4")
biplot(p, x="PC5", y="PC6")


#normalisation with ComBat
NPC_cts <- txi_griesi_NPC$counts

filter <- apply(NPC_cts, 1, function(x) length(x[x>10])>=12)
filtered <- NPC_cts[filter,]
NPC_cts2 <- apply(filtered, 2, function(x) as.integer(x))
rownames(NPC_cts2) <-rownames(filtered)

#normalisation for neuron_proportion
count_matrix <- as.matrix(NPC_cts2)
group <- c("B", "B", "C", "C", "C", "C", "B", "C", "D", "A", "B", "C", "A", "C", "C", "C", "D", "D", "B", "C", "B","B", "B", "D", "C", "C", "C", "C", "B")
condition <- c(rep("ASD", 16), rep("CTL", 13))
adjusted_counts_NPC_griesi <- ComBat_seq(count_matrix, batch=group, group=condition)

dds_NPC <- DESeqDataSetFromMatrix(countData = adjusted_counts_NPC_griesi,
                                  colData = NPC_coldata,
                                  design = ~library_batch + group + condition)

dds_NPC <- estimateSizeFactors(dds_NPC)
vsd_NPC <- vst(dds_NPC,blind=TRUE)

pcaData <- plotPCA(vsd_NPC, intgroup = c("condition", "group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs after CombatSeq- condition and neuron proportion")

#normalisation for library batch
count_matrix2 <- as.matrix(NPC_cts2)[,-23]
batch <- c(1,2,2,1,2,2,1,2,2,1,1,2,2,1,2,2,2,2,1,2,2,2,2,1,2,2,1,1)
condition <- c(rep("ASD", 16), rep("CTL", 12))
adjusted_counts_NPC_griesi2 <- ComBat_seq(count_matrix2, batch=batch, group=condition)

dds_NPC2 <- DESeqDataSetFromMatrix(countData = adjusted_counts_NPC_griesi2,
                                  colData = NPC_coldata[-23,],
                                  design = ~library_batch + group + condition)

dds_NPC2 <- estimateSizeFactors(dds_NPC2)
vsd_NPC2 <- vst(dds_NPC2,blind=TRUE)
pcaData <- plotPCA(vsd_NPC2, intgroup = c("condition", "library_batch"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =library_batch )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs after CombatSeq- condition and library")

#normalisation for library batch after neuron_proportion - not better!
count_matrix2 <-adjusted_counts_NPC_griesi
batch <- c(1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1)
condition <- c(rep("ASD", 16), rep("CTL", 13))
adjusted_counts2 <- ComBat_seq(count_matrix2, batch=batch, group=condition)

dds_NPC3 <- DESeqDataSetFromMatrix(countData = adjusted_counts2,
                                   colData = NPC_coldata,
                                   design = ~library_batch + group + condition)

dds_NPC3 <- estimateSizeFactors(dds_NPC3)

vsd_NPC3 <- vst(dds_NPC3,blind=TRUE)
pcaData3 <- plotPCA(vsd_NPC3, intgroup = c("condition", "group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData3, "percentVar"))
ggplot(pcaData3, aes(x = PC1, y = PC2, color = condition, shape =group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs after CombatSeq- condition and neuron proportion")


# Specify parallel processing parameters
#necessary for dream()
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

geneExpr <-  DGEList(adjusted_counts_NPC_griesi)
geneExpr <- calcNormFactors( geneExpr)
#categorical variable to be modeled as a random effect
# The variable to be tested must be a fixed effect
form3 <- ~ 0 + condition +(1|individual) + library_batch

# estimate weights using linear mixed model of dream
#Transform count data to log2-counts per million (logCPM), 
#estimate the mean-variance relationship 
#and use this to compute appropriate observation-level weights. 
#The data are then ready for linear mixed modelling with dream(). 
#This method is the same as limma::voom(), except that it allows random effects in the formula
vobjDream = voomWithDreamWeights(geneExpr, form3, NPC_coldata)
L = getContrast( vobjDream, form3, NPC_coldata, c("conditionASD", "conditionCTL"))

# fit dream model with contrasts
#Fit linear mixed model for differential expression 
#and preform hypothesis test on fixed effects as specified in the contrast matrix L
fit = dream( vobjDream, form3, NPC_coldata, L)
fitm = eBayes( fit )

top.table <- topTable( fitm, number= 100, coef ="L1")
top.table$SYMBOL = mapIds(org.Hs.eg.db,
                                  keys=rownames(top.table), 
                                  column= "SYMBOL",
                                  keytype="ENTREZID",
                                  multiVals="first")
top.table_NPC_Combat <-top.table[, c(8,1:7)]
write.csv(top.table_NPC_Combat,"DE_NPC_griesi_dream_Combat.csv",row.names=FALSE)



###### normalisation with RUVSeq
genes <- rownames(NPC_cts2)

#store the data in an object of S4 class SeqExpressionSet from the EDASeq package
b<- c(rep(c("ASD"), 16), rep(c("CTL"),13))
x <- as.factor(b)
set <- newSeqExpressionSet(as.matrix(NPC_cts2), phenoData = data.frame(x, row.names=colnames(NPC_cts2)))

#RLE = log-ratio of read count to median read count across sample
plotRLE(set, outline=FALSE, ylim=c(-4, 4))
plotPCA(set, cex=1)

#normalise the data using upper-quartile (UQ) normalization
set <- betweenLaneNormalization(set, which="upper")
plotPCA(set, cex=0.7)
plotRLE(set, outline=FALSE, ylim=c(-4, 4))
#After upper-quartile normalization, some samples still show extra variability 


#obtain a set of “in-silico empirical” negative controls, 
#e.g., least significantly DE genes based on a first-pass DE analysis 
#performed prior to RUVg normalization.
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

##
#plot pca for k= 1 to k=5; do the same analysis and keep the value of k that leads to the best clustering of ASD samples
set1 <- RUVg(set, empirical, k=1)
new_data <-pData(set1)
plotPCA(set1, cex=0.6)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4))

counts1 <-counts(set1)
normalized_counts <- normCounts(set1)

individual <-c("F2613", "F2613", "F2613", "F2688", "F2688", "F2688", "F2735", "F2735", "F2735","F6136", "F7511", "F7511", "F7511", "F6281", "F6281", "F6281", "F5541", "F5541", "F5541", "F5618", "F5618", "F6119", "F6119","F7315", "F7315", "F7647", "F7647", "F7647","F7672", "F2613","F2613","F2735", "F2735", "F2735","F6281", "F6281", "F6281", "F7511", "F5541", "F5541", "F5618", "F5618", "F7647", "F7647","F7672","F7672","F7672","F7672")   
individual_NPC <-individual[1:29]
new_data$individual <- individual_NPC 
#pca
dds_count_norm <- DESeqDataSetFromMatrix(countData = normalized_counts,
                                         colData = NPC_coldata,
                                         design = ~library_batch + group + condition)

keep <- rowSums(counts(dds_count_norm)) >= 10
table(keep)
dds_count_norm <- dds_count_norm[keep,]
dds_count_norm <- estimateSizeFactors(dds_count_norm)
vsd_NPC_norm <- vst(dds_count_norm,blind=TRUE)

pcaData <- plotPCA(vsd_NPC_norm, intgroup = c("condition", "group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for NPCs after normalizing with RUVSeq")


#we could not find a value of k leading to a "good" PCA => no assessement of DE with dream

#analysis for neurons and mixed cells
txi_griesi_neuron <- lapply(txi_griesi, function(x) if(is.matrix(x)) return(x[,30:48]) else return(x))

neurons_coldata <-subset(coldata2, subset = cell_type !="NPC")

neurons_coldata$condition <-factor(neurons_coldata$condition, levels = c("CTL","ASD"))
neurons_coldata$group <-factor(neurons_coldata$group, levels = c("F", "E", "D"))
neurons_coldata$cell_type <- factor(neurons_coldata$cell_type, levels = c("neuron","mixed"))

colnames(txi_griesi_neuron$counts) <-rownames(neurons_coldata)

dds_neurons <- DESeqDataSetFromTximport(txi_griesi_neuron,
                                        colData = neurons_coldata,
                                        design = ~ group + cell_type + condition)

keep <- rowSums(counts(dds_neurons)) >= 10
table(keep)
#FALSE  TRUE 
#4987 25659  
dds_neurons <- dds_neurons[keep,]

#normalize for library size
dds_neurons <- estimateSizeFactors(dds_neurons)

# Get log2 counts
vsd_neurons <- vst(dds_neurons,blind=TRUE)
# Check distributions of samples using boxplots
boxplot(assay(vsd_neurons), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(assay(vsd_neurons)), col="blue")

#heatmap with the 50 genes with most variance
mat_neurons <- assay(vsd_neurons)
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

mat <- counts(dds_neurons, normalized=T)[variancedatExpr_order_50$gene,]
mat.z <- t(apply(mat, 1, scale))

colnames(mat.z) <-c("ASD_1_F4", "ASD_1_F5", 
                    "ASD_3_E4", "ASD_3_F4", "ASD_3_E5", 
                    "ASD_6_E4", "ASD_6_D4", "ASD_6_F5", "ASD_5_D4",
                    "CTL_7_D4", "CTL_7_F5",
                    "CTL_8_E4", "CTL_8_F5", "CTL_11_E4", "CTL_11_D4",
                    "CTL_12_E4", "CTL_12_D4", "CTL_12_F5", "CTL_12_F5")

h_neurons_ASD <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), 
             name ="Z-score", row_labels =variancedatExpr_order_50$SYMBOL, row_names_gp = gpar(fontsize = 10))
h_neurons_ASD

#PCA before normalisation
pcaData <- plotPCA(vsd_neurons, intgroup = c( "group", "condition"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons")

pcaData <- plotPCA(vsd_neurons, intgroup = c( "cell_type", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = cell_type )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons")

#analysis with pcatools
#removing the lower 10% of variables based on variance
vsd2 <-assay(vsd_neurons)
p<-pca(vsd2, metadata=neurons_coldata, removeVar=0.1)

#replace ENTREZ ids with genes symbols
newnames <- mapIds(org.Hs.eg.db,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype = 'ENTREZID')
newnames <- ifelse(is.na(newnames) | duplicated(newnames),names(newnames), newnames)
rownames(p$loadings) <- newnames

biplot(p, x="PC3", y="PC4")
biplot(p, x="PC5", y="PC6")

#normalization with ComBat_Seq
neurons_cts <-txi_griesi_neuron$counts

filter <- apply(neurons_cts, 1, function(x) length(x[x>10])>=12)
filtered <- neurons_cts[filter,]
neurons_cts2 <- apply(filtered, 2, function(x) as.integer(x))
rownames(neurons_cts2) <-rownames(filtered)

#normalise for neuron proportion
count_matrix <- as.matrix(neurons_cts2)
group <- c("F", "F", "E", "F", "F", "E", "D", "F", "D", "D", "F", "E", "F", "E", "D", "E", "D", "F", "F")
condition <- c(rep("ASD", 9), rep("CTL", 10))
adjusted_counts_neurons_griesi <- ComBat_seq(count_matrix, batch=group, group=condition)

dds_neurons1 <- DESeqDataSetFromMatrix(countData = adjusted_counts_neurons_griesi,
                                       colData = neurons_coldata,
                                       design = ~cell_type + group + condition)

dds_neurons1 <- estimateSizeFactors(dds_neurons1)
vsd_neurons <- vst(dds_neurons1,blind=TRUE)

pcaData <- plotPCA(vsd_neurons, intgroup = c("condition", "group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons after normalizing with ComBat")

#normalisation for library batch
batch <- neurons_coldata$cell_type
condition <- c(rep("ASD", 9), rep("CTL", 10))
adjusted_counts3 <- ComBat_seq(count_matrix, batch=batch, group=condition)

dds_neurons1 <- DESeqDataSetFromMatrix(countData = adjusted_counts3,
                                       colData = neurons_coldata,
                                       design = ~cell_type + group + condition)

dds_neurons1 <- estimateSizeFactors(dds_neurons1)
vsd_neurons <- vst(dds_neurons1,blind=TRUE)

pcaData <- plotPCA(vsd_neurons, intgroup = c("condition", "cell_type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =cell_type )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons after normalizing with ComBat")

#normalise for both
batch <- neurons_coldata$cell_type
condition <- c(rep("ASD", 9), rep("CTL", 10))
adjusted_counts_both <- ComBat_seq(adjusted_counts_neurons_griesi, batch=batch, group=condition)

dds_neurons <- DESeqDataSetFromMatrix(countData = adjusted_counts_both ,
                                      colData = neurons_coldata,
                                      design = ~cell_type + group + condition)

dds_neurons <- estimateSizeFactors(dds_neurons)
vsd_neurons <- vst(dds_neurons,blind=TRUE)

pcaData <- plotPCA(vsd_neurons, intgroup = c("condition", "group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons after normalizing with ComBat")


## RUVSEq normalization
#store the data in an object of S4 class SeqExpressionSet from the EDASeq package
b<- c(rep(c("ASD"), 9), rep(c("CTL"),10))
x <- as.factor(b)
set <- newSeqExpressionSet(as.matrix(neurons_cts2), phenoData = data.frame(x, row.names=colnames(neurons_cts2)))

#RLE = log-ratio of read count to median read count across sample
plotRLE(set, outline=FALSE, ylim=c(-4, 4))
plotPCA(set, cex=1)

#normalize the data using upper-quartile (UQ) normalization
set <- betweenLaneNormalization(set, which="upper")
plotPCA(set, cex=0.7)
plotRLE(set, outline=FALSE, ylim=c(-4, 4))

#obtain a set of “in-silico empirical” negative controls, 
#e.g., least significantly DE genes based on a first-pass DE analysis 
#performed prior to RUVg normalization.
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=4)
new_data <-pData(set2)
plotPCA(set2, cex=0.6)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4))

counts2neurons <-counts(set2)

#remove outlier
counts2neurons <-counts2neurons[,-9]

#get matrix for neurons
#remove outlier for k=4- SRR10775642
neurons_coldata<-subset(coldata2, subset = cell_type !="NPC")
neurons_coldata2 <-neurons_coldata[-9,]

normalized_counts_neurons <- normCounts(set2)
normalized_counts2 <-normalized_counts_neurons [,-9]


#pca on normalised counts
dds_count_norm <- DESeqDataSetFromMatrix(countData = normalized_counts_neurons,
                                         colData = neurons_coldata,
                                         design = ~cell_type + group + condition)

dds_count_norm <- estimateSizeFactors(dds_count_norm)
vsd_neurons_norm <- vst(dds_count_norm,blind=TRUE)

pcaData <- plotPCA(vsd_neurons_norm, intgroup = c("condition", "group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape =group.1 )) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data for neurons after normalizing with RUVSeq")

#############
isexpr = rowSums(cpm(counts2neurons)>0.1) >= 5

# Create DGEList object that stores counts
geneExpr <- DGEList( counts2neurons[isexpr,] )
geneExpr <- calcNormFactors( geneExpr )
colnames(geneExpr) <-rownames(neurons_coldata2)
#add df genes to geneExpr
geneid <-rownames(geneExpr)
genes <-AnnotationDbi::select(org.Hs.eg.db, keys = geneid, columns = c("SYMBOL", "ENSEMBL"), keytype = "ENTREZID")
genes <-genes[!duplicated(genes$ENTREZID),]
geneExpr$genes <-genes

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

#categorical variable to be modeled as a random effect
# The variable to be tested must be a fixed effect
new_data2<- new_data[-9,]
#use the coefficients produced by RUVSeq in the design
form <- ~x + W_1 + W_2 + W_3 + W_4 + (1|individual)
new_data2$individual <-neurons_coldata2$individual 

# estimate weights using linear mixed model of dream
#Transform count data to log2-counts per million (logCPM), 
#estimate the mean-variance relationship 
#and use this to compute appropriate observation-level weights. 
#The data are then ready for linear mixed modelling with dream(). 
#This method is the same as limma::voom(), except that it allows random effects in the formula
vobjDream <- voomWithDreamWeights(geneExpr, form, new_data2 )

# fit dream model with contrasts
#Fit linear mixed model for differential expression 
#and preform hypothesis test on fixed effects as specified in the contrast matrix L
fit <- dream( vobjDream, form, new_data2 )
fitm = eBayes(fit)
top.table_neurons <- topTable( fitm, number= 1000, coef="xCTL")
View(top.table_neurons)

#analysis using normalised counts
isexpr = rowSums(cpm(normalized_counts2)>0.1) >= 5
geneExpr_norm <- DGEList(normalized_counts2[isexpr,] )
geneExpr_norm <- calcNormFactors( geneExpr_norm)
colnames(geneExpr_norm) <-rownames(neurons_coldata2)
geneid <-rownames(geneExpr_norm)
genes <-AnnotationDbi::select(org.Hs.eg.db, keys = geneid, columns = c("SYMBOL", "ENSEMBL"), keytype = "ENTREZID")
genes <-genes[!duplicated(genes$ENTREZID),]
geneExpr_norm$genes <-genes


form2 <- ~condition + (1|individual)
vobjDream_neurons = voomWithDreamWeights(geneExpr_norm, form2, neurons_coldata2 )
fit_neurons = dream( vobjDream_neurons, form2, neurons_coldata2)
fitm_neurons =eBayes(fit_neurons)
top.table_neurons1 <- topTable( fitm_neurons, number= 1000, coef="conditionCTL")
top.table_neurons05 <- subset(top.table_neurons1, adj.P.Val <0.05)
top.table_neurons05_ordered<- top.table_neurons05[order(abs(top.table_neurons05$logFC), decreasing=TRUE),]
View(top.table_neurons05_ordered)

write.csv(top.table_neurons05_ordered,"DE_neurons_griesi_RUVSeq.csv",row.names=FALSE)
writeLines(top.table_neurons05_ordered$ENTREZID, "DElist_griesi_neurons.txt")

# list of genes neurons, k= 4, outlier removed, use normalized counts
list_DE_neurons_norm_K4 <- top.table_neurons05_ordered$SYMBOL
list_DE_neurons_norm_K4
# [1] "LOC112268044" "PCDHGB6"      "RNF139-DT"    "LYNX1-SLURP2" "ACSL6"        "TRIM7"        "PDE1A"        "NUP88"       
#[9] "MANEAL"       "INSRR"        "NME7"         "ANKRD20A1"    "LOC112268030" "SEMA3C"       "SPDYE2"       "DNPEP"       
#[17] "TMEM59L"      "LINC00402"    "MEIS1"        "KLKB1"        "ZNF467"       "LBX1-AS1"     "AVPR1A"       "SFMBT2"      
#[25] "ARHGEF16"     "MAPT-AS1"     "DACH2"        "CTSH"         "EPHX2"        "AK2"          "GNG7"         "TTYH2"       
#[33] "ABCC8"        "MYH14"        "LINC01750"    "SLC39A1"      "CNTNAP5-DT"   "TNNI1"        "GAD1"         "KCNS1"       
#[41] "SMTNL2"       "PCSK4"        "GASK1A"       "HLA-DQB1"     "OSBPL11"      "ADAM9"        "GPC4"         "SSPN"        
#[49] "ZCCHC12"      "AP1B1"        "RGPD6"        "CD151"        "POU4F1"       "C2orf72"      "PTPRN2"       "CFAP210"     
#[57] "TBC1D3I"      "CTNNA2"       "VWA5B2"       "SLC26A11"     "SLC32A1"      "CFAP43"       "CFAP74"       "CACNA1I"     
#[65] "ZIM2-AS1"     "DLEC1"        "TFDP2"        "AFF1-AS1"     "IPO8"         "KCNA2"        "SHANK3"       "PPP1R17"     
#[73] "LEFTY2"       "TTC39C"       "SLC34A2"      "LZTS3"        "SLC43A2"      "BPNT1"        "THSD7A"       "PTH1R"       
#[81] "ERICH2-DT"    "CPS1"         "KCTD20"       "FBXO2"        "MOK"          "SPEF1"


DE_list_griesi_neurons_RUVSeq <-top.table_neurons05_ordered$ENTREZID
upregul_griesi_neurons_RUVSeq <-subset(top.table_neurons05_ordered, logFC>0)
downregul_griesi_neurons_RUVSeq <-subset(top.table_neurons05_ordered, logFC<0)

list_upregul_griesi_neurons_RUVSeq <-upregul_griesi_neurons_RUVSeq$ENTREZID
list_downregul_griesi_neurons_RUVSeq <-downregul_griesi_neurons_RUVSeq$ENTREZID

#GO category enrichment with clusterProfiler
ego2 <- clusterProfiler::enrichGO(gene          = list_downregul_griesi_neurons_RUVSeq,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)

#no enriched term found...

## GSEA
geneList <- top.table_neurons05_ordered$logFC
names(geneList) <- top.table_neurons05_ordered$ENTREZID
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

#gene set enrichment; determines whether a pre-defined set of genes 
#shows statistically significant, concordant differences between two biological states
gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPermSimple = 10000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

#no term enriched under specific pvalueCutoff...

#import spreadsheet with sfari genes
Sfari_list <- read.csv("~/Desktop/MRes/20220611/SFARI.csv", header =TRUE, stringsAsFactors = FALSE)

Sfari_high_score <- subset(Sfari_list, gene.score %in% c(1,2))

list_SFARI <- Sfari_high_score$gene.symbol
list_SFARI_all <- Sfari_list$gene.symbol

go.obj <-newGeneOverlap(list_SFARI, list_DE_neurons_norm_K4, 25659)
go.obj <-testGeneOverlap(go.obj)
go.obj
getIntersection((go.obj)) #"AVPR1A"   "CACNA1I"  "GPC4"     "MAPT-AS1" "SHANK3"

#GeneOverlap object:
#listA size=895
#listB size=86
#Intersection size=5
#Overlapping p-value=0.18
#Jaccard Index=0.0


#overlap downregul IFN neurons (from IFN DE analysis)- any DE ASD 
go.obj <-newGeneOverlap(downregul_neurons, list_DE_neurons_norm_K4, 25659)
go.obj <-testGeneOverlap(go.obj)
go.obj

getIntersection((go.obj)) #"SEMA3C" "INSRR""
#GeneOverlap object:
#listA size=250
#listB size=86
#Intersection size=2
#Overlapping p-value=0.2
#Jaccard Index=0.0

#overlap downregul IFN NPC - ASD 
go.obj <-newGeneOverlap(downregul_NPC, list_DE_neurons_norm_K4, 25659)
go.obj <-testGeneOverlap(go.obj)
go.obj
getIntersection((go.obj)) #"ACSL6"    "PPP1R17"  "ARHGEF16"

#GeneOverlap object:
#listA size=547
#listB size=86
#Intersection size=3
#Overlapping p-value=0.28
#Jaccard Index=0.0

#overlap with genes found DE by the authors
#import spreadsheet .
griesi_neurons_DE_authors <- read_excel("griesi_DE_neurons.xlsx", col_names =TRUE)
#convert to df
df_griesi_neurons_DE_authors <-as.data.frame(griesi_neurons_DE_authors, col.names = colnames(griesi_neurons_DE_authors), row.names = rownames(griesi_neurons_DE_authors))
colnames(griesi_neurons_DE_authors) <-griesi_neurons_DE_authors[2,]
griesi_neurons_DE_authors <-griesi_neurons_DE_authors[3:15028,]
griesi_neurons_DE_authors <-griesi_neurons_DE_authors[1:20,]
griesi_neurons_DE_authors_list <-griesi_neurons_DE_authors$Gene
griesi_neurons_DE_authors_list

go.obj <-newGeneOverlap(griesi_neurons_DE_authors_list,list_DE_neurons_norm_K4, 25659)
go.obj <-testGeneOverlap(go.obj) 
go.obj

#GeneOverlap object:
# listA size=20
#listB size=86
#Intersection size=7
#Overlapping p-value=2.8e-13
#Jaccard Index=0.1

getIntersection((go.obj))
#"NUP88"   "SEMA3C"  "NME7"    "ACSL6"   "CPS1"    "ZCCHC12" "C2orf72"

#intersection authors with IFN
go.obj <-newGeneOverlap(downregul_NPC, griesi_neurons_DE_authors_list, 20565)
go.obj <-testGeneOverlap(go.obj)
go.obj #0.098
getIntersection((go.obj)) #"ACSL6" "KCNB1"

#intersection authors with IFN
go.obj <-newGeneOverlap(downregul_neurons, griesi_neurons_DE_authors_list, 20890)
go.obj <-testGeneOverlap(go.obj)
go.obj #0.024
getIntersection((go.obj)) #""SEMA3C"  "SLC12A5"
