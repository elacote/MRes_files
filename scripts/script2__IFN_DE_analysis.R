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

#BiocManager::install("ggnewscale")
library("ggnewscale")

#BiocManager::install("readxl")
library("readxl")

#BiocManager::install("GeneOverlap")
library("GeneOverlap")

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

#BiocManager::install("variancePartition")
library('variancePartition')

#BiocManager::install("edgeR")
library('edgeR')

#BiocManager::install("BiocParallel")
library('BiocParallel')

########################
#directory with the files containing the abundance data
dir <- "~/Desktop/MRes/20210402/IFN_from_kallisto/"
list.files(dir)

#names of the folders containing the abundance data
samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)

#list of files to be used to get abundance data for all the samples
files <- file.path(dir,  samples$V1, "abundance.tsv")
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
#select NPCs only
txi_NPC <- lapply(txi, function(x) if(is.matrix(x)) return(x[,1:6]) else return(x))

#import matrix with experiment data
deepak_coldata <- readRDS("~/Desktop/MRes/20210214/coldata1.rds")
#select NPCs only and modify the matrix
NPC_coldata <-subset(deepak_coldata, subset = cell_type =="NPC")
NPC_coldata$D30 <-NULL
NPC_coldata$cell_type <-NULL
rownames(NPC_coldata) <-colnames(txi_NPC)

#convert the data as factors
NPC_coldata$cell_line <- factor(NPC_coldata$cell_line, levels = c("M1","M2", "M3"))
NPC_coldata$D18 <- factor(NPC_coldata$D18, levels = c("untreated","treated"))

#construct a DESeqDataSet 
dds_NPC <- DESeqDataSetFromTximport(txi_NPC,
                                    colData = NPC_coldata,
                                    design = ~  cell_line + D18)

#prefiltering: only keep genes with at least 10 reads aligned
keep <- rowSums(counts(dds_NPC)) >= 10
table(keep)
#FALSE  TRUE 
#10081 20565 
dds_NPC <- dds_NPC[keep,]
dds_NPC <- estimateSizeFactors(dds_NPC)

#Differential expression
dds_NPC <-DESeq(dds_NPC)
#set argument alpha for FDR cutoff to 0.05
res05_NPC <- results(dds_NPC, alpha=0.05)

tot_genes_NPC <- as.data.frame(results(dds_NPC)) #for overlap with literature

#obtain a list of genes with significant DE, padj <0.05
DE_list_NPC <- subset(res05_NPC, padj <0.05)
# Add gene symbol column
DE_list_NPC$SYMBOL = mapIds(org.Hs.eg.db,
                            keys=rownames(DE_list_NPC), 
                            column="SYMBOL",
                            keytype="ENTREZID",
                            multiVals="first")

#list of non DE genes for IsoformSwitchAnalyzeR analysis
non_DE_list_NPC <-subset(res05_NPC, padj >0.05)
non_DE_list_NPC$SYMBOL = mapIds(org.Hs.eg.db,
                            keys=rownames(non_DE_list_NPC), 
                            column="SYMBOL",
                            keytype="ENTREZID",
                            multiVals="first")

#order DE genes by Log2Fold Change
DE_list_NPC_IFN <-DE_list_NPC[order(abs(DE_list_NPC$log2FoldChange), decreasing=TRUE),]

#genes overexpressed
DE_list_NPC_IFN_over <- subset(DE_list_NPC_IFN, log2FoldChange > 0)
df_DE_list_NPC_IFN_over <- as.data.frame(DE_list_NPC_IFN_over)
rownames(df_DE_list_NPC_IFN_over) <- df_DE_list_NPC_IFN_over$SYMBOL
df_DE_list_NPC_IFN_over2 <- format(df_DE_list_NPC_IFN_over, digits=2)
View(df_DE_list_NPC_IFN_over2[1:20,c(1,2,6)])
length(DE_list_NPC_IFN_over$SYMBOL) #1154

#genes underexpressed
DE_list_NPC_IFN_under <- subset(DE_list_NPC_IFN, log2FoldChange < 0)
df_DE_list_NPC_IFN_under <- as.data.frame(DE_list_NPC_IFN_under)
rownames(df_DE_list_NPC_IFN_under) <- df_DE_list_NPC_IFN_under$SYMBOL
df_DE_list_NPC_IFN_under2 <- format(df_DE_list_NPC_IFN_under, digits=2)
View(df_DE_list_NPC_IFN_under2[1:20,c(1,2,6)])
length(DE_list_NPC_IFN_under$SYMBOL) #547

#convert DE info into a dataframe
df_DE_list_NPC_IFN <- as.data.frame(DE_list_NPC_IFN)
View(df_DE_list_NPC_IFN)

write.csv(df_DE_list_NPC_IFN, file= "DE_List_NPC_IFN20220702.csv", row.names = T)

#to assess GO enrichment with DAVID tool
writeLines(rownames(df_DE_list_NPC_IFN), "DE_list_IFN_NPC20220702.txt")
writeLines(rownames(DE_list_NPC_IFN_under ), "DE_list_down__IFN_NPC.txt")
writeLines(rownames(DE_list_NPC_IFN_over), "DE_list_up__IFN_NPC.txt")

#overlap with DE obtained by the authors of the publication
DGE_IFN_lit_NPC <- read_excel("~/Desktop/MRes/20220629/IFN_DE_authors.xlsx", sheet ="NPC", col_names =TRUE)
df_DGE_IFN_lit_NPC <-as.data.frame(DGE_IFN_lit_NPC, col.names = colnames(DGE_IFN_lit_NPC))
View(df_DGE_IFN_lit_NPC)
df_DGE_IFN_lit_NPC2 <- df_DGE_IFN_lit_NPC[2:1835, c(2,4,8)]
colnames(df_DGE_IFN_lit_NPC2) <- c("gene_name", "log2FC", "adjusted_p_value")
df_DGE_IFN_lit_NPC2$log2FC <- as.numeric(df_DGE_IFN_lit_NPC2$log2FC)
df_DGE_IFN_lit_NPC3 <-df_DGE_IFN_lit_NPC2[order(abs(df_DGE_IFN_lit_NPC2$log2FC), decreasing=TRUE),]
View(df_DGE_IFN_lit_NPC3[1:500,])

go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL[1:100], df_DGE_IFN_lit_NPC3$gene_name[1:100], length(rownames(res05_NPC))  )
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap_NPC <- getIntersection(go.obj)

#look for a particular gene in the DE list
gene_studied <-"RELN"
line_gene <-grep(gene_studied,DE_list_NPC_IFN$SYMBOL)
line_gene

#volcano plot
EnhancedVolcano(DE_list_NPC_IFN,
                lab = DE_list_NPC_IFN$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'NPCs  (18U vs 18T)',
                pCutoff = 10e-32,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)

# GO analysis with clusterprofiler
OrgDb <- org.Hs.eg.db
# enrichment in term of BP of the set of significantly DE genes
ego2_NPC <- clusterProfiler::enrichGO(gene          = rownames(DE_list_NPC_IFN),
                                      OrgDb         = OrgDb,
                                      ont           = "BP",
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.01, 
                                      readable      = TRUE)

ego2_NPC_1 <- clusterProfiler::enrichGO(gene          = rownames(DE_list_NPC_IFN_over),
                                      OrgDb         = OrgDb,
                                      ont           = "BP",
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.01, 
                                      readable      = TRUE)

ego2_NPC_2 <- clusterProfiler::enrichGO(gene          = rownames(DE_list_NPC_IFN_under),
                                        OrgDb         = OrgDb,
                                        ont           = "BP",
                                        pAdjustMethod = "fdr",
                                        pvalueCutoff  = 0.05,
                                        qvalueCutoff  = 0.01, 
                                        readable      = TRUE)


x2_NPC <- pairwise_termsim(ego2_NPC)
emapplot(x2_NPC, showCategory = 15, node_label = "category")

x2_NPC_over <- pairwise_termsim(ego2_NPC_1)
emapplot(x2_NPC_over, showCategory = 15, node_label = "category")
x2_NPC_under <- pairwise_termsim(ego2_NPC_2)
emapplot(x2_NPC_under, showCategory = 15, node_label = "category")

## GSEA
#Gene list with Fold Change data
geneList_NPC <- DE_list_NPC_IFN$log2FoldChange
names(geneList_NPC) <- rownames(DE_list_NPC_IFN)
geneList_NPC <- na.omit(geneList_NPC)
geneList_NPC <- sort(geneList_NPC, decreasing = TRUE)

#gene set enrichment; determines whether a pre-defined set of genes 
#shows statistically significant, concordant differences between two biological states
gse <- gseGO(geneList=geneList_NPC, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPermSimple = 10000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=20, split=".sign", font.size = 6) + facet_grid(.~.sign)

# analysis on neurons
txi_neurons <- lapply(txi, function(x) if(is.matrix(x)) return(x[,7:18]) else return(x))

neurons_coldata <-subset(deepak_coldata, subset = cell_type !="NPC")
neurons_coldata$cell_type <-NULL

# test the effect of the late treatment
neurons_UUvsUT <- subset(neurons_coldata, subset = D18 =="untreated")

neurons_UUvsUT$D18 <-NULL
neurons_UUvsUT$D30 <- factor(neurons_UUvsUT$D30, levels = c("untreated", "treated"))
neurons_UUvsUT$cell_line <-factor(neurons_UUvsUT$cell_line, levels =c("M1", "M2", "M3"))

txi_neurons_late <- lapply(txi, function(x) if(is.matrix(x)) return(x[,13:18]) else return(x))

rownames(neurons_UUvsUT) <- c("sample13", "sample14", "sample15", "sample16", "sample17", "sample18")

dds_late <- DESeqDataSetFromTximport(txi_neurons_late,
                                     colData = neurons_UUvsUT,
                                     design = ~ cell_line + D30)

#prefiltering: only keep genes with at least 10 reads aligned
keep <- rowSums(counts(dds_late)) >= 10
table(keep)
#late: 20890
dds_late<- dds_late[keep,]
dds_late<- estimateSizeFactors(dds_late)

#Differential expression
dds_late<- DESeq(dds_late)
#set argument alpha for FDR cutoff to 0.05
res05_late <- results(dds_late, alpha=0.05)
tot_genes_late <- as.data.frame(results(dds_late)) # for overlap with literature

DE_list_late_IFN <- subset(res05_late , padj <0.05)

#non DE list for IsoformSwitchAnalyzeR
nonDE_list_late_IFN <- subset(res05_late , padj >0.05)
nonDE_list_late_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                 keys=rownames(nonDE_list_late_IFN), 
                                 column="SYMBOL",
                                 keytype="ENTREZID",
                                 multiVals="first")

# Add gene symbol column to DE matrix
DE_list_late_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                 keys=rownames(DE_list_late_IFN), 
                                 column="SYMBOL",
                                 keytype="ENTREZID",
                                 multiVals="first")

DE_list_late_IFN <-DE_list_late_IFN[order(abs(DE_list_late_IFN$log2FoldChange), decreasing=TRUE),]
df_DE_list_late_IFN <-as.data.frame(DE_list_late_IFN )


#genes overexpressed
DE_list_over_late <- subset(DE_list_late_IFN, log2FoldChange > 0)
length(DE_list_over_late$SYMBOL) #late: 519
df_DE_list_over_late <-as.data.frame(DE_list_over_late )
rownames(df_DE_list_over_late ) <- df_DE_list_over_late$SYMBOL
df_DE_list_over_late2 <- format(df_DE_list_over_late, digits=2)
View(df_DE_list_over_late2[1:20,c(1,2,6)])

#genes underexpressed
DE_list_under_late <- subset(DE_list_late_IFN , log2FoldChange < 0)
length(DE_list_under_late$SYMBOL) #late: 205
df_DE_list_under_late <-as.data.frame(DE_list_under_late )
rownames(df_DE_list_under_late ) <- df_DE_list_under_late$SYMBOL
df_DE_list_under_late2 <- format(df_DE_list_under_late, digits=2)
View(df_DE_list_under_late2[1:20,c(1,2,6)])

#import spreadsheet with DGE from the authors
DGE_IFN_lit_late <- read_excel("~/Desktop/MRes/20220629/IFN_DE_authors.xlsx", sheet ="late", col_names =TRUE)
#convert to df
df_DGE_IFN_lit_late <-as.data.frame(DGE_IFN_lit_late, col.names = colnames(DGE_IFN_lit_late))
df_DGE_IFN_lit_late3 <-df_DGE_IFN_lit_late[order(abs(df_DGE_IFN_lit_late$log2FoldChange), decreasing=TRUE),]
View(df_DGE_IFN_lit_late3[1:500,])

go.obj <-newGeneOverlap(DE_list_late_IFN$SYMBOL[1:100], df_DGE_IFN_lit_late3$`Gene name`[1:100], length(rownames(res05_late))  )
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap <- getIntersection(go.obj)

#look for a particular gene in the DE list
gene_studied <-"RELN"
line_gene <-grep(gene_studied,df_DE_list_late_IFN$SYMBOL)
line_gene

write.csv(df_DE_list_late_IFN, file= "DE_List_Nneurons_late_IFN20220702.csv", row.names = T)

#to assess GO enrichment with DAVID tool
writeLines(rownames(DE_list_late_IFN), "DE_list_IFN_neuronslate20220702.txt")
writeLines(rownames(DE_list_under_late), "DE_list_down__IFN_neuronslate.txt")
writeLines(rownames(DE_list_over_late), "DE_list_up__IFN_neuronslate.txt")

#volcano plot
EnhancedVolcano(DE_list_late_IFN,
                lab = DE_list_late_IFN$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'neurons, acute effect of IFN-gamma treatment',
                pCutoff = 10e-32,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)

# GO enrichment with clusterprofiler
OrgDb <- org.Hs.eg.db
ego2 <- clusterProfiler::enrichGO(gene          = rownames(df_DE_list_late_IFN),
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2_late <- pairwise_termsim(ego2)
emapplot(x2_late, showCategory = 15, node_label = "category")

ego2_over <- clusterProfiler::enrichGO(gene          = rownames(DE_list_over_late),
                                       OrgDb         = OrgDb,
                                       ont           = "BP",
                                       pAdjustMethod = "fdr",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.01, 
                                       readable      = TRUE)
x2_late_over <- pairwise_termsim(ego2_over)
emapplot(x2_late_over, showCategory = 15, node_label = "category")

ego2_under <- clusterProfiler::enrichGO(gene          = rownames(DE_list_under_late),
                                        OrgDb         = OrgDb,
                                        ont           = "BP",
                                        pAdjustMethod = "fdr",
                                        pvalueCutoff  = 0.05,
                                        qvalueCutoff  = 0.01, 
                                        readable      = TRUE)
x2_late_under <- pairwise_termsim(ego2_under)
emapplot(x2_late_under, showCategory = 15, node_label = "category")

## GSEA
geneList <- df_DE_list_late_IFN$log2FoldChange
names(geneList) <- rownames(df_DE_list_late_IFN)
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

#gene set enrichment; determines whether a pre-defined set of genes 
#shows statistically significant, concordant differences between two biological states
gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPermSimple = 100000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=10, split=".sign", font.size = 6) + facet_grid(.~.sign)


# to test effect of the early treatment
neurons_UUvsTU <- subset(neurons_coldata, subset = D30 =="untreated")

neurons_UUvsTU$D30 <-NULL
neurons_UUvsTU$D18 <- factor(neurons_UUvsTU$D18, levels = c("untreated", "treated"))
neurons_UUvsTU$cell_line <-factor(neurons_UUvsTU$cell_line, levels =c("M1", "M2", "M3"))

rownames(neurons_UUvsTU) <- c("sample10", "sample11", "sample12", "sample16", "sample17", "sample18")
txi_neurons_early <- lapply(txi, function(x) if(is.matrix(x)) return(x[,c(10,11,12,16,17,18)]) else return(x))

dds_early <- DESeqDataSetFromTximport(txi_neurons_early,
                                      colData = neurons_UUvsTU,
                                      design = ~ cell_line + D18)


#prefiltering: only keep genes with at least 10 reads aligned
keep <- rowSums(counts(dds_early)) >= 10
table(keep)
#early: 21015
dds_early<- dds_early[keep,]
dds_early<- estimateSizeFactors(dds_early)

#Differential expression
dds_early<- DESeq(dds_early)
#set argument alpha for FDR cutoff to 0.05
res05_early <- results(dds_early, alpha=0.05)

tot_genes_early <- as.data.frame(results(dds_early)) # for overlap literature

DE_list_early_IFN <- subset(res05_early, padj <0.05)
DE_list_early_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                 keys=rownames(DE_list_early_IFN), 
                                 column="SYMBOL",
                                 keytype="ENTREZID",
                                 multiVals="first")

DE_list_early_IFN <-DE_list_early_IFN[order(abs(DE_list_early_IFN$log2FoldChange), decreasing=TRUE),]
df_DE_list_early_IFN <-as.data.frame(DE_list_early_IFN )

#non DE list for DTU analysis
nonDE_list_early_IFN <- subset(res05_early, padj >0.05)
nonDE_list_early_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                     keys=rownames(nonDE_list_early_IFN), 
                                     column="SYMBOL",
                                     keytype="ENTREZID",
                                     multiVals="first")

#genes overexpressed
DE_list_over_early <- subset(DE_list_early_IFN, log2FoldChange > 0)
length(DE_list_over_early$SYMBOL) #early: 60
df_DE_list_over_early <- as.data.frame(DE_list_over_early)
rownames(df_DE_list_over_early) <- df_DE_list_over_early$SYMBOL
df_DE_list_over_early2 <- format(df_DE_list_over_early, digits=2)
View(df_DE_list_over_early2[1:20,c(1,2,6)])

#genes underexpressed
DE_list_under_early<- subset(DE_list_early_IFN , log2FoldChange < 0)
length(DE_list_under_early$SYMBOL) #early: 35
df_DE_list_under_early <- as.data.frame(DE_list_under_early)
rownames(df_DE_list_under_early) <- df_DE_list_under_early$SYMBOL
df_DE_list_under_early2 <- format(df_DE_list_under_early, digits=2)
View(df_DE_list_under_early2[1:20,c(1,2,6)])


#import spreadsheet with DGE from the authors
DGE_IFN_lit_early <- read_excel("~/Desktop/MRes/20220629/IFN_DE_authors.xlsx", sheet ="early", col_names =TRUE)
df_DGE_IFN_lit_early <-as.data.frame(DGE_IFN_lit_early, col.names = colnames(DGE_IFN_lit_early))
df_DGE_IFN_lit_early2 <-df_DGE_IFN_lit_early[order(abs(df_DGE_IFN_lit_early$log2FoldChange), decreasing=TRUE),]
View(df_DGE_IFN_lit_early2[1:500,])

go.obj <-newGeneOverlap(DE_list_early_IFN$SYMBOL[1:25], df_DGE_IFN_lit_early2$`Gene name`[1:25], length(rownames(res05_early))  )
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap <- getIntersection(go.obj)

#look for a particular gene in the DE list
gene_studied <-"LOC112268044"
line_gene <-grep(gene_studied,df_DE_list_early_IFN$SYMBOL)
line_gene

write.csv(df_DE_list_early_IFN, file= "DE_List_Nneurons_early_IFN20220702.csv", row.names = T)

#to assess GO enrichment with DAVID tool
writeLines(rownames(DE_list_early_IFN), "DE_list_IFN_neuronsearly20220702.txt")
writeLines(rownames(DE_list_under_early), "DE_list_down__IFN_neuronsearly.txt")
writeLines(rownames(DE_list_over_early), "DE_list_up__IFN_neuronsearly.txt")

#volcano plot
EnhancedVolcano(DE_list_early_IFN,
                lab = DE_list_early_IFN$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'neurons, early effect of IFN-gamma treatment',
                pCutoff = 10e-3,
                FCcutoff = 0,
                pointSize = 3.0,
                labSize = 6.0)

# GO enrichment with clusterprofiler
OrgDb <- org.Hs.eg.db
# enrichment in term of BP of the set of significantly DE genes
ego2 <- clusterProfiler::enrichGO(gene          = rownames(df_DE_list_early_IFN),
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2_early <- pairwise_termsim(ego2)
emapplot(x2_early, showCategory = 15, node_label = "category")

ego2_over <- clusterProfiler::enrichGO(gene          = rownames(DE_list_over_early),
                                       OrgDb         = OrgDb,
                                       ont           = "BP",
                                       pAdjustMethod = "fdr",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.01, 
                                       readable      = TRUE)

x2_early_over <- pairwise_termsim(ego2_over)
emapplot(x2_early_over, showCategory = 15, node_label = "category")

ego2_under <- clusterProfiler::enrichGO(gene          = rownames(DE_list_under_early),
                                        OrgDb         = OrgDb,
                                        ont           = "BP",
                                        pAdjustMethod = "fdr",
                                        pvalueCutoff  = 0.05,
                                        qvalueCutoff  = 0.01, 
                                        readable      = TRUE)
x2_early_under <- pairwise_termsim(ego2_under)
emapplot(x2_early_under, showCategory = 15, node_label = "category")

## GSEA
geneList <- df_DE_list_early_IFN$log2FoldChange
names(geneList) <- rownames(df_DE_list_early_IFN)
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

#gene set enrichment; determines whether a pre-defined set of genes 
#shows statistically significant, concordant differences between two biological states
gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPermSimple = 100000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=10, split=".sign", font.size = 6) + facet_grid(.~.sign)


# to test effect of repeated treatment vs late
neurons_UTvsTT <- subset(neurons_coldata, subset = D30 =="treated")
neurons_UTvsTT$D18 <- factor(neurons_UTvsTT$D18, levels = c("untreated", "treated"))
neurons_UTvsTT$cell_line <-factor(neurons_UTvsTT$cell_line, levels =c("M1", "M2", "M3"))

txi_neurons_repeated <- lapply(txi, function(x) if(is.matrix(x)) return(x[,c(7,8,9,13,14,15)]) else return(x))

rownames(neurons_UTvsTT) <- c("sample7", "sample8", "sample9", "sample13", "sample14", "sample15")

dds_repvslate <- DESeqDataSetFromTximport(txi_neurons_repeated,
                                    colData = neurons_UTvsTT,
                                    design = ~ cell_line + D18)

#prefiltering: only keep genes with at least 10 reads aligned
keep <- rowSums(counts(dds_repvslate)) >= 10
table(keep)
#rep vs late 21153
dds_repvslate<- dds_repvslate[keep,]
dds_repvslate <- estimateSizeFactors(dds_repvslate)

#Differential expression
dds_repvslate <- DESeq(dds_repvslate)
#set argument alpha for FDR cutoff to 0.05
res05_repvslate <- results(dds_repvslate, alpha=0.05)

tot_genes_rep<- as.data.frame(results(dds_repvslate)) # for overlap literature

DE_list_repvslate_IFN <- subset(res05_repvslate, padj <0.05)
# Add gene symbol column
DE_list_repvslate_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                 keys=rownames(DE_list_repvslate_IFN), 
                                 column="SYMBOL",
                                 keytype="ENTREZID",
                                 multiVals="first")

DE_list_repvslate_IFN <-DE_list_repvslate_IFN[order(abs(DE_list_repvslate_IFN$log2FoldChange), decreasing=TRUE),]
df_DE_list_repvslate_IFN <-as.data.frame(DE_list_repvslate_IFN )

#non DE list for DTU analysis
nonDE_list_repvslate_IFN <- subset(res05_repvslate, padj >0.05)
nonDE_list_repvslate_IFN$SYMBOL = mapIds(org.Hs.eg.db,
                                         keys=rownames(nonDE_list_repvslate_IFN), 
                                         column="SYMBOL",
                                         keytype="ENTREZID",
                                         multiVals="first")

#genes overexpressed
DE_list_over_repvslate <- subset(DE_list_repvslate_IFN, log2FoldChange > 0)
length(DE_list_over_repvslate$SYMBOL)  #147
df_DE_list_over_repvslate <-as.data.frame(DE_list_over_repvslate )
rownames(df_DE_list_over_repvslate) <- df_DE_list_over_repvslate$SYMBOL
df_DE_list_over_repvslate2 <- format(df_DE_list_over_repvslate, digits=2)
View(df_DE_list_over_repvslate2[1:20,c(1,2,6)])

#genes underexpressed
DE_list_under_repvslate <- subset(DE_list_repvslate_IFN , log2FoldChange < 0)
length(DE_list_under_repvslate$SYMBOL) # 44
df_DE_list_under_repvslate <-as.data.frame(DE_list_under_repvslate )
rownames(df_DE_list_under_repvslate) <- df_DE_list_under_repvslate$SYMBOL
df_DE_list_under_repvslate2 <- format(df_DE_list_under_repvslate, digits=2)
View(df_DE_list_under_repvslate2[1:20,c(1,2,6)])


#look for a particular gene in the DE list
gene_studied <-"LOC112268044"
line_gene <-grep(gene_studied,df_DE_list_repvslate_IFN$SYMBOL)
line_gene

write.csv(df_DE_list_repvslate_IFN, file= "DE_List_Nneurons_repvslate_IFN20220702.csv", row.names = T)

#to assess GO enrichment with DAVID tool
writeLines(rownames(DE_list_repvslate_IFN), "DE_list_IFN_neuronsrepvslate20220702.txt")
writeLines(rownames(DE_list_under_repvslate), "DE_list_down__IFN_neuronsrepvslate.txt")
writeLines(rownames(DE_list_over_repvslate), "DE_list_up__IFN_neuronsrepvslate.txt")

# Volcano plot
EnhancedVolcano(DE_list_repvslate_IFN,
                lab = DE_list_repvslate_IFN$SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'neurons, repeated vs late effect of IFN-gamma treatment',
                pCutoff = 10e-3,
                FCcutoff = 0,
                pointSize = 3.0,
                labSize = 6.0)

# GO enrichment with clusterprofiler
OrgDb <- org.Hs.eg.db
# enrichment in term of BP of the set of significantly DE genes
ego2 <- clusterProfiler::enrichGO(gene          = rownames(df_DE_list_repvslate_IFN),
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2_rep <- pairwise_termsim(ego2)
emapplot(x2_rep, showCategory = 15, node_label = "category")

ego2_over <- clusterProfiler::enrichGO(gene          = rownames(DE_list_over_repvslate),
                                       OrgDb         = OrgDb,
                                       ont           = "BP",
                                       pAdjustMethod = "fdr",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.01, 
                                       readable      = TRUE)
x2_rep_over <- pairwise_termsim(ego2_over)
emapplot(x2_rep_over , showCategory = 15, node_label = "category")

ego2_under <- clusterProfiler::enrichGO(gene          = rownames(DE_list_under_repvslate),
                                        OrgDb         = OrgDb,
                                        ont           = "BP",
                                        pAdjustMethod = "fdr",
                                        pvalueCutoff  = 0.05,
                                        qvalueCutoff  = 0.01, 
                                        readable      = TRUE)
x2_rep_under <- pairwise_termsim(ego2_under)
emapplot(x2_rep_under , showCategory = 15, node_label = "category")

## GSEA
geneList <- df_DE_list_repvslate_IFN$log2FoldChange
names(geneList) <- rownames(df_DE_list_repvslate_IFN)
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

dotplot(gse, showCategory=10, split=".sign", font.size = 6) + facet_grid(.~.sign)

#overlaps
#number of genes filtered both for neurons and NPCs
dfNPC <- as.data.frame(results(dds_NPC))

dflate <- as.data.frame(results(dds_late))

dfearly <- as.data.frame(results(dds_early))

dfrepvslate <- as.data.frame(results(dds_repvslate))

#overlap NPC- early
#number of genes to consider for the overlap
genes <-NULL
for (name in rownames(dfNPC) ) {
  if (name %in% rownames(dfearly)){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #19572

go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL, DE_list_early_IFN$SYMBOL, 19572 )
go.obj <-testGeneOverlap(go.obj)
go.obj
overlap_NPC_early <- getIntersection(go.obj)
#GeneOverlap object:
#listA size=1701
#listB size=95
#Intersection size=42
#Overlapping p-value=3e-20
#Jaccard Index=0.0
writeLines(overlap_NPC_early, "overlap_IFN_DE_NPC_early20220702.txt")

DGE_IFNg_NPCs_early <- list(
  NPCs<- DE_list_NPC_IFN$SYMBOL,
  early <-DE_list_early_IFN$SYMBOL
)

p <- ggVennDiagram(DGE_IFNg_NPCs_early, category.names = c("NPCs","neurons early treatment"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)

#overlap between early and late
#number of genes filtered both for neurons and NPCs
genes <-NULL
for (name in rownames(dflate) ) {
  if (name %in% rownames(dfearly)){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #20576
go.obj <-newGeneOverlap(DE_list_late_IFN$SYMBOL, DE_list_early_IFN$SYMBOL, 20576 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
# listA size=724
#listB size=95
#Intersection size=43
#Overlapping p-value=3.5e-37
#Jaccard Index=0.1

overlap_early_late <- getIntersection(go.obj)
writeLines(overlap_early_late, "overlap_IFN_DE_late_early20220702.txt") 

DGE_IFNg_early_late <- list(
  late<- DE_list_late_IFN$SYMBOL,
  early <-DE_list_early_IFN$SYMBOL
)

p <- ggVennDiagram(DGE_IFNg_early_late, category.names = c("acute treatment","early treatment"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)

#overlap between NPCs and late
#number of genes filtered both for neurons and NPCs
genes <-NULL
for (name in rownames(dflate) ) {
  if (name %in% rownames(dfNPC)){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #19541

go.obj <-newGeneOverlap(DE_list_late_IFN$SYMBOL,DE_list_NPC_IFN$SYMBOL, 19541 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
# listA size=724
#listB size=1701
#Intersection size=398
#Overlapping p-value=4.5e-238
#Jaccard Index=0.2

overlap_NPC_late <- getIntersection(go.obj)
writeLines(overlap_NPC_late, "overlap_IFN_DE_late_NPC20220702.txt")

DE_NPC_late <- subset(DE_list_late_IFN, DE_list_late_IFN$SYMBOL %in% overlap_NPC_late)
DE_NPC_early <- subset(DE_list_early_IFN, DE_list_early_IFN$SYMBOL %in%  overlap_NPC_early)
DE_early_late <-subset(DE_list_late_IFN, DE_list_late_IFN$SYMBOL %in% overlap_early_late)

list_for_enrichment <-rownames(DE_early_late)
list_for_enrichment1 <-rownames(DE_NPC_late)
list_for_enrichment2 <-rownames(DE_NPC_early)
writeLines(list_for_enrichment, "overlap_IFN_DE_late_early_E_20220702.txt")
writeLines(list_for_enrichment1, "overlap_IFN_DE_NPC_late_E_20220702.txt")
writeLines(list_for_enrichment2, "overlap_IFN_DE_NPC_early_E_20220702.txt")

#enrichment of overlap early-late
ego2 <- clusterProfiler::enrichGO(gene          = list_for_enrichment,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, node_label = "category")

#enrichment of overlap NPC-late
ego2 <- clusterProfiler::enrichGO(gene          = list_for_enrichment1,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, node_label = "category")

#enrichment of overlap NPC-early
ego2 <- clusterProfiler::enrichGO(gene          = list_for_enrichment2,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, node_label = "category")

#Venn diagram
DGE_IFNg_NPC_late <- list(
  NPCs<- DE_list_NPC_IFN$SYMBOL,
  late <-DE_list_late_IFN$SYMBOL
)

p <- ggVennDiagram(DGE_IFNg_NPC_late, category.names = c("NPCs","acute treatment"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)


#overlap between NPCs and SFARI
#import spreadsheet with sfari genes
table_Sfari <- read.csv("~/Desktop/MRes/20220611/SFARI.csv", header =TRUE, stringsAsFactors = FALSE)
Sfari_high_score <- subset(table_Sfari, gene.score %in% c(1,2))
list_SFARI <- Sfari_high_score$gene.symbol

go.obj <-newGeneOverlap(list_SFARI, DE_list_NPC_IFN$SYMBOL, 20068 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
# listA size=8955
#listB size=1701
#Intersection size=93
#Overlapping p-value=0.023
#Jaccard Index=0.0

overlap_NPC_sfari <- getIntersection(go.obj)
writeLines(overlap_NPC_sfari, "overlap_IFN_DE_NPC_sfari20220702.txt")

DE_NPC_Sfari <- subset(DE_list_NPC_IFN, DE_list_NPC_IFN$SYMBOL %in% overlap_NPC_sfari)
write.csv(DE_NPC_Sfari, file= "overlap_IFN_NPCDE_sfari_20220702.csv", row.names = T)
writeLines(rownames(DE_NPC_Sfari), "overlap_IFN_DE_NPC_sfari_E_20220702.txt")

#GO enrichment of the overlap genes DE by IFN-g- Sfari
OrgDb <- org.Hs.eg.db
# enrichment in term of BP of the set of significantly DE genes
ego2 <- clusterProfiler::enrichGO(gene          = rownames(DE_NPC_Sfari),
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)

x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

# overlap between sfari genes and genes DE by a late IFN-g treatment
go.obj <-newGeneOverlap(list_SFARI, DE_list_late_IFN$SYMBOL, 20068 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
 # listA size=895
#listB size=724
#Intersection size=58
#Overlapping p-value=1.2e-05
#Jaccard Index=0.0

overlap_late_sfari <- getIntersection(go.obj)
writeLines(overlap_late_sfari, "overlap_IFN_DE_late_sfari20220702.txt")

DE_late_Sfari <- subset(DE_list_late_IFN, DE_list_late_IFN$SYMBOL %in% overlap_late_sfari)
write.csv(DE_late_Sfari, file= "overlap_IFN_lateDE_sfari_20220702.csv", row.names = T)
writeLines(rownames(DE_late_Sfari), "overlap_IFN_DE_late_sfari_E_20220702.txt")

#GO enrichment of the overlap late-Sfari
ego2 <- clusterProfiler::enrichGO(gene          = rownames(DE_late_Sfari),
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

#Venn diagram
overlap_NPC_sfari <- list(
  NPCs<- DE_list_NPC_IFN$SYMBOL,
  late <-list_SFARI
)

p <- ggVennDiagram(overlap_NPC_sfari, category.names = c("IFN-gamma on NPCs","SFARI"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)

overlap_late_sfari <- list(
  late<- DE_list_late_IFN$SYMBOL,
  sfari <-list_SFARI
)

p <- ggVennDiagram(overlap_late_sfari, category.names = c("IFN-gamma acute on neurons","SFARI"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)

#list of downregulated genes in neurons (for overlap with WGCNA modules)
downregul_neurons <-c(DE_list_under_late$SYMBOL, DE_list_under_early$SYMBOL, DE_list_under_repvslate$SYMBOL)
downregul_neurons <-unique(downregul_neurons)
length(downregul_neurons) #250

#list of genes downregul in NPCs (for overlap with WGCNA modules)
downregul_NPC <- DE_list_NPC_IFN_under$SYMBOL
length(DE_list_NPC_IFN_under$SYMBOL) #547

#DE analysis with limma/voom
#select neurons only
txi_neurons <- lapply(txi, function(x) if(is.matrix(x)) return(x[,7:18]) else return(x))
coldata1 <- readRDS("~/Desktop/MRes/20210214/coldata1.rds")
coldata <-subset(coldata1, subset = cell_type !="NPC")
coldata$cell_type <-NULL
coldata$D18 <-NULL
coldata$D30 <-NULL

counts_from_kallisto <- txi_neurons$counts
colnames(counts_from_kallisto) <- rownames(coldata)

# filter genes by number of counts
isexpr = rowSums(counts_from_kallisto) >= 10

# Create DGEList object that stores counts
geneExpr = DGEList( counts_from_kallisto[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

#add df genes to geneExpr
geneid <-rownames(geneExpr)
genes <-AnnotationDbi::select(org.Hs.eg.db, keys = geneid, columns = c("SYMBOL", "ENSEMBL"), keytype = "ENTREZID")
genes <-genes[!duplicated(genes$ENTREZID),]
geneExpr$genes <-genes #22491


treatment <- factor(
  x= c(rep("TT",3), rep("TU", 3), rep("UT", 3), rep("UU", 3)),
  levels =c("UU", "UT", "TU", "TT")
)
cell_line <- factor(
  x= c("M1", "M2", "M3", "M1", "M2", "M3", "M1", "M2", "M3", "M1", "M2", "M3"),
  levels = c("M1", "M2", "M3")
)

#reference UU
coldata$treatment <-treatment
coldata$cell_line <-cell_line

design <- model.matrix(~treatment, coldata)
vobj_temp <- voom(geneExpr, design, plot=TRUE)
dupcor <-duplicateCorrelation(vobj_temp, design, block =coldata$cell_line)
vobj <- voom(geneExpr, design, plot=FALSE, block =coldata$cell_line, correlation = dupcor$consensus.correlation )
dupcor <-duplicateCorrelation(vobj, design, block =coldata$cell_line)
fitdupcor <-lmFit(vobj,design, block =coldata$cell_line, correlation = dupcor$consensus.correlation ) 
fitdupcor <-eBayes(fitdupcor)

head(fitdupcor$coefficients,5)

#late treatment
top.table_neurons_IFN_late_limma  <- topTable(
  fitdupcor,
  n=Inf,
  adjust.method = "BH",
  coef ='treatmentUT'
)
top.table_neurons_IFN_05_late_limma <- subset(top.table_neurons_IFN_late_limma, adj.P.Val <0.05 )
top.table_neurons_IFN_05_late_limma_ordered<- top.table_neurons_IFN_05_late_limma[order(abs(top.table_neurons_IFN_05_late_limma$logFC), decreasing=TRUE),]
View(top.table_neurons_IFN_05_late_limma_ordered)
write.csv(top.table_neurons_IFN_05_late_limma_ordered, file= "limma_DE_List_neurons_late_IFN20220702.csv", row.names = T)

non_DE_neurons_IFN_05_late_limma <- subset(top.table_neurons_IFN_late_limma, adj.P.Val >0.05 )

gene_studied <-"MOV10"
line_gene <-grep(gene_studied,top.table_neurons_IFN_05_late_limma_ordered$SYMBOL)
line_gene

#early treatment
top.table_neurons_IFN_early_limma  <- topTable(
  fitdupcor,
  n=Inf,
  adjust.method = "BH",
  coef ='treatmentTU'
)
top.table_neurons_IFN_05_early_limma <- subset(top.table_neurons_IFN_early_limma, adj.P.Val <0.05 )
top.table_neurons_IFN_05_early_limma_ordered<- top.table_neurons_IFN_05_early_limma[order(abs(top.table_neurons_IFN_05_early_limma$logFC), decreasing=TRUE),]
View(top.table_neurons_IFN_05_early_limma_ordered)
View(top.table_neurons_IFN_05_early_limma_ordered[, c(2,4,8)])

write.csv(top.table_neurons_IFN_05_early_limma_ordered, file= "limma_DE_List_neurons_early_IFN20220702.csv", row.names = T)

##
#repeated treatment
top.table_neurons_IFN_rep_limma  <- topTable(
  fitdupcor,
  n=Inf,
  adjust.method = "BH",
  coef ='treatmentTT'
)
top.table_neurons_IFN_05_rep_limma <- subset(top.table_neurons_IFN_rep_limma, adj.P.Val <0.05 )
top.table_neurons_IFN_05_rep_limma_ordered<- top.table_neurons_IFN_05_rep_limma[order(abs(top.table_neurons_IFN_05_rep_limma$logFC), decreasing=TRUE),]
View(top.table_neurons_IFN_05_rep_limma_ordered)
write.csv(top.table_neurons_IFN_05_rep_limma_ordered, file= "limma_DE_List_neurons_rep_IFN20220702.csv", row.names = T)

gene_studied <-"SORC3"
line_gene <-grep(gene_studied,top.table_neurons_IFN_05_early_limma_ordered$SYMBOL)
line_gene

#reference late for comparison repeated vs late
treatment2 <- factor(
  x= c(rep("TT",3), rep("TU", 3), rep("UT", 3), rep("UU", 3)),
  levels =c("UT", "UU", "TU", "TT")
)
coldata2 <-coldata
coldata2$treatment <-treatment2

design2 <- model.matrix(~treatment, coldata2)
vobj_temp2 <- voom(geneExpr, design2, plot=TRUE)
dupcor2 <-duplicateCorrelation(vobj_temp2, design2, block =coldata2$cell_line)
vobj2 <- voom(geneExpr, design2, plot=FALSE, block =coldata2$cell_line, correlation = dupcor2$consensus.correlation )
dupcor2 <-duplicateCorrelation(vobj2, design2, block =coldata2$cell_line)
fitdupcor2 <-lmFit(vobj2,design2, block =coldata2$cell_line, correlation = dupcor2$consensus.correlation ) 
fitdupcor2 <-eBayes(fitdupcor2)

head(fitdupcor2$coefficients,5)

#repeated vslate treatment
top.table_neurons_IFN_repvslate_limma  <- topTable(
  fitdupcor2,
  n=Inf,
  adjust.method = "BH",
  coef ='treatmentTT'
)
top.table_neurons_IFN_05_repvslate_limma <- subset(top.table_neurons_IFN_repvslate_limma, adj.P.Val <0.05 )
top.table_neurons_IFN_05_repvslate_limma_ordered<- top.table_neurons_IFN_05_repvslate_limma[order(abs(top.table_neurons_IFN_05_repvslate_limma$logFC), decreasing=TRUE),]
View(top.table_neurons_IFN_05_repvslate_limma_ordered)
View(top.table_neurons_IFN_05_repvslate_limma_ordered[, c(2,4,8)])
write.csv(top.table_neurons_IFN_05_repvslate_limma_ordered, file= "limma_DE_List_neurons_repvslate_IFN20220702.csv", row.names = T)

gene_studied <-"MOV10"
line_gene <-grep(gene_studied,top.table_neurons_IFN_05_repvslate_limma_ordered$SYMBOL)
line_gene


#Go enrichment
late_upregul <- subset(top.table_neurons_IFN_05_late_limma_ordered, logFC >0)
late_downregul <- subset(top.table_neurons_IFN_05_late_limma_ordered, logFC < 0)
View(late_upregul[1:20, c(2,4,8)])
late_downregul <- subset(top.table_neurons_IFN_05_late_limma_ordered, logFC < 0)
View(late_downregul[1:20, c(2,4,8)])
late_downregul$SYMBOL

early_upregul <- subset(top.table_neurons_IFN_05_early_limma_ordered, logFC >0)
early_downregul <- subset(top.table_neurons_IFN_05_early_limma_ordered, logFC < 0)

repvslate_upregul <- subset(top.table_neurons_IFN_05_repvslate_limma_ordered, logFC >0)
repvslate_downregul <- subset(top.table_neurons_IFN_05_repvslate_limma_ordered, logFC < 0)

genes <- late_upregul$ENTREZID
genes <- late_downregul$ENTREZID
genes <- early_upregul$ENTREZID
genes <- early_downregul$ENTREZID
genes <- repvslate_upregul$ENTREZID
genes <- repvslate_downregul$ENTREZID

#to assess GO enrichment with DAVID tool
writeLines(rownames(late_upregul), "limma_DE_list_IFN_neuronslate_up_20220702.txt")
writeLines(rownames(late_downregul), "limma_DE_list_IFN_neuronslate_down_20220702.txt")
writeLines(rownames(early_upregul), "limma_DE_list_IFN_neuronsearly_up_20220702.txt")
writeLines(rownames(early_downregul), "limma_DE_list_IFN_neuronsearly_down_20220702.txt")
writeLines(rownames(repvslate_upregul), "limma_DE_list_IFN_neuronsrepvslate_up_20220702.txt")
writeLines(rownames(repvslate_downregul), "limma_DE_list_IFN_neuronsrepvslate_down_20220702.txt")

# enrichment in term of BP of the set of significantly DE genes
OrgDb <- org.Hs.eg.db
# get the "genes" variable for the treatment you want to test
ego2 <- clusterProfiler::enrichGO(gene          = genes,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15)

## GSEA
geneList <- top.table_neurons_IFN_05_early_limma$logFC
names(geneList) <- rownames(top.table_neurons_IFN_05_early_limma)
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

geneList <- top.table_neurons_IFN_05_late_limma$logFC
names(geneList) <- rownames(top.table_neurons_IFN_05_late_limma)
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

geneList <- top.table_neurons_IFN_05_repvslate_limma$logFC
names(geneList) <- rownames(top.table_neurons_IFN_05_repvslate_limma)
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

#gene set enrichment; use the function with the list corresponding to the treatment you want to study
gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
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

#####
#select NPCs only
txi_NPC <- lapply(txi, function(x) if(is.matrix(x)) return(x[,1:6]) else return(x))

NPC_coldata <-subset(coldata1, subset = cell_type =="NPC")
NPC_coldata$D30 <-NULL
NPC_coldata$cell_type <-NULL

#convert the data as factors
NPC_coldata$cell_line <- factor(NPC_coldata$cell_line, levels = c("M1","M2", "M3"))
NPC_coldata$treatment <- factor(NPC_coldata$D18, levels = c("untreated","treated"))
colnames(txi_NPC$counts) <-rownames(NPC_coldata)

counts_from_kallisto_NPC <- txi_NPC$counts
colnames(counts_from_kallisto_NPC) <- rownames(NPC_coldata)

# filter genes by number of counts
isexpr <- rowSums(counts_from_kallisto_NPC) >= 10
#isexpr = rowSums(cpm(counts_from_kallisto_NPC)>0.1) >= 5

# Create DGEList object that stores counts
geneExpr_NPC = DGEList( counts_from_kallisto_NPC[isexpr,] )
geneExpr_NPC = calcNormFactors( geneExpr_NPC )

#add df genes to geneExpr
geneid <-rownames(geneExpr_NPC)
genes <-AnnotationDbi::select(org.Hs.eg.db, keys = geneid, columns = c("SYMBOL", "ENSEMBL"), keytype = "ENTREZID")
genes <-genes[!duplicated(genes$ENTREZID),]
geneExpr_NPC$genes <-genes #20526

design_NPC <- model.matrix(~treatment, NPC_coldata)
vobj_temp_NPC <- voom(geneExpr_NPC, design_NPC, plot=TRUE)
dupcor_NPC <-duplicateCorrelation(vobj_temp_NPC, design_NPC, block =NPC_coldata$cell_line)
vobj_NPC <- voom(geneExpr_NPC, design_NPC, plot=FALSE, block =NPC_coldata$cell_line, correlation = dupcor_NPC$consensus.correlation )
dupcor_NPC <-duplicateCorrelation(vobj_NPC, design_NPC, block =NPC_coldata$cell_line)
fitdupcor_NPC <-lmFit(vobj_NPC,design_NPC, block =NPC_coldata$cell_line, correlation = dupcor_NPC$consensus.correlation ) 
fitdupcor_NPC <-eBayes(fitdupcor_NPC)

head(fitdupcor_NPC$coefficients,5)

top.table_NPC_IFN_limma  <- topTable(
  fitdupcor_NPC,
  n=Inf,
  adjust.method = "BH",
  coef ='treatmenttreated'
)
top.table_NPC_IFN_limma_05 <- subset(top.table_NPC_IFN_limma, adj.P.Val <0.05 )
top.table_NPC_IFN_limma_05_ordered<- top.table_NPC_IFN_limma_05[order(abs(top.table_NPC_IFN_limma_05$logFC), decreasing=TRUE),]

write.csv(top.table_NPC_IFN_limma_05_ordered, file= "limma_DE_List_NPC_IFN20220702.csv", row.names = T)

gene_studied <-"SORCS3"
line_gene <-grep(gene_studied,top.table_NPC_IFN_limma_05_ordered$SYMBOL)
line_gene

NPC_IFN_nonDEgenes <- subset(top.table_NPC_IFN_limma, adj.P.Val > 0.05 )
writeLines(NPC_IFN_nonDEgenes$SYMBOL, "NPC_IFN_nonDEgenes_limma_20220702.txt")

#Go enrichment
NPC_upregul <- subset(top.table_NPC_IFN_limma_05_ordered, logFC >0)
NPC_downregul <- subset(top.table_NPC_IFN_limma_05_ordered, logFC < 0)

View(NPC_upregul[1:20,c(2,4,7)])
View(NPC_downregul[1:20,c(2,4,7)])

genes <- NPC_upregul$ENTREZID
genes <- NPC_downregul$ENTREZID

#to assess GO enrichment with DAVID tool
writeLines(rownames(NPC_upregul), "limma_DE_list_IFN_NPC_up_20220702.txt")
writeLines(rownames(NPC_downregul), "limma_DE_list_IFN_NPC_down_20220702.txt")


# enrichment in term of BP of the set of significantly DE genes
#define the "genes" variable for up- or -downregulated genes
ego2 <- clusterProfiler::enrichGO(gene          = genes,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15)

## GSEA
geneList_NPC <- top.table_NPC_IFN_limma_05$logFC
names(geneList_NPC) <- rownames(top.table_NPC_IFN_limma_05)
geneList_NPC <- na.omit(geneList_NPC)
geneList_NPC <- sort(geneList_NPC, decreasing = TRUE)

#gene set enrichment
gse <- gseGO(geneList=geneList_NPC, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPermSimple = 10000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=15, split=".sign", font.size = 6) + facet_grid(.~.sign)

######
#overlap avec list DESeq2
genes <-NULL
for (name in rownames(top.table_NPC_IFN_limma_05) ) {
  if (name %in% rownames(DE_list_NPC_IFN )){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #1351


go.obj <-newGeneOverlap(top.table_NPC_IFN_limma_05$SYMBOL,DE_list_NPC_IFN$SYMBOL, 17883 )
go.obj
go.obj <-testGeneOverlap(go.obj)
#GeneOverlap object:
#listA size=2218
#listB size=1701
#Intersection size=1351
#Overlapping p-value=0e+00
#Jaccard Index=0.5

overlap_NPC_DESeq2_limma <- getIntersection(go.obj)
writeLines(overlap_NPC_DESeq2_limma , "overlap_NPC_DESeq2_limma20220702.txt")

#overlap  late
genes <-NULL
DEseq_late <-rownames(results(dds_late))
for (name in rownames(top.table_neurons_IFN_late_limma)) {
  if (name %in% DEseq_late ){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #


go.obj <-newGeneOverlap(top.table_neurons_IFN_05_late_limma_ordered$SYMBOL,DE_list_late_IFN$SYMBOL, 20426 )
go.obj
go.obj <-testGeneOverlap(go.obj)
#GeneOverlap object:
#GeneOverlap object:
#  listA size=549
#listB size=724
#Intersection size=374
#Overlapping p-value=0e+00
#Jaccard Index=0.4

overlap_late_DESeq2_limma <- getIntersection(go.obj)
writeLines(overlap_late_DESeq2_limma , "overlap_late_DESeq2_limma20220702.txt")

#overlap  early
genes <-NULL
DEseq_early <-rownames(results(dds_early))
for (name in rownames(top.table_neurons_IFN_early_limma)) {
  if (name %in% DEseq_early){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #


go.obj <-newGeneOverlap(top.table_neurons_IFN_05_early_limma_ordered$SYMBOL,DE_list_early_IFN$SYMBOL, 20401 )
go.obj
go.obj <-testGeneOverlap(go.obj)
#GeneOverlap object:
#listA size=4
#listB size=95
#Intersection size=3
#Overlapping p-value=3.9e-07
#Jaccard Index=0.0

overlap_early_DESeq2_limma <- getIntersection(go.obj)
writeLines(overlap_early_DESeq2_limma, "overlap_early_DESeq2_limma20220702.txt")
overlap_early_DESeq2_limma

#overlap rep vs late
genes <-NULL
DEseq_rep <-rownames(results(dds_repvslate))
for (name in rownames(top.table_neurons_IFN_repvslate_limma)) {
  if (name %in% DEseq_rep){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #21153
#GeneOverlap object:
#listA size=22
#listB size=191
#Intersection size=20
#Overlapping p-value=3.1e-38
#Jaccard Index=0.1

go.obj <-newGeneOverlap(top.table_neurons_IFN_05_repvslate_limma_ordered$SYMBOL,DE_list_repvslate_IFN$SYMBOL, 21153 )
go.obj
go.obj <-testGeneOverlap(go.obj)
#GeneOverlap object:


overlap_repvslate_DESeq2_limma <- getIntersection(go.obj)
writeLines(overlap_repvslate_DESeq2_limma , "overlap_repvslate_DESeq2_limma20220702.txt")

overlap_repvslate_DESeq2_limma 

#overlap with DESeq2 genes downreguated by late treatment
overlap_down_late <- NULL
for (gene in DE_list_under_late$SYMBOL[1:100]){
  if (gene %in% late_downregul$SYMBOL)
    overlap_down_late <-c(overlap_down_late, gene)
}
overlap_down_late
length(overlap_down_late)

#overlap in downregulated geneset  by IFN-g in NPCs found by limma and deseq2
top100_down_limma <-NPC_downregul$SYMBOL[1:100]
overlap_down_NPC <- NULL
for (gene in DE_list_NPC_IFN_under$SYMBOL[1:100]){
  if (gene %in% top100_down_limma)
    overlap_down_NPC <-c(overlap_down_NPC, gene)
}
overlap_down_NPC
length(overlap_down_NPC)

