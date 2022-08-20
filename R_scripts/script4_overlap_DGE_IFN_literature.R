#install.packages("readxl")
library(readxl)

#install.packages("GeneOverlap")
library(GeneOverlap)

#install.packages("GenomicFeatures")
library("GenomicFeatures")

#install.packages("AnnotationDbi")
library("AnnotationDbi")

#install.packages("org.Hs.eg.db")
library("org.Hs.eg.db")

#BiocManager::install("clusterProfiler")
library("clusterProfiler")
library(enrichplot)

#BiocManager::install("ggplot2")
library("ggplot2")

#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

#install.packages("gplots")
library(gplots)

########################################
#Gandal's dataset- import spreadsheet with DGE
Gandal_DE <- read_excel("~/Desktop/MRes/20220224/NIHMS1007359-supplement-Table_S1.xlsx", sheet ="DGE", col_names =TRUE)
colnames(Gandal_DE)
#convert to df
df_DGE_Gandal <-as.data.frame(Gandal_DE, col.names = colnames(Gandal_DE))
View(df_DGE_Gandal)

#only keep relevant columns
df_DGE_Gandal_ASD <- df_DGE_Gandal[,c(1, 6, 8,14,18,19)]
length(df_DGE_Gandal_ASD$ensembl_gene_id)

#order the data by log2FC
df_DGE_Gandal_ASD2 <- df_DGE_Gandal_ASD[order(abs(df_DGE_Gandal_ASD$ASD.log2FC),na.last=TRUE, decreasing=TRUE),]
length(df_DGE_Gandal_ASD2$ensembl_gene_id)

#keep only isoforms when FDR-adjusted pvalue <0.05
df_DGE_Gandal_ASD3 <- df_DGE_Gandal_ASD2[which(df_DGE_Gandal_ASD2$ASD.fdr < 0.05),]

genes_DE_gandal <-df_DGE_Gandal_ASD3$gene_name
genes_DE_gandal
length(unique(genes_DE_gandal))

#look for a particular gene in the list
gene_studied <-"SYT13"
line_gene <-grep(gene_studied,df_DGE_Gandal_ASD3$gene_name)
line_gene

df_DGE_Gandal_ASD$ENTREZID = mapIds(org.Hs.eg.db,
                                     keys=df_DGE_Gandal_ASD$ensembl_gene_id, 
                                     column="ENTREZID",
                                     keytype="ENSEMBL",
                                     multiVals="first")

df_DGE_Gandal_ASD <-df_DGE_Gandal_ASD[!duplicated(df_DGE_Gandal_ASD$ENTREZID),]

df_DGE_Gandal_ASD3$ENTREZID = mapIds(org.Hs.eg.db,
                                     keys=df_DGE_Gandal_ASD3$ensembl_gene_id, 
                                     column="ENTREZID",
                                     keytype="ENSEMBL",
                                     multiVals="first")

df_DGE_Gandal_ASD3 <-df_DGE_Gandal_ASD3[!duplicated(df_DGE_Gandal_ASD3$ENTREZID),]
df_DGE_Gandal_ASD3 <-na.omit(df_DGE_Gandal_ASD3)

df_DGE_Gandal_ASD4 <- format(df_DGE_Gandal_ASD3, digits=2)
View(df_DGE_Gandal_ASD4[1:20,c(3,4,6)])

Gandal_overexpressed <- subset(df_DGE_Gandal_ASD3, ASD.log2FC >= 0)
Gandal_underexpressed <- subset(df_DGE_Gandal_ASD3, ASD.log2FC < 0)

#clusterprofiler: overrepresentation analysis
OrgDb <- org.Hs.eg.db
genes <- df_DGE_Gandal_ASD3$ENTREZID
Gandal_genes_over <- Gandal_overexpressed$ENTREZID
Gandal_genes_under <- Gandal_underexpressed$ENTREZID
#for david assessment
writeLines(Gandal_genes_over, "Gandal_genes_over_E_20220709.txt")
writeLines(Gandal_genes_under, "Gandal_genes_under_E_20220709.txt")

# enrichment in term of BP of the set of significantly DE genes
ego2 <- clusterProfiler::enrichGO(gene          = Gandal_genes_under,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

ego2 <- clusterProfiler::enrichGO(gene          = Gandal_genes_over,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

## GSEA
#Gene list with Fold Change data
geneList <- df_DGE_Gandal_ASD3$ASD.log2FC
names(geneList) <- as.character(unique(df_DGE_Gandal_ASD3$ENTREZID))
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

#gene set enrichment; determines whether a pre-defined set of genes 
#shows statistically significant, concordant differences between two biological states
gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPermSimple = 1000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=10, split=".sign", font.size = 6) + facet_grid(.~.sign)

#enrichment in Sfari genes
Sfari_genes <- read.csv("~/Desktop/MRes/20220611/SFARI.csv", header =TRUE, stringsAsFactors = FALSE)
Sfari_high_score <- subset(Sfari_genes, gene.score %in% c(1,2))
list_SFARI <- Sfari_high_score$gene.symbol

go.obj <-newGeneOverlap(list_SFARI,genes_DE_gandal , 25774 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
# listA size=895
#listB size=1611
#Intersection size=96
#Overlapping p-value=1.6e-07
#Jaccard Index=0.0

overlap_Gandal_sfari <- getIntersection(go.obj)
writeLines(overlap_Gandal_sfari, "overlap_Gandal_sfari20220709.txt")
overlap_Gandal_sfari
#"ABCA10"   "ANK3"     "AHNAK"    "ANKS1B"   "APBB1"    "ARHGAP5"  "ARHGEF9"  "ATP1A1"   "ATP2B2"   "BCL11A"   "BRAF"     "C4B"      "CACNA2D3"
#[14] "CADM1"    "CADPS"    "CADPS2"   "CAMK4"    "CCNG1"    "CHD3"     "CD99L2"   "CHRNB3"   "CLCN4"    "CTTNBP2"  "CYFIP1"   "DAGLA"    "DDHD2"   
#[27] "DLGAP1"   "DPP10"    "DPYSL3"   "EIF4E"    "FABP5"    "GABRB3"   "FGFR1"    "GALNT10"  "FRG1"     "GABRB2"   "GPD2"     "GPR37"    "GPR85"   
#[40] "GRID2"    "GRIN2A"   "HIVEP2"   "HLA-B"    "HMGN1"    "HS3ST5"   "HSD11B1"  "ICA1"     "ITGB3"    "KCNB1"    "KCNK7"    "KCNS3"    "KDM1B"   
#[53] "LAMB1"    "LZTR1"    "MAOB"     "MAPK3"    "MYO16"    "NBEA"     "NEGR1"    "NIPA2"    "NRXN3"    "NUAK1"    "NUP133"   "OPHN1"    "P4HA2"   
#[66] "PCDH15"   "PLAUR"    "PLCB1"    "PRICKLE1" "PRKAR1B"  "PRKCB"    "PSMD12"   "PTCHD1"   "RBFOX1"   "RGS7"     "SAE1"     "SCN1A"    "SCN2A"   
#[79] "SCN8A"    "SCN9A"    "SDC2"     "SERPINE1" "SEZ6L2"   "SLC24A2"  "SLC25A12" "SLC25A39" "SLC27A4"  "SMG6"     "SNAP25"   "SLC22A15" "SYT17"   
#[92] "TBL1X"    "UBE3A"    "XPO1"     "ZNF385B"  "ZNF711"  

#overlap between Gandal's dataset and IFN
#NPCs-number of genes
genes <-rownames(tot_genes_NPC)
for (name in df_DGE_Gandal_ASD$ENTREZID ) {
  if (!name %in% rownames(tot_genes_NPC)){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #23377

go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL, genes_DE_gandal, 23377 )
go.obj <-testGeneOverlap(go.obj)
go.obj

overlap_NPC_gandal <- getIntersection(go.obj)
overlap_NPC_gandal
#GeneOverlap object:
#  listA size=1701
#listB size=1611
#Intersection size=222
#Overlapping p-value=2.5e-21
#Jaccard Index=0.1

#overlap neurons gandal
DE_genes_neurons_IFN <- c(DE_list_late_IFN$SYMBOL, DE_list_early_IFN$SYMBOL, DE_list_repvslate_IFN$SYMBOL) #from IFN DE analysis
DE_genes_neurons_IFN <- unique(DE_genes_neurons_IFN)
DE_genes_neurons_IFN 

#total_genes 
total_genes_neurons <-rownames(tot_genes_late) #from DE analysis IFN
for (name in rownames(tot_genes_early )) {
  if (!name %in% total_genes_neurons){
    total_genes_neurons <- c(total_genes_neurons, name)
  }
}
for (name in rownames(tot_genes_rep )) {
  if (!name %in% total_genes_neurons){
    total_genes_neurons <- c(total_genes_neurons, name)
  }
}

total_genes_neurons<-unique(total_genes_neurons)
length(total_genes_neurons) #21581

total_genes_gandal <-total_genes_neurons
for (name in df_DGE_Gandal$ENTREZID ) {
  if (!name %in% total_genes_neurons){
    total_genes_gandal <- c(total_genes_gandal, name)
  }
}
total_genes_gandal<-unique(total_genes_gandal)
length(total_genes_gandal) #21642

go.obj <-newGeneOverlap(DE_genes_neurons_IFN,genes_DE_gandal , 21642 )
go.obj <-testGeneOverlap(go.obj)
go.obj
#GeneOverlap object:
#  listA size=840
#listB size=1611
#Intersection size=130
#Overlapping p-value=6.3e-16
#Jaccard Index=0.1
overlap_neuron_gandal <- getIntersection(go.obj)

# GO analysis with clusterprofiler
OrgDb <- org.Hs.eg.db
#overlap DE Gandal- DE following IFN-g in NPCs
overlap_NPC_Gandal_df <- data.frame(
  symbol <- overlap_NPC_gandal
)
colnames(overlap_NPC_Gandal_df) <- "SYMBOL"

overlap_NPC_Gandal_df$ENTREZID = mapIds(org.Hs.eg.db,
                                    keys=overlap_NPC_Gandal_df$SYMBOL, 
                                    column="ENTREZID",
                                    keytype="SYMBOL",
                                    multiVals="first")

# enrichment in term of BP of the set of overlapping DE genes by IFN in NPC and in Gandal's dataset
ego2_NPC <- clusterProfiler::enrichGO(gene          = overlap_NPC_Gandal_df$ENTREZID,
                                      OrgDb         = OrgDb,
                                      ont           = "BP",
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.01, 
                                      readable      = TRUE)

x2 <- pairwise_termsim(ego2_NPC)
emapplot(x2, showCategory = 15, node_label = "category")

#overlap DE Gandal- DE following IFN-g in neurons
overlap_neuron_Gandal_df <- data.frame(
  symbol <- overlap_neuron_gandal
)
colnames(overlap_neuron_Gandal_df) <- "SYMBOL"

overlap_neuron_Gandal_df$ENTREZID = mapIds(org.Hs.eg.db,
                                        keys=overlap_neuron_Gandal_df$SYMBOL, 
                                        column="ENTREZID",
                                        keytype="SYMBOL",
                                        multiVals="first")

# enrichment in term of BP of the set of the overlap
ego2_neurons <- clusterProfiler::enrichGO(gene          = overlap_neuron_Gandal_df$ENTREZID,
                                      OrgDb         = OrgDb,
                                      ont           = "BP",
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.01, 
                                      readable      = TRUE)

x2 <- pairwise_termsim(ego2_neurons)
emapplot(x2, showCategory = 15, node_label = "category")


#import spreadsheet with DGE from Haney et al.
DGE_Haney  <- read_excel("~/Desktop/MRes/20220313/Haney_table3.xlsx", sheet ="DEGene_Statistics", col_names =TRUE)
#convert to df
df_DGE_Haney <-as.data.frame(DGE_Haney , col.names = colnames(DGE_Haney ), row.names = rownames(DGE_Haney ))
colnames(df_DGE_Haney)
df_DGE_Haney$ENTREZID = mapIds(org.Hs.eg.db,
                                     keys=df_DGE_Haney$ensembl_gene_id, 
                                     column="ENTREZID",
                                     keytype="ENSEMBL",
                                     multiVals="first")

#only keep relevant columns
df_DGE_Haney_2 <- df_DGE_Haney[,c(1,2,10,11, 12)]

#keep only isoforms when FDR <0.05
df_DGE_Haney_3 <- df_DGE_Haney_2[which(df_DGE_Haney_2$WholeCortex_ASD_FDR < 0.05),]

#order the data by log2FC
df_DGE_Haney_4 <- df_DGE_Haney_3[order(abs(df_DGE_Haney_3$WholeCortex_ASD_logFC),na.last=TRUE, decreasing=TRUE),]
#be more strict on Log2FC cutoff
df_DGE_Haney_5 <- subset(df_DGE_Haney_4,abs(WholeCortex_ASD_logFC) > 0.3)

length(df_DGE_Haney_4$external_gene_name)

df_DGE_Haney_6 <- format(df_DGE_Haney_4, digits=2)
View(df_DGE_Haney_7[1:20,c(2,3,4)])

#look for a particular gene in the list
gene_studied <-"SYT13"
line_gene <-grep(gene_studied,df_DGE_Haney_4$external_gene_name)
line_gene

genes_DE_haney <-df_DGE_Haney_4$external_gene_name
genes_DE_haney2 <-df_DGE_Haney_5$external_gene_name
genes_DE_haney
length(genes_DE_haney2) #2077

DE_haney_over <- subset(df_DGE_Haney_4, WholeCortex_ASD_logFC > 0 )
DE_haney_over2 <- subset(df_DGE_Haney_5, WholeCortex_ASD_logFC > 0 )
length(DE_haney_over$external_gene_name ) #1942

DE_haney_under <- subset(df_DGE_Haney_4, WholeCortex_ASD_logFC < 0 )
DE_haney_under2 <- subset(df_DGE_Haney_5, WholeCortex_ASD_logFC < 0 )
length(DE_haney_under$external_gene_name ) #2277

#for david assessment
writeLines(DE_haney_over$ENTREZID, "Haney_genes_over_E_20220709.txt")
writeLines(DE_haney_under$ENTREZID, "Haney_genes_under_E_20220709.txt")

# enrichment in term of BP of the set of significantly DE genes
#clusterprofiler: overrepresentation analysis
OrgDb <- org.Hs.eg.db
genes <- DE_haney_over$ENTREZID
genes <- DE_haney_under$ENTREZID

#run with over- and underexpressed genes
ego2 <- clusterProfiler::enrichGO(gene          = genes,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)
x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

## GSEA
geneList <- df_DGE_Haney_4$WholeCortex_ASD_logFC
names(geneList) <- as.character(unique(df_DGE_Haney_4$ENTREZID))
geneList2 <- df_DGE_Haney_5$WholeCortex_ASD_logFC
names(geneList2) <- as.character(unique(df_DGE_Haney_5$ENTREZID))
geneList <- na.omit(geneList)
geneList2 <- na.omit(geneList2)
geneList <- sort(geneList, decreasing = TRUE)
geneList2 <- sort(geneList2, decreasing = TRUE)


#gene set enrichment
#also assess with geneList2
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

dotplot(gse, showCategory=7, split=".sign", font.size = 6) + facet_grid(.~.sign)

#overlap Haney and Sfari
#go.obj <-newGeneOverlap(list_SFARI,genes_DE_haney, 24836 )
go.obj <-newGeneOverlap(list_SFARI,genes_DE_haney2, 24836 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
# listA size=895
#listB size=2075
#Intersection size=105
#Overlapping p-value=2.4e-04
#Jaccard Index=0.0

overlap_Haney_sfari <- getIntersection(go.obj)
writeLines(overlap_Haney_sfari, "overlap_Haney001_sfari20220709.txt")
overlap_Haney_sfari

#overlap between Haney's dataset and IFN
#NPCs-number of genes

genes <-rownames(tot_genes_NPC)
for (name in df_DGE_Haney$ENTREZID ) {
  if (!name %in% rownames(tot_genes_NPC)){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #23311

go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL, genes_DE_haney, 23311 )
go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL, genes_DE_haney2, 23311 )
go.obj <-testGeneOverlap(go.obj)
go.obj

overlap_NPC_haney <- getIntersection(go.obj)

#GeneOverlap object:
#  listA size=1701
#listB size=2075
#Intersection size=266
#Overlapping p-value=7.4e-21
#Jaccard Index=0.1

#overlap_NPC_haney2 <- getIntersection(go.obj)
#overlap_NPC_haney2
#GeneOverlap object:
#listA size=1701     listB size=4216     Intersection size=476     Overlapping p-value=1.1e-25       Jaccard Index=0.1

#overlap neurons haney
DE_genes_neurons_IFN <- c(DE_list_late_IFN$SYMBOL, DE_list_early_IFN$SYMBOL, DE_list_repvslate_IFN$SYMBOL)
DE_genes_neurons_IFN <- unique(DE_genes_neurons_IFN)
DE_genes_neurons_IFN 

#total_genes 
tot_genes_late <- as.data.frame(results(dds_late))
tot_genes_early <- as.data.frame(results(dds_early))
tot_genes_rep<- as.data.frame(results(dds_repvslate))

total_genes_neurons <-rownames(tot_genes_late) 
for (name in rownames(tot_genes_early )) {
  if (!name %in% total_genes_neurons){
    total_genes_neurons <- c(total_genes_neurons, name)
  }
}
for (name in rownames(tot_genes_rep )) {
  if (!name %in% total_genes_neurons){
    total_genes_neurons <- c(total_genes_neurons, name)
  }
}

total_genes_neurons<-unique(total_genes_neurons)
length(total_genes_neurons) #21581

total_genes_haney <-total_genes_neurons
for (name in df_DGE_Haney$ENTREZID ) {
  if (!name %in% total_genes_neurons){
    total_genes_haney <- c(total_genes_haney, name)
  }
}
total_genes_haney<-unique(total_genes_haney)
length(total_genes_haney) #24181

go.obj <-newGeneOverlap(DE_genes_neurons_IFN,genes_DE_haney2 , 24181 )
go.obj <-testGeneOverlap(go.obj)
go.obj
#GeneOverlap object:
#  listA size=840
#listB size=2075
#Intersection size=152
#Overlapping p-value=4.1e-19
#Jaccard Index=0.1
overlap_neuron_haney <- getIntersection(go.obj)

#go.obj <-newGeneOverlap(DE_genes_neurons_IFN,genes_DE_haney2 , 24181 )
#go.obj <-testGeneOverlap(go.obj)
#go.obj
#GeneOverlap object:         listA size=840     listB size=4216    Intersection size=259     Overlapping p-value=3e-22     Jaccard Index=0.1

# GO analysis of the overlaps with clusterprofiler
OrgDb <- org.Hs.eg.db

#overlap genes DE in Haney's- genes DE by IFN in NPCs
overlap_NPC_haney_df <- data.frame(
  symbol <- overlap_NPC_haney
)
colnames(overlap_NPC_haney_df) <- "SYMBOL"

overlap_NPC_haney_df$ENTREZID = mapIds(org.Hs.eg.db,
                                        keys=overlap_NPC_haney_df$SYMBOL, 
                                        column="ENTREZID",
                                        keytype="SYMBOL",
                                        multiVals="first")

# enrichment in term of BP of the overlap haney-NPC (genes DE by IFN-g)
ego2_NPC <- clusterProfiler::enrichGO(gene          = overlap_NPC_haney_df$ENTREZID,
                                      OrgDb         = OrgDb,
                                      ont           = "BP",
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.01, 
                                      readable      = TRUE)

x2 <- pairwise_termsim(ego2_NPC)
emapplot(x2, showCategory = 15, node_label = "category")

#overlap genes DE in Haney's- genes DE by IFN in neurons
overlap_neuron_Haney_df <- data.frame(
  symbol <- overlap_neuron_haney
)
colnames(overlap_neuron_Haney_df) <- "SYMBOL"

overlap_neuron_Haney_df$ENTREZID = mapIds(org.Hs.eg.db,
                                           keys=overlap_neuron_Haney_df$SYMBOL, 
                                           column="ENTREZID",
                                           keytype="SYMBOL",
                                           multiVals="first")

# enrichment in term of BP of the overlap
ego2_neurons <- clusterProfiler::enrichGO(gene          = overlap_neuron_Haney_df$ENTREZID,
                                          OrgDb         = OrgDb,
                                          ont           = "BP",
                                          pAdjustMethod = "fdr",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.01, 
                                          readable      = TRUE)

x2 <- pairwise_termsim(ego2_neurons)
emapplot(x2, showCategory = 15, node_label = "category")

#import spreadsheet with DGE from Gupta et al.
table_1 <- read_excel("~/Desktop/MRes/20220629/from_literature/Gupta_DE.xlsx", col_names =TRUE)
colnames(table_1) <-table_1[1,]
table_1 <- table_1[2:750,]
colnames(table_1) <- c("gene_name", "p_adj", "Log2FC", "beta", "std_error")

table_1$Log2FC <- as.numeric(table_1$Log2FC)
table_1$p_adj <- as.numeric(table_1$p_adj)
table_1 <- table_1[order(abs(table_1$Log2FC), na.last=TRUE,decreasing=TRUE),]

table_2 <- as.data.frame(table_1)
table_3 <- format(table_2, digits=2)
View(table_3[1:20, c(1,2,3)])

list_DE_Gupta <- table_1$gene_name

#look for a particular gene in the list
gene_studied <-"FERMT2"
line_gene <-grep(gene_studied, table_1$gene_name)
line_gene

table_1$ENTREZID = mapIds(org.Hs.eg.db,
                          keys=table_1$gene_name, 
                          column="ENTREZID",
                          keytype="SYMBOL",
                          multiVals="first")

table_1 <-table_1[!duplicated(table_1$ENTREZID),]
table_1 <-na.omit(table_1)
genes <- table_1$ENTREZID
writeLines(table_1$ENTREZID, "Gupta_genes_over_E_20220709.txt")

# enrichment in term of BP of the set of significantly DE genes
ego2 <- clusterProfiler::enrichGO(gene          = genes,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)

x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

## GSEA
geneList <- table_1$Log2FC
names(geneList) <- as.character(unique(table_1$ENTREZID))
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPermSimple = 1000,
             minGSSize = 3, 
             maxGSSize = 800, 
             eps = 0,
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory=10, split=".sign", font.size = 6) + facet_grid(.~.sign)

#overlap with sfari high score
go.obj <-newGeneOverlap(list_SFARI, list_DE_Gupta, 24000 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
# listA size=895
#listB size=749
#Intersection size=47
#Overlapping p-value=3.8e-04
#Jaccard Index=0.0

overlap_Gupta_sfari <- getIntersection(go.obj)
writeLines(overlap_Gupta_sfari, "overlap_Gupta_sfari20220709.txt")
overlap_Gupta_sfari
#[1] "AHDC1"    "ARX"      "AP2S1"    "CELF4"    "CNTNAP4"  "CUX1"     "DAGLA"    "CPEB4"    "DPP6"     "DPYSL2"   "ELAVL3"   "FOXP2"    "GRIA2"   
#[14] "HOMER1"   "HS3ST5"   "IL1RAPL1" "KCND3"    "KIRREL3"  "LRRC1"    "MCM6"     "MDGA2"    "MAP1A"    "NCOR1"    "NFIA"     "NFIB"     "NINL"    
#[27] "PDE1C"    "POMT1"    "PRKCB"    "PSD3"     "PXDN"     "SCP2"     "SETD2"    "SEZ6L2"   "SHANK3"   "SLC12A5"  "SMG6"     "SYAP1"    "SYN2"    
#[40] "SYT17"    "TBCK"     "TCF4"     "UIMC1"    "YTHDC2"   "ZBTB20"   "ZMYND8"   "ZNF827"  

#overlap between Gupta's dataset and IFN
#NPCs-number of genes

genes <-rownames(tot_genes_NPC)
for (name in table_1$Log2FC ) {
  if (!name %in% rownames(tot_genes_NPC)){
    genes <- c(genes, name)
  }
}
genes<-unique(genes)
length(genes) #21210

go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL,list_DE_Gupta, 21210 )
go.obj <-testGeneOverlap(go.obj)
go.obj

overlap_NPC_gupta <- getIntersection(go.obj)
overlap_NPC_gupta
#GeneOverlap object:
#  listA size=1701
#listB size=749
#Intersection size=70
#Overlapping p-value=0.1
#Jaccard Index=0.1

#overlap neurons 
DE_genes_neurons_IFN <- c(DE_list_late_IFN$SYMBOL, DE_list_early_IFN$SYMBOL, DE_list_repvslate_IFN$SYMBOL)
DE_genes_neurons_IFN <- unique(DE_genes_neurons_IFN)
DE_genes_neurons_IFN 

#total_genes
total_genes_neurons <-rownames(tot_genes_late) 
for (name in rownames(tot_genes_early )) {
  if (!name %in% total_genes_neurons){
    total_genes_neurons <- c(total_genes_neurons, name)
  }
}
for (name in rownames(tot_genes_rep )) {
  if (!name %in% total_genes_neurons){
    total_genes_neurons <- c(total_genes_neurons, name)
  }
}

total_genes_neurons<-unique(total_genes_neurons)
length(total_genes_neurons) #21581

total_genes_gupta <-total_genes_neurons
for (name in table_1$ENTREZID ) {
  if (!name %in% total_genes_neurons){
    total_genes_gupta <- c(total_genes_gupta, name)
  }
}
total_genes_gupta<-unique(total_genes_gupta)
length(total_genes_gupta) #21669

go.obj <-newGeneOverlap(DE_genes_neurons_IFN,list_DE_Gupta, 21669 )
go.obj <-testGeneOverlap(go.obj)
go.obj
#GeneOverlap object:
#  listA size=840
#listB size=749
#Intersection size=43
#Overlapping p-value=6.9e-03
#Jaccard Index=0.0
overlap_neuron_gupta <- getIntersection(go.obj)

overlap_neuron_Gupta_df <- data.frame(
  symbol <- overlap_neuron_gupta
)
colnames(overlap_neuron_Gupta_df) <- "SYMBOL"
overlap_neuron_Gupta_df$ENTREZID = mapIds(org.Hs.eg.db,
                            keys=overlap_neuron_Gupta_df$SYMBOL, 
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")

# enrichment in term of BP of the set of the overlap
ego2_neurons <- clusterProfiler::enrichGO(gene          = overlap_neuron_Gupta_df$ENTREZID,
                                          OrgDb         = OrgDb,
                                          ont           = "BP",
                                          pAdjustMethod = "fdr",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.01, 
                                          readable      = TRUE)
x2 <- pairwise_termsim(ego2_neurons)
emapplot(x2, showCategory = 10, node_label = "category")


#overlap NPCs
overlap <- NULL
overlap1 <- NULL
overlap2 <- NULL
overlap3 <- NULL

for (name in DE_list_NPC_IFN$SYMBOL) {
  if (name %in% top.table_NPC_IFN_limma_05_ordered$SYMBOL){
    overlap1 <- c(overlap1, name)
  }
}

for (name in overlap1) {
  if (name %in% genes_DE_gandal){
    overlap2 <- c(overlap2, name)
  }
}

for (name in overlap2) {
  if (name %in% list_DE_Gupta){
    overlap3 <- c(overlap3, name)
  }
}

for (name in overlap3) {
  if (name %in% genes_DE_haney){
    overlap <- c(overlap, name)
  }
}

overlap
#""BST2"   "IFI6"   "C1R"    "UBE2L6" "MOB3C"  "CLU"    "FERMT2" "LRP10"  "ASAP3"  

DGE_IFNg_ASD <- list(
  Gandal<- genes_DE_gandal,
  Gupta <- list_DE_Gupta,
  Gandal<-genes_DE_haney,
  IFN_DESEq2 <- DE_list_NPC_IFN$SYMBOL,
  IFN_limmavoom <- top.table_NPC_IFN_limma_05_ordered$SYMBOL
)

p <- ggVennDiagram(DGE_IFNg_ASD, category.names = c("Gandal","Gupta","Haney", "IFN-NPC-DESEq2", "IFN_NPC-limma"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)

v.table <- venn(DGE_IFNg_ASD)
print(v.table)

#overlap neurons
DE_limma_neurons <- unique(c(top.table_neurons_IFN_05_early_limma_ordered$SYMBOL, top.table_neurons_IFN_05_late_limma_ordered$SYMBOL,
                             top.table_neurons_IFN_05_repvslate_limma_ordered$SYMBOL))
overlap_N <- NULL
overlap1_N <- NULL
overlap2_N <- NULL
overlap3_N <-NULL

for (name in DE_genes_neurons_IFN ) {
  if (name %in% DE_limma_neurons){
    overlap1_N <- c(overlap1_N, name)
  }
}

for (name in overlap1_N) {
  if (name %in% genes_DE_gandal){
    overlap2_N <- c(overlap2_N, name)
  }
}

for (name in overlap2_N) {
  if (name %in% list_DE_Gupta){
    overlap3_N <- c(overlap3_N, name)
  }
}

for (name in overlap3_N) {
  if (name %in% genes_DE_haney){
    overlap_N <- c(overlap_N, name)
  }
}

overlap_N
#"SLC6A6"  "IFI6"    "BST2"    "UBE2L6"  "SYNPO2"  "FERMT2"  "FBXW7"   "SYT13"   "PRKCB"   "S100A10"
#"IFI6"   "UBE2L6" "SYNPO2" "FERMT2"

DGE_IFNg_ASD_neurons <- list(
  Gandal<- genes_DE_gandal,
  Gupta <- list_DE_Gupta,
  Gandal<-genes_DE_haney,
  IFN_DESeq<- DE_genes_neurons_IFN,
  IFN_limma <- DE_limma_neurons
)

p <- ggVennDiagram(DGE_IFNg_ASD_neurons, category.names = c("Gandal","Gupta","Haney", "IFN-neurons_DESEq2", "IFN_neurons_limma"), label = "count")

p + scale_fill_distiller(palette = "Blues", direction = 1)

v.table <- venn(DGE_IFNg_ASD)
print(v.table)


###### same with IFN list from limma
#overlap NPCs
overlap <- NULL
overlap1 <- NULL
overlap2 <- NULL
for (name in top.table_NPC_IFN_limma_05_ordered$SYMBOL) {
  if (name %in% genes_DE_gandal){
    overlap1 <- c(overlap1, name)
  }
}

for (name in overlap1) {
  if (name %in% list_DE_Gupta){
    overlap2 <- c(overlap2, name)
  }
}

for (name in overlap2) {
  if (name %in% genes_DE_haney){
    overlap <- c(overlap, name)
  }
}

overlap
#" "BST2"   "IFI6"   "C1R"    "UBE2L6" "MOB3C"  "RAB20"  "SYNPO2" "CLU"    "FERMT2" "MDK"    "ASAP3"  "LRP10"  

DGE_IFNg_ASD2 <- list(
  Gandal<- genes_DE_gandal,
  Gupta <- list_DE_Gupta,
  Gandal<-genes_DE_haney,
  IFN<- top.table_NPC_IFN_limma_05_ordered$SYMBOL
)

p <- ggVennDiagram(DGE_IFNg_ASD2, category.names = c("Gandal","Gupta","Haney", "IFN-g NPCs"), label = "count")
p + scale_fill_distiller(palette = "Blues", direction = 1)


#overlap neurons
overlap_N <- NULL
overlap1_N <- NULL
overlap2_N <- NULL
for (name in DE_limma_neurons ) {
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
  if (name %in% genes_DE_haney){
    overlap_N <- c(overlap_N, name)
  }
}

overlap_N
#"IFI6"   "MOB3C"  "UBE2L6" "SYNPO2" "FERMT2"

DGE_IFNg_ASD <- list(
  Gandal<- genes_DE_gandal,
  Gupta <- list_DE_Gupta,
  Gandal<-genes_DE_haney,
  IFN<- DE_genes_neurons_IFN
)


p <- ggVennDiagram(DGE_IFNg_ASD, category.names = c("Gandal","Gupta","Haney", "IFN-neurons"), label = "count")

p + scale_fill_distiller(palette = "Blues", direction = 1)


print(v.table)
v.table <- venn(DGE_IFNg_ASD)


#overlap between genes downregulated by interferon gamma and genes dysregulated in ASD post-mortem brain
#lists from DE analysis IFN DEseq2
IFN_downregul_genes <- c(rownames(df_DE_list_NPC_IFN_under),rownames(df_DE_list_under_late ),rownames(df_DE_list_under_early), rownames(df_DE_list_under_repvslate))
IFN_downregul_genes <- unique(IFN_downregul_genes)

literature_dysregul_genes <- c(genes_DE_gandal, list_DE_Gupta,genes_DE_haney)
literature_dysregul_genes <- unique(literature_dysregul_genes)

overlap_IFN_down_literature_dys <- NULL
for (name in IFN_downregul_genes) {
  if (name %in% literature_dysregul_genes){
    overlap_IFN_down_literature_dys <- c(overlap_IFN_down_literature_dys, name)
  }
}

overlap_IFN_down_literature_dys

df_overlap_IFN_down_literature_dys <- data.frame(
  symbol <-overlap_IFN_down_literature_dys
)
colnames(df_overlap_IFN_down_literature_dys) <- "SYMBOL"

df_overlap_IFN_down_literature_dys$ENTREZID = mapIds(org.Hs.eg.db,
                          keys=df_overlap_IFN_down_literature_dys$SYMBOL, 
                          column="ENTREZID",
                          keytype="SYMBOL",
                          multiVals="first")
#GO enrichment
OrgDb <- org.Hs.eg.db
genes <-df_overlap_IFN_down_literature_dys$ENTREZID

ego2 <- clusterProfiler::enrichGO(gene          = genes,
                                  OrgDb         = OrgDb,
                                  ont           = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.01, 
                                  readable      = TRUE)

x2 <- pairwise_termsim(ego2)
emapplot(x2, showCategory = 15, font.size = 6)

#ASD authors (obtained with normalised values given by RUVSeq)
#import spreadsheet with DGE from Griesi et al.
table_ASD_authors <- read_excel("~/Desktop/MRes/20220629/from_literature/griesi1.xlsx", col_names =TRUE)
colnames(table_ASD_authors)
df_DGE_Griesi <-as.data.frame(table_ASD_authors, col.names = colnames(table_ASD_authors), row.names = rownames(table_ASD_authors))
colnames(df_DGE_Griesi) <-df_DGE_Griesi[2,]
df_DGE_Griesi <-df_DGE_Griesi[3:22,]
genes_DE_griesi <- df_DGE_Griesi$Gene

#overlap with down IFN
go.obj <-newGeneOverlap(IFN_downregul_genes,genes_DE_griesi , 21642 )
go.obj <-testGeneOverlap(go.obj)
go.obj #Overlapping p-value=4.6e-03
overlap<- getIntersection(go.obj)
overlap #"ACSL6"   "KCNB1"   "SEMA3C"  "SLC12A5"

go.obj <-newGeneOverlap(DE_genes_neurons_IFN,genes_DE_griesi , 21642 )
go.obj <-testGeneOverlap(go.obj)
go.obj #Overlapping p-value=0.18
overlap<- getIntersection(go.obj)
overlap #"SEMA3C"  "SLC12A5 sfari"

go.obj <-newGeneOverlap(DE_list_NPC_IFN$SYMBOL,genes_DE_griesi , 21642 )
go.obj <-testGeneOverlap(go.obj)
go.obj #Overlapping p-value=0.47
overlap<- getIntersection(go.obj)
overlap #"ACSL6" "KCNB1 sfari"
