#install.packages("readxl")
library(readxl)

#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(ggplot2)

#install.packages("gplots")
library(gplots)

#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

########################
# use data from the literature ASD DTE-DTU
#read the table from Gandal's paper, DTU data
table_1 <- read_excel("~/Desktop/MRes/20220224/NIHMS1007359-supplement-Table_S1.xlsx", sheet ="DTU", col_names =TRUE)
df_DTU_Gandal <-as.data.frame(table_1, col.names = colnames(table_1))
#only keep relevant columns
df_DTU_Gandal_ASD <- df_DTU_Gandal[,c(1,2,3,4,6,13, 28)]
#order the data by beta-value (effect size) in ASD
df_DTU_Gandal_ASD <- df_DTU_Gandal_ASD[order(abs(df_DTU_Gandal_ASD$DTU.ASD.Value),na.last=TRUE, decreasing=TRUE),]

#keep only isoforms when pvalue <0.05
df_DTU_Gandal_ASD <- df_DTU_Gandal_ASD[which(df_DTU_Gandal_ASD$DTU.ASD.FDR < 0.05),]

#look for a particular gene in the list
gene_studied <-"ACTN1"
line_gene <-grep(gene_studied, df_DTU_Gandal_ASD$external_gene_id)
line_gene

list_Gandal_DTU <- df_DTU_Gandal_ASD$external_gene_id #609

#finding overlap between Gandal DTU et DTU IFN NPCs
overlap <- NULL
for (name in DTU_non_DE_IFN_NPC) {
  if (name %in% list_Gandal_DTU){
    overlap <- c(overlap, name)
  }
}
overlap_NPC_Gandal <-overlap   #"IL1RAP" "ISG20"  "ALG13"  "VMP1"   "SPHK2" 

#finding overlap between Gandal DTU et DTU IFN neurons
overlap_late <- NULL
for (name in list_Gandal_DTU) {
  if (name %in% DTU_non_DE_IFN_Late){
    overlap_late <- c(overlap_late, name)
  }
}
overlap_late_Gandal <- overlap_late
#P4HA1"  

overlap_early <- NULL
for (name in list_Gandal_DTU) {
  if (name %in% DTU_non_DE_IFN_early){
    overlap_early <- c(overlap_early, name)
  }
}
overlap_early_Gandal <- overlap_early
#"SLC30A6" "RPL22L1" "FAM13A"  "PLS3"    "ADAM22"  "SPCS1"   "SPTBN4"  "ADAM22" 

overlap_rep <- NULL
for (name in list_Gandal_DTU) {
  if (name %in% DTU_non_DE_IFN_rep){
    overlap_rep <- c(overlap_rep, name)
  }
}
overlap_rep_Gandal <- overlap_rep
####

#import spreadsheet with DGE and convert to df
table_2 <- read_excel("~/Desktop/MRes/20220224/NIHMS1007359-supplement-Table_S1.xlsx", sheet ="DGE", col_names =TRUE)
df_DGE_Gandal <-as.data.frame.matrix(table_2, col.names = colnames(table_2))

#only keep relevant columns
df_DGE_Gandal_ASD <- df_DGE_Gandal[,c(1,2,8,14,18,19)]
View(df_DGE_Gandal_ASD)
#order the data by log2FC
df_DGE_Gandal_ASD <- df_DGE_Gandal_ASD[order(abs(df_DGE_Gandal_ASD$ASD.log2FC),na.last=TRUE, decreasing=TRUE),]

#keep only isoforms when pvalue <0.05
df_DGE_Gandal_ASD <- df_DGE_Gandal_ASD[which(df_DGE_Gandal_ASD$ASD.fdr < 0.05),]

# Finds rows in DTU whose ensembl_gene_id is NOT in df_DGE_Gandal_ASD$ensembl_gene_id
Gandal_DTU_not_DGE <- !(df_DTU_Gandal_ASD$ensembl_gene_id[[1]] %in%  df_DGE_Gandal_ASD$ensembl_gene_id[[1]] )
Gandal_DTU_not_DGE

df_Gandal_DTU_not_DGE <-NULL

df_Gandal_DTU_not_DGE <- data.frame(
  ensembl_gene_id <- df_DTU_Gandal_ASD$ensembl_gene_id[[1]][Gandal_DTU_not_DGE],
  external_gene_id <- df_DTU_Gandal_ASD$external_gene_id[[1]][Gandal_DTU_not_DGE],
  ensembl_transcript_id <-df_DTU_Gandal_ASD$ensembl_transcript_id[[1]][Gandal_DTU_not_DGE],
  external_transcript_id <- df_DTU_Gandal_ASD$external_transcript_id[[1]][Gandal_DTU_not_DGE],
  DTU.ASD.Value <- df_DTU_Gandal_ASD$DTU.ASD.Value[[1]][Gandal_DTU_not_DGE],
  DTU.ASD.FDR <- df_DTU_Gandal_ASD$DTU.ASD.FDR[[1]][Gandal_DTU_not_DGE]
)

colnames(df_Gandal_DTU_not_DGE) <- c("ensembl_gene_id", "external_gene_id", "ensembl_transcript_id", "external_transcript_id", "DTU.ASD.Value", "DTU.ASD.FDR") 
View(df_Gandal_DTU_not_DGE)

#list of genes DTU non DGE in decreasing order of effect size
list_Gandal_DTUnonDGE <- df_Gandal_DTU_not_DGE$external_gene_id
list_Gandal_DTUnonDGE 

#finding overlap between the 2 lists
overlap <- NULL
for (name in DTU_non_DE_IFN_NPC) {
  if (name %in% list_Gandal_DTUnonDGE){
    overlap <- c(overlap, name)
  }
}

overlap
#"IL1RAP" "ISG20"  "VMP1"   "SPHK2"

#read the table from Haney's paper, only the sheet with DTE
table_3 <- read_excel("~/Desktop/MRes/20220313/Haney_table3.xlsx", sheet ="DEIsoform_Statistics", col_names =TRUE)
df_DTE_Haney <-as.data.frame(table_3, col.names = colnames(table_3))
#only keep relevant columns
df_DTE_Haney_2 <- df_DTE_Haney[,c(1,2,3,10,11,12)]
#keep only isoforms when FDR <0.05
df_DTE_Haney_3 <- df_DTE_Haney_2[which(df_DTE_Haney_2$WholeCortex_ASD_FDR < 0.05),]
#order the data by beta-value (effect size) in ASD
df_DTE_Haney_3 <- df_DTE_Haney_3[order(abs(df_DTE_Haney_3$WholeCortex_ASD_logFC),na.last=TRUE, decreasing=TRUE),]

list_Haney_DTU <-df_DTE_Haney_3$external_gene_name

#import spreadsheet with DGE
table_4 <- read_excel("~/Desktop/MRes/20220313/Haney_table3.xlsx", sheet ="DEGene_Statistics", col_names =TRUE)

#convert to df
df_DGE_Haney <-as.data.frame.matrix(table_4, col.names = colnames(table_4))
#only keep relevant columns
df_DGE_Haney_2 <- df_DGE_Haney[,c(1,2,10,11)]
#keep only isoforms when FDR <0.05
df_DGE_Haney_3 <- df_DGE_Haney_2[which(df_DGE_Haney_2$WholeCortex_ASD_FDR < 0.05),]
#order the data by log2FC
df_DGE_Haney_3 <- df_DGE_Haney_3[order(abs(df_DGE_Haney_3$WholeCortex_ASD_logFC),na.last=TRUE, decreasing=TRUE),]


# Finds rows in DTU whose ensembl_gene_id is NOT in DGE
df_Haney_DTE_not_DGE <- !(df_DTE_Haney_3$ensembl_gene_id[[1]] %in%  df_DGE_Haney_3$ensembl_gene_id[[1]] )

df_Haney_DTE_not_DGE <- data.frame(
  ensembl_gene_id <- df_DTE_Haney_3$ensembl_gene_id[[1]][df_Haney_DTE_not_DGE],
  external_gene_name <- df_DTE_Haney_3$external_gene_name[[1]][df_Haney_DTE_not_DGE],
  ensembl_transcript_id <-df_DTE_Haney_3$ensembl_transcript_id[[1]][df_Haney_DTE_not_DGE],
  description <- df_DTE_Haney_3$description[[1]][df_Haney_DTE_not_DGE],
  WholeCortex_ASD_logFC <- df_DTE_Haney_3$WholeCortex_ASD_logFC[[1]][df_Haney_DTE_not_DGE],
  WholeCortex_ASD_FDR <- df_DTE_Haney_3$WholeCortex_ASD_FDR[[1]][df_Haney_DTE_not_DGE]
)

colnames(df_Haney_DTE_not_DGE) <- c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id", "description",  "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR") 

#list of genes DTU non DGE in decreasing order of effect size
list_Haney_DTE_not_DGE <- df_Haney_DTE_not_DGE$external_gene_name

#finding overlap between the 2 lists
overlap2 <- NULL
for (name in DTU_non_DE_IFN_NPC) {
  if (name %in% list_Haney_DTE_not_DGE){
    overlap2 <- c(overlap2, name)
  }
}

overlap2
#[1] "CCSER1"     "ID1"        "IQSEC2"     "ZFR2"       "IL1RAP"     "PDZRN3"     "NRBP2"      "AUNIP"      "AUNIP"      "NRBP2"      "ISG20"      "KRTAP5-AS1" "CACNA1B"   
#[14] "SH2D5"      "NDRG1"      "CYP46A1"    "PPP1R10"    "HERPUD2"    "HERPUD2"  

#overlap between the 2 datasets from the literature, genes with DTU and no DGE
overlap3 <- NULL
for (name in list_Gandal_DTUnonDGE) {
  if (name %in% list_Haney_DTE_not_DGE){
    overlap3 <- c(overlap3, name)
  }
}
overlap3

DTU_only <- list(
  Gandal<- list_Gandal_DTUnonDGE,
  Haney <-list_Haney_DTE_not_DGE,
  INF_NPC<- DTU_non_DE_IFN_NPC
)


p <- ggVennDiagram(DTU_only, category.names = c("Gandal","Haney","IFNg-NPC"), label = "count")
p + scale_fill_distiller(palette = "Reds", direction = 1)

#information about overlaps
v.table <- venn(DTU_only)
print(v.table)

#find the characteristics of DTU in IFN dataset when DTU alos found in Gandal's OR Haney's dataset (always no DGE)
overlap_late <- NULL
for (name in list_Gandal_DTUnonDGE) {
  if (name %in% DTU_non_DE_IFN_Late){
    overlap_late <- c(overlap_late, name)
  }
}
for (name in list_Haney_DTE_not_DGE) {
  if (name %in% DTU_non_DE_IFN_Late){
    overlap_late <- c(overlap_late, name)
  }
}
overlap_late
overlap_late <-unique(overlap_late)

# MYO7A"    "POLD4"    "CNKSR2"   "CCNB2"    "PRKCQ"    "C21orf62" "HACL1"    "SPOCK3"   "GSN"      "MEF2C"    "PAK6"     "PACSIN1"  "ARNTL2"   "DCLRE1A"  "TULP3"   
#"LYSMD2"  

overlap_early <- NULL
for (name in list_Gandal_DTUnonDGE) {
  if (name %in% DTU_non_DE_IFN_early){
    overlap_early <- c(overlap_early, name)
  }
}
for (name in list_Haney_DTE_not_DGE) {
  if (name %in% DTU_non_DE_IFN_early){
    overlap_early <- c(overlap_early, name)
  }
}

overlap_early
overlap_early <-unique(overlap_early)
# "SLC30A6" "RPL22L1" "FAM13A"  "PLS3"    "SPCS1"   "SPTBN4"  "POLD4"   "SLC29A2" "NR4A1"   "COL26A1" "HACL1"   "SPOCK3"  "WFDC2"   "GSN"     "SORBS3"  "PPFIA4"  "PRDX4"  
# "TULP3"   "ZNF329"  "SLC22A5"

overlap_rep <- NULL
for (name in list_Gandal_DTUnonDGE) {
  if (name %in% DTU_non_DE_IFN_rep){
    overlap_rep <- c(overlap_rep, name)
  }
}

for (name in list_Haney_DTE_not_DGE) {
  if (name %in% DTU_non_DE_IFN_rep){
    overlap_rep <- c(overlap_rep, name)
  }
}
overlap_rep <- unique(overlap_rep)
overlap_rep
# "ACTN1"       "TIPIN"       "PLD2"        "TSPAN4"      "ADCYAP1"     "FBLN5"       "HELQ"        "KCNH3"       "HMGCL"       "SLC25A5-AS1" "SNHG3"       "SNCA"       
#[13] "KCNH1"       "SLC7A5"      "MAGI2"       "PIP4K2C"     "PCDHB9"      "CGRRF1"      "SNRNP25"     "MXI1"        "DNAJC21"     "ZNF689"

#############
#read the table from Okay's paper
Okay<- read_excel("~/Desktop/MRes/20220629/from_literature/Okay_DAS.xlsx", col_names =F)
df_Okay<-as.data.frame(Okay)
colnames(df_Okay) <-df_Okay[3,]
df_Okay2 <- df_Okay[4:186, c(1,2,21, 25)]
df_Okay2 <- df_Okay[4:186,]
colnames(df_Okay2) <- c("isoform", "geneID", "log2FC(parents/children)", "p-val")

ens <- df_Okay2$geneID
ens2 <- NULL
for (id in df_Okay2$geneID) {
  new_id <-gsub("\\..*","", id)
  ens2 <- c(ens2, new_id)
}
ens2 
df_Okay2$geneID <- ens2

df_Okay2$SYMBOL = mapIds(org.Hs.eg.db,
                         keys=df_Okay2$geneID, 
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")

df_Okay2 <-df_Okay2[,c(1,2,5,3,4)]
list_Okay_DTU <-df_Okay2$SYMBOL
list_Okay_DTU 

#look for a particular gene in the list
gene_studied <-"SAMD4"
line_gene <-grep(gene_studied, df_Okay2$SYMBOL)
line_gene

df_Okay2[53, 18]
#CLCN4
#[1] "ENST00000421085.6;ENST00000380833.9;ENST00000380829.5" "ENST00000421085.6;ENST00000380833.9;ENST00000380829.5"
#MYO15B
#[1] "ENST00000621743.4;ENST00000610510.4;ENST00000645453.2" "ENST00000621743.4;ENST00000610510.4;ENST00000645453.2" "ENST00000582561.2;ENST00000584516.5"
#NR4A1
#"ENST00000478250.1;ENST00000550557.1;ENST00000564201.1"
#TBC1D4
#[1] "ENST00000493487.1;ENST00000648194.1;ENST00000377625.6;ENST00000377636.7;ENST00000431480.6"
#[2] "ENST00000493487.1;ENST00000648194.1;ENST00000377625.6;ENST00000377636.7;ENST00000431480.6"
#SAMD4
#"ENST00000251091.9;ENST00000631086.2;ENST00000554335.5;ENST00000555091.1;ENST00000557013.1;ENST00000555112.1;ENST00000557692.2;ENST00000392067.7"

#finding overlap between Okay DTU et DTU IFN NPCs
overlap <- NULL
for (name in DTU_non_DE_IFN_NPC) {
  if (name %in% list_Okay_DTU ){
    overlap <- c(overlap, name)
  }
}
overlap_NPC_Okay <-overlap  

#finding overlap between Gandal DTU et DTU IFN neurons
overlap_late <- NULL
for (name in list_Okay_DTU) {
  if (name %in% DTU_non_DE_IFN_Late){
    overlap_late <- c(overlap_late, name)
  }
}
overlap_late_Okay <- overlap_late

overlap_early <- NULL
for (name in list_Okay_DTU ) {
  if (name %in% DTU_non_DE_IFN_early){
    overlap_early <- c(overlap_early, name)
  }
}
overlap_early_Okay <- overlap_early

overlap_rep <- NULL
for (name in list_Okay_DTU ) {
  if (name %in% DTU_non_DE_IFN_rep){
    overlap_rep <- c(overlap_rep, name)
  }
}
overlap_rep_Okay <-overlap_rep

#Voineagu's dataset
Voineagu <- read_excel("~/Desktop/MRes/20220629/from_literature/DAS_ASD_CTL_Voineagu.xlsx", col_names =F)
df_Voineagu <-as.data.frame(Voineagu)

colnames (df_Voineagu) <- c("gene_id", "gene_symbol", "chromosome", "C1", "A", "C2", "%inc difference")
df_Voineagu <- df_Voineagu[4:214,]
list_Voineagu_DTU <- as.list(df_Voineagu$gene_symbol)

#finding overlap between Voineagu DTU et DTU IFN NPCs
overlap <- NULL
for (name in DTU_non_DE_IFN_NPC) {
  if (name %in% list_Voineagu_DTU ){
    overlap <- c(overlap, name)
  }
}
overlap_NPC_Voineagu <- overlap  

#finding overlap between Voineagu DTU et DTU IFN neurons
overlap_late <- NULL
for (name in list_Voineagu_DTU) {
  if (name %in% DTU_non_DE_IFN_Late){
    overlap_late <- c(overlap_late, name)
  }
}
overlap_late_Voineagu <- overlap_late


overlap_early <- NULL
for (name in list_Voineagu_DTU) {
  if (name %in% DTU_non_DE_IFN_early){
    overlap_early <- c(overlap_early, name)
  }
}
overlap_early_Voineagu <- overlap_early

overlap_rep <- NULL
for (name in list_Voineagu_DTU) {
  if (name %in% DTU_non_DE_IFN_rep){
    overlap_rep <- c(overlap_rep, name)
  }
}
overlap_rep_Voineagu <- overlap_rep

#Stamova, from paper
list_Stamova_DTU <- c("CDK13", "USP48", "SFPQ", "FXR1", "C19orf2", 'TARS2', 'LRPPRC', 'PIK3C3','CLTB', 'SOD2', 'OS9', 'ACPT', 'PPP2R2A', 'C14orf159', 'FGR', 'GSN', 'EAPP',
                      'PIP4K2A', 'TADA3', 'PRSS36', 'HELQ', 'EMD', 'C1orf175', 'AEBP2', 'R3HDM1', 'SRPK1', 'LEF1', 'MPHOSPH10', 'PRDX1', 'MPP1', 'CNOT2', 'GOLGA7', 'WDR67',
                      'AAMP', 'KLHL9', 'CC2D1A', 'STAT4', 'DHX29', 'MGST3', 'TEPP', 'UTRN', 'PUM2', 'CHID1', 'GFER', 'RPGR', 'SUCLA2', 'ZNF512B', 'MORN2', 'DNAJC17', 'FGD3')


#finding overlap between Stamova DTU et DTU IFN NPCs
overlap <- NULL
for (name in DTU_non_DE_IFN_NPC) {
  if (name %in% list_Stamova_DTU){
    overlap <- c(overlap, name)
  }
}
overlap_NPC_Stamova <- overlap  

#finding overlap between Stamova DTU et DTU IFN neurons
overlap_late <- NULL
for (name in list_Stamova_DTU) {
  if (name %in% DTU_non_DE_IFN_Late){
    overlap_late <- c(overlap_late, name)
  }
}
overlap_late_Stamova <- overlap_late


overlap_early <- NULL
for (name in list_Stamova_DTU) {
  if (name %in% DTU_non_DE_IFN_early){
    overlap_early <- c(overlap_early, name)
  }
}
overlap_early_Stamova  <- overlap_early

overlap_rep <- NULL
for (name in list_Stamova_DTU) {
  if (name %in% DTU_non_DE_IFN_rep){
    overlap_rep <- c(overlap_rep, name)
  }
}
overlap_rep_Stamova  <- overlap_rep

####  ASD genes DTU from the literature- overlap with NPCs IFN
#Switchlist from DTU analysis for IFN
NPC_DTU_ASD_lit <- c(overlap_NPC_Gandal, overlap_NPC_Okay, overlap_NPC_Stamova, overlap_NPC_Voineagu)
DTU_IFN_NPC_features <-SwitchList_NPCs$isoformFeatures
ASD_lit_genes_IFN <- subset(DTU_IFN_NPC_features, gene_id %in% NPC_DTU_ASD_lit)
View(ASD_lit_genes_IFN)
ASD_lit_genes_IFN2 <-ASD_lit_genes_IFN[, c(3,4,8,9,27,28,32, 37)]
View(ASD_lit_genes_IFN2)

DTU_IFN_NPC_altern_splicing <- SwitchList_NPCs$AlternativeSplicingAnalysis
View(DTU_IFN_NPC_altern_splicing)

altern_NPC_IFN_ASD_lit <- subset(DTU_IFN_NPC_altern_splicing, DTU_IFN_NPC_altern_splicing$isoform_id %in% ASD_lit_genes_IFN2$isoform_id)
View(altern_NPC_IFN_ASD_lit)
altern_NPC_IFN_ASD_lit$gene_id <-ASD_lit_genes_IFN2$gene_id
altern_NPC_IFN_ASD_lit <- altern_NPC_IFN_ASD_lit[, c(26, 1, 2, 5, 8, 11, 14, 17, 20, 23)]

altern_NPC_IFN_ASD_lit$iso_biotype <- ASD_lit_genes_IFN2$iso_biotype
altern_NPC_IFN_ASD_lit$dIF <- ASD_lit_genes_IFN2$dIF
altern_NPC_IFN_ASD_lit$isoform_switch_q_value <- ASD_lit_genes_IFN2$isoform_switch_q_value
altern_NPC_IFN_ASD_lit$switchConsequencesGene <- ASD_lit_genes_IFN2$switchConsequencesGene

####  ASD genes DTU from the literature- overlap with neurons late treatment IFN
##Switchlist from DTU analysis for IFN
late_DTU_ASD_lit <- c(overlap_late_Gandal, overlap_late_Okay, overlap_late_Stamova, overlap_late_Voineagu)
DTU_IFN_late_features <-SwitchList_neurons_late$isoformFeatures
ASD_lit_genes_IFN_late <- subset(DTU_IFN_late_features, gene_id %in% late_DTU_ASD_lit)
ASD_lit_genes_IFN_late2 <-ASD_lit_genes_IFN_late[, c(3,4,8,9,27,28,32, 37)]
View(ASD_lit_genes_IFN_late2)

altern_late_IFN_ASD_lit <- subset(DTU_IFN_late_altern_splicing, isoform_id %in% ASD_lit_genes_IFN_late2$isoform_id)
altern_late_IFN_ASD_lit$gene_id <-ASD_lit_genes_IFN_late2$gene_id
altern_late_IFN_ASD_lit <- altern_late_IFN_ASD_lit[, c(26, 1, 2, 5, 8, 11, 14, 17, 20, 23)]

altern_late_IFN_ASD_lit$iso_biotype <- ASD_lit_genes_IFN_late2$iso_biotype
altern_late_IFN_ASD_lit$dIF <- ASD_lit_genes_IFN_late2$dIF
altern_late_IFN_ASD_lit$isoform_switch_q_value <- ASD_lit_genes_IFN_late2$isoform_switch_q_value
altern_late_IFN_ASD_lit$switchConsequencesGene <- ASD_lit_genes_IFN_late2$switchConsequencesGene

####  ASD genes DTU from the literature- overlap with neurons early treatment IFN
##Switchlist from DTU analysis for IFN
early_DTU_ASD_lit <- c(overlap_early_Gandal, overlap_early_Okay, overlap_early_Stamova, overlap_early_Voineagu)
DTU_IFN_early_features <-SwitchList_neurons_early$isoformFeatures
ASD_lit_genes_IFN_early <- subset(DTU_IFN_early_features, gene_id %in% early_DTU_ASD_lit)
View(ASD_lit_genes_IFN_early)
ASD_lit_genes_IFN_early2 <-ASD_lit_genes_IFN_early[, c(3,4,8,9,27,28,32, 37)]
View(ASD_lit_genes_IFN_early2)

DTU_IFN_early_altern_splicing <- SwitchList_neurons_early$AlternativeSplicingAnalysis
View(DTU_IFN_early_altern_splicing)

altern_early_IFN_ASD_lit <- subset(DTU_IFN_early_altern_splicing, isoform_id %in% ASD_lit_genes_IFN_early2$isoform_id)
View(altern_early_IFN_ASD_lit)
altern_early_IFN_ASD_lit$gene_id <-ASD_lit_genes_IFN_early2$gene_id
altern_early_IFN_ASD_lit <- altern_early_IFN_ASD_lit[, c(26, 1, 2, 5, 8, 11, 14, 17, 20, 23)]

altern_early_IFN_ASD_lit$iso_biotype <- ASD_lit_genes_IFN_early2$iso_biotype
altern_early_IFN_ASD_lit$dIF <- ASD_lit_genes_IFN_early2$dIF
altern_early_IFN_ASD_lit$isoform_switch_q_value <- ASD_lit_genes_IFN_early2$isoform_switch_q_value
altern_early_IFN_ASD_lit$switchConsequencesGene <- ASD_lit_genes_IFN_early2$switchConsequencesGene
########

####  ASD genes DTU from the literature- overlap with neurons repeated treatment IFN
#l#Switchlist from DTU analysis for IFN
rep_DTU_ASD_lit <- c(overlap_rep_Gandal, overlap_rep_Okay, overlap_rep_Stamova, overlap_rep_Voineagu)
DTU_IFN_rep_features <-SwitchList_neurons_rep$isoformFeatures
ASD_lit_genes_IFN_rep <- subset(DTU_IFN_rep_features, gene_id %in% rep_DTU_ASD_lit)
View(ASD_lit_genes_IFN_rep)
ASD_lit_genes_IFN_rep2 <-ASD_lit_genes_IFN_rep[, c(3,4,8,9,27,28,32, 37)]
View(ASD_lit_genes_IFN_rep2)

DTU_IFN_rep_altern_splicing <- SwitchList_neurons_rep$AlternativeSplicingAnalysis
View(DTU_IFN_rep_altern_splicing)

altern_rep_IFN_ASD_lit <- subset(DTU_IFN_rep_altern_splicing, isoform_id %in% ASD_lit_genes_IFN_rep2$isoform_id)
View(altern_rep_IFN_ASD_lit)
altern_rep_IFN_ASD_lit$gene_id <-ASD_lit_genes_IFN_rep2$gene_id
altern_rep_IFN_ASD_lit <- altern_rep_IFN_ASD_lit[, c(26, 1, 2, 5, 8, 11, 14, 17, 20, 23)]


altern_rep_IFN_ASD_lit$iso_biotype <- ASD_lit_genes_IFN_rep2$iso_biotype
altern_rep_IFN_ASD_lit$dIF <- ASD_lit_genes_IFN_rep2$dIF
altern_rep_IFN_ASD_lit$isoform_switch_q_value <- ASD_lit_genes_IFN_rep2$isoform_switch_q_value
altern_rep_IFN_ASD_lit$switchConsequencesGene <- ASD_lit_genes_IFN_rep2$switchConsequencesGene


