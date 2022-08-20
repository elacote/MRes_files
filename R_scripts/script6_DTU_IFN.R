#BiocManager::install("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)

#BiocManager::install("clusterProfiler")
library("clusterProfiler")
library(enrichplot)

#install.packages("readxl")
library(readxl)

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

#BiocManager::install("GeneOverlap")
library("GeneOverlap")

### Import quantifications
kallistoQuant <- importIsoformExpression(
  parentDir = "~/Desktop/MRes/20210402/IFN_from_kallisto/")
 #result= list with count and abundance estimates for each isoform in each sample 
#=count and abundance matrix

#analysis for NPCs
myDesign2 <-data.frame(
  sampleID = c("D18_T_M1", "D18_T_M2", "D18_T_M3", "D18_U_M1", "D18_U_M2", "D18_U_M3"),
  condition = c(rep("T",3), rep("U",3)),
  cell_line = c("M1", "M2", "M3", "M1", "M2", "M3")
)
#untreated as reference
myDesign2$condition <-factor(myDesign2$condition, levels = c("U", "T"))

### Create switchAnalyzeRlist
NPC_IFN_SwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix         = myDesign2,
  isoformExonAnnoation = "~/Desktop/MRes/20210402/gencode.v34.annotation.gtf",
  isoformNtFasta       = "~/Desktop/MRes/20210402/gencode.v34.transcripts_norRNA.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

NPC_IFN_switchList1 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = NPC_IFN_SwitchList,
  pathToOutput = "~/Desktop/MRes/20220226/NPCs/",
  outputSequences      = TRUE,  
  prepareForWebServers = TRUE  
)

extractSwitchSummary(NPC_IFN_switchList1)

SwitchList_NPCs2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = NPC_IFN_switchList1, 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  codingCutoff = 0.725,
  removeNoncodinORFs = TRUE,
  pathToCPATresultFile      = "~/Desktop/MRes/20210516/NPCs/result_CPAT_NPCs.txt",
  pathToPFAMresultFile      = "~/Desktop/MRes/20210516/NPCs/Pfam_NPCs.txt",
  pathToIUPred2AresultFile  = "~/Desktop/MRes/20210516/NPCs/iuPred2A_NPCs.result",
  pathToSignalPresultFile   = "~/Desktop/MRes/20210516/NPCs/signalP_NPCs.txt",
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'IDR_identified',
    'IDR_type',
    'signal_peptide_identified'),
  pathToOutput = "~/Desktop/MRes/20220226/NPCs/",
  outputPlots               = TRUE,  # keeps the function from outputting the plots from this example
  quiet = FALSE
 )

SwitchList_NPCs <-SwitchList_NPCs2

extractTopSwitches(SwitchList_NPCs2, n=15)
#            gene_ref    gene_id  gene_name condition_1 condition_2 gene_switch_q_value switchConsequencesGene Rank
#1  geneComp_00019848     IL18BP     IL18BP           U           T        3.092802e-89                   TRUE    1
#2  geneComp_00019698      IFI27      IFI27           U           T        1.135668e-82                   TRUE    2
#3  geneComp_00027733     PSMB10     PSMB10           U           T        7.924566e-71                   TRUE    3
#4  geneComp_00035067      WARS1      WARS1           U           T        9.656119e-66                   TRUE    4
#5  geneComp_00007638       ADAR       ADAR           U           T        2.502528e-57                  FALSE    5
#6  geneComp_00028632       RMI2       RMI2           U           T        1.754526e-45                   TRUE    6
#7  geneComp_00010934      ALG13      ALG13           U           T        1.918284e-39                  FALSE    7
#8  geneComp_00033316     TIMM10     TIMM10           U           T        9.431934e-35                  FALSE    8
#9  geneComp_00019702      IFI35      IFI35           U           T        4.143783e-28                   TRUE    9
#10 geneComp_00020082      ISG15      ISG15           U           T        1.655137e-25                  FALSE   10
#11 geneComp_00011195     ANXA11     ANXA11           U           T        5.372066e-25                   TRUE   11
#12 geneComp_00027775      PSME2      PSME2           U           T        8.020882e-25                   TRUE   12
#13 geneComp_00032652      SUMF2      SUMF2           U           T        6.508695e-23                   TRUE   13
#14 geneComp_00000981 AC007326.5 AC007326.5           U           T        1.787101e-19                  FALSE   14
#15 geneComp_00018735      GSDMD      GSDMD           U           T        7.167628e-19                   TRUE   15

#switchplot for a given gene
gene_with_DTU <-"CLCN4"
switchPlot(SwitchList_NPCs, gene = gene_with_DTU, condition1 = "U", condition2 ="T")

#order both table in same order (isoform_id)
SwitchList_NPCs$isoformFeatures <- SwitchList_NPCs$isoformFeatures[order(SwitchList_NPCs$isoformFeatures$isoform_id),]
SwitchList_NPCs$isoformSwitchAnalysis <- SwitchList_NPCs$isoformSwitchAnalysis[order(SwitchList_NPCs$isoformSwitchAnalysis$isoform_id),]
SwitchList_NPCs$AlternativeSplicingAnalysis <- SwitchList_NPCs$AlternativeSplicingAnalysis[order(SwitchList_NPCs$AlternativeSplicingAnalysis$isoform_id),]

#keep only isoform with no dge from DESqe2
selected_isoforms3 <-SwitchList_NPCs$isoformFeatures$gene_id %in% non_DE_list_NPC$SYMBOL
SwitchList_NPCs$isoformFeatures <- SwitchList_NPCs$isoformFeatures[selected_isoforms3,]
SwitchList_NPCs$orfAnalysis <- SwitchList_NPCs$orfAnalysis[selected_isoforms3,]
SwitchList_NPCs$isoformSwitchAnalysis <- SwitchList_NPCs$isoformSwitchAnalysis[selected_isoforms3,]
SwitchList_NPCs$AlternativeSplicingAnalysis <- SwitchList_NPCs$AlternativeSplicingAnalysis[selected_isoforms3,]

#volcano plot
EnhancedVolcano(SwitchList_NPCs$isoformFeatures,
                lab = SwitchList_NPCs$isoformFeatures$gene_name,
                x = 'dIF',
                y = 'isoform_switch_q_value',
                title = 'NPCs, IFN vs untreated',
                xlim= c(-1,1),
                ylim = c(-1,50),
                xlab="dIF",
                ylab="-log10(isoform_switch_q_value)", 
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = T,
                legendLabels =c("NS", "dIF", "q-value", "dIF and q-value"))

#keep only the isoforms with padj <0.05
selected_isoforms <- SwitchList_NPCs$isoformFeatures$isoform_switch_q_value <0.05
SwitchList_NPCs$isoformFeatures <- SwitchList_NPCs$isoformFeatures[selected_isoforms,]
SwitchList_NPCs$orfAnalysis <- SwitchList_NPCs$orfAnalysis[selected_isoforms,]
SwitchList_NPCs$isoformSwitchAnalysis <- SwitchList_NPCs$isoformSwitchAnalysis[selected_isoforms,]
SwitchList_NPCs$AlternativeSplicingAnalysis <- SwitchList_NPCs$AlternativeSplicingAnalysis[selected_isoforms,]

View(SwitchList_NPCs$isoformFeatures)
View(SwitchList_NPCs$isoformSwitchAnalysis)
View(SwitchList_NPCs$AlternativeSplicingAnalysis)
list_gene <- unique(SwitchList_NPCs$isoformFeatures$gene_name)
length(list_gene) #129
length(SwitchList_NPCs$isoformFeatures$gene_name) #154

#order by abs(dIF)
SwitchList_NPCs$isoformSwitchAnalysis <- SwitchList_NPCs$isoformSwitchAnalysis[order(-abs(SwitchList_NPCs$isoformFeatures$dIF)),]
SwitchList_NPCs$AlternativeSplicingAnalysis <- SwitchList_NPCs$AlternativeSplicingAnalysis[order(-abs(SwitchList_NPCs$isoformFeatures$dIF)),]
SwitchList_NPCs$isoformFeatures <- SwitchList_NPCs$isoformFeatures[order(-abs(SwitchList_NPCs$isoformFeatures$dIF)),]
write.csv(SwitchList_NPCs$isoformFeatures, file= "SwitchList_NPCs$isoformFeatures_IFN20220702.csv", row.names = T)
write.csv(SwitchList_NPCs$AlternativeSplicingAnalysis, file= "SwitchList_NPCs$AlternativeSplicingAnalysis_IFN20220702.csv", row.names = T)
write.csv(SwitchList_NPCs$isoformSwitchAnalysis, file= "SwitchList_NPCs$isoformSwitchAnalysis_IFN20220702.csv", row.names = T)
writeLines(SwitchList_NPCs$isoformFeatures$gene_name, "DTUnonDE_IFN_NPC_genes_20220703.txt")

DTU_non_DE_IFN_NPC <- SwitchList_NPCs$isoformFeatures$gene_id

SwitchList_NPCs$isoformFeatures$gene_id
# "HIVEP3"      "ADGRF3"      "TMEM225B"    "RDM1P5"      "RDM1P5"      "FRMD3"       "LINC01431"   "SYS1-DBNDD2" "FSIP2-AS1"  
# "FSIP2-AS1"   "ZFAT"        "DYNC1I1"     "PAX8-AS1"    "ZBTB7C"      "DOC2GP"      "LMCD1-AS1"   "H4C15"       "LHX6"       
# "FRMD3"       "H4C15"       "MFSD6"       "NPEPL1"      "PAG1"        "SGK3"        "H2BC5"       "H2BC5"       "HEXIM2"     
#[28] "BCDIN3D-AS1" "CCSER1"      "SYS1-DBNDD2" "THORLNC"     "LINC01881"   "ZNF438"      "EDEM3"       "CDK6-AS1"    "GPR143"     
#[37] "PTPRD-AS1"   "CFAP299"     "ZNF213"      "GAREM1"      "HEXIM2"      "ID1"         "PRSS53"      "ZFAT"        "SNX10"      
#[46] "IQSEC2"      "SPHK1"       "HLA-DMA"     "TLCD1"       "IER3"        "IER3"        "SSH3"        "ZFR2"        "SUMF2"      
#[55] "CLCN4"       "CLCN4"       "PCDHA7"      "PCDHA7"      "EGLN3"       "CYP27A1"     "SUMF2"       "PPT2-EGFL8"  "IL1RAP"     
#[64] "ANXA11"      "THSD7B"      "WASHC2A"     "TRAM2-AS1"   "LRRC46"      "MYRF"        "ZNF815P"     "HTATIP2"     "PDZRN3"     
#[73] "ETFBKMT"     "ZNF473"      "BMPR1B"      "SNX10"       "UNC13A"      "GAREM1"      "TADA2A"      "NRBP2"       "ROCK2"      
#[82] "PDCD2L"      "SFXN4"       "MAN1C1"      "IGSF8"       "PUS10"       "AUNIP"       "AUNIP"       "ST6GAL1"     "ADGRF3"     
#[91] "FOLH1"       "CNTN2"       "C3orf33"     "AATK"        "DHRS11"      "NRBP2"       "ISG20"       "DDC"         "KRTAP5-AS1" 
#[100] "CACNA1B"     "MYO3A"       "QPCTL"       "SPINT2"      "BRF2"        "ALG13"       "SH2D5"       "ARPIN"       "ROCK2"      
#[109] "ZNF221"      "ST6GAL1"     "FAM110A"     "ING2"        "ING2"        "TMEM63B"     "ZNF189"      "NDRG1"       "FAM183A"    
#[118] "GABRD"       "ZBTB25"      "ZNF594"      "CYP46A1"     "DHRS4L2"     "ZNF32"       "FAM104A"     "H4C14"       "VMP1"       
#[127] "SPHK2"       "AS3MT"       "AS3MT"       "LIG4"        "ARPIN"       "COL18A1"     "MYO15B"      "MVP"         "ANKRD44"    
#[136] "DCC"         "SEMA5B"      "ADGRA2"      "KLHL20"      "RMI2"        "MAML1"       "ETNK2"       "CACNB4"      "TRMT5"      
#[145] "PPP1R10"     "CIART"       "RBM38"       "ZNF32"       "FDX2"        "AKIP1"       "HERPUD2"     "ZNF774"      "EGLN3"      
#[154] "HERPUD2"  

#look for a particular gene
gene_studied <-"ZNF774"
line_gene <-grep(gene_studied,SwitchList_NPCs$isoformFeatures$gene_id)
line_gene

#genes DTU with dIF>0.3
DTUnonDE_IFN_NPC_top <-subset(SwitchList_NPCs$isoformFeatures, abs(dIF) >0.3)
View(DTUnonDE_IFN_NPC_top)

#consequence summary
extractConsequenceSummary(SwitchList_NPCs, consequencesToAnalyze = 'all', plotGenes= FALSE, asFractionTotal = FALSE)

extractSplicingSummary(SwitchList_NPCs, returnResult = FALSE)
View(SwitchList_NPCs$AlternativeSplicingAnalysis)

splicingEnrichment <- extractSplicingEnrichment(
  SwitchList_NPCs,
  splicingToAnalyze='all',
  alpha = 0.05,
  returnResult=TRUE,
  returnSummary=TRUE
)
View(splicingEnrichment)

summary(SwitchList_NPCs)
#This switchAnalyzeRlist list contains:
#154 isoforms from 129 genes
#1 comparison from 2 conditions (in total 6 samples)

#Switching features:
#  Comparison Isoforms Switches Genes
#1     U vs T      152       23   129

#Feature analyzed:
#  [1] "Isoform Switch Identification, ORFs, ntSequence, aaSequence, Protein Domains, IDR, Signal Peptides, Alternative splicing, Switch Consequences, Coding Potential"

#number of isoform switches with predicted functional consequences
extractSwitchSummary(SwitchList_NPCs)
top_switches_NPCs <- extractTopSwitches(SwitchList_NPCs, filterForConsequences = TRUE, n = Inf, sortByQvals = FALSE)
genes_with_switches_NPCs <- top_switches_NPCs$gene_name
genes_with_switches_NPCs



extractSplicingGenomeWide(
  SwitchList_NPCs,
  featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
  splicingToAnalyze = "all", 
  plot=TRUE,
  returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
)

#GO enrichment
DTU_IFN_NPC_ENTREZID <-AnnotationDbi::select(org.Hs.eg.db, keys = SwitchList_NPCs$isoformFeatures$gene_id, columns = "ENTREZID", keytype = "SYMBOL")
ego2_NPC_IFN_DTU <- clusterProfiler::enrichGO(gene          = DTU_IFN_NPC_ENTREZID$ENTREZID,
                                        OrgDb         = OrgDb,
                                        ont           = "BP",
                                        pAdjustMethod = "fdr",
                                        pvalueCutoff  = 0.05,
                                        qvalueCutoff  = 0.01, 
                                        readable      = TRUE)

#no enriched term found...


#overlap between DTU NPCs and SFARI
#import spreadsheet with sfari genes
Sfari_genes <- read.csv("~/Desktop/MRes/20220611/SFARI.csv", header =TRUE, stringsAsFactors = FALSE)
Sfari_high_score <- subset(Sfari_genes , gene.score %in% c(1,2))
list_SFARI <- Sfari_high_score$gene.symbol

go.obj <-newGeneOverlap(list_SFARI, unique(SwitchList_NPCs$isoformFeatures$gene_id), 20068 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
#listA size=895
#listB size=129
#Intersection size=10
#Overlapping p-value=0.063
#Jaccard Index=0.0

overlap_DTU_NPC_sfari <- getIntersection(go.obj)
writeLines(overlap_DTU_NPC_sfari, "overlap_IFN_DTU_NPC_sfari20220703.txt")
overlap_DTU_NPC_sfari

#select only Sfari genes
ASD_genes_NPC_DTU <- subset(SwitchList_NPCs$isoformFeatures, SwitchList_NPCs$isoformFeatures$gene_id %in% c("HIVEP3", "CCSER1", "CLCN4", "PCDHA7", "DDC", "DCC", "UNC13A", "ZNF774", "CACNA1B", "IQSEC2"))

altern_NPC <- subset(SwitchList_NPCs$AlternativeSplicingAnalysis, isoform_id %in% ASD_genes_NPC_DTU$isoform_id)


#neurons effect of late treatment
myDesign3 <-data.frame(
  sampleID = c("D30_UT_M1", "D30_UT_M2", "D30_UT_M3", "D30_UU_M1", "D30_UU_M2", "D30_UU_M3"),
  condition = c(rep("UT",3), rep("UU",3)),
  cell_line = c("M1", "M2", "M3", "M1", "M2", "M3")
)
myDesign3$condition <-factor(myDesign3$condition, levels = c("UU", "UT"))

neurons_late_SwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix         = myDesign3,
  isoformExonAnnoation = "~/Desktop/MRes/20210402/gencode.v34.annotation.gtf",
  isoformNtFasta       = "~/Desktop/MRes/20210402/gencode.v34.transcripts_norRNA.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

################
switchList_neurons_IFN_late <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = neurons_late_SwitchList,
  pathToOutput = "~/Desktop/MRes/20220306/neurons_late_treatment/",
  outputSequences      = TRUE,  
  prepareForWebServers = TRUE  
)

SwitchList_neurons_late2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = switchList_neurons_IFN_late,
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  codingCutoff = 0.725,
  removeNoncodinORFs = TRUE,
  pathToCPATresultFile      = "~/Desktop/MRes/20210516/neurons_late_treatment/result_neurons_late.txt",
  pathToPFAMresultFile      = "~/Desktop/MRes/20210516/neurons_late_treatment/Pfam_neurons_late.txt",
  pathToIUPred2AresultFile  = "~/Desktop/MRes/20210516/neurons_late_treatment/IUPred2A_neurons_late.result",
  pathToSignalPresultFile   = "~/Desktop/MRes/20210516/neurons_late_treatment/signalP_neurons_late.txt",
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'IDR_identified',
    'IDR_type',
    'signal_peptide_identified'),
  pathToOutput = "~/Desktop/MRes/20210516/neurons_late_treatment/",
  outputPlots               = TRUE,  # keeps the function from outputting the plots from this example
  quiet = FALSE
)

SwitchList_neurons_late <-SwitchList_neurons_late2

View(SwitchList_neurons_late$isoformFeatures)
View(SwitchList_neurons_late$isoformSwitchAnalysis)
View(SwitchList_neurons_late$AlternativeSplicingAnalysis)

extractTopSwitches(SwitchList_neurons_late, n=15)
#gene_ref    gene_id  gene_name condition_1 condition_2 gene_switch_q_value switchConsequencesGene Rank
#1  geneComp_00036831      WARS1      WARS1          UT          UU       1.476203e-172                   TRUE    1
#2  geneComp_00029108     PSMB10     PSMB10          UT          UU       2.169318e-133                   TRUE    2
#3  geneComp_00013292        B2M        B2M          UT          UU        8.484244e-76                   TRUE    3
#4  geneComp_00029151      PSME2      PSME2          UT          UU        8.686170e-71                   TRUE    4
#5  geneComp_00032161      SBNO2      SBNO2          UT          UU        1.726667e-39                  FALSE    5
#6  geneComp_00028416      PMPCA      PMPCA          UT          UU        8.352120e-25                   TRUE    6
#7  geneComp_00020962      ISG15      ISG15          UT          UU        6.499145e-24                   TRUE    7
#8  geneComp_00025827        MX1        MX1          UT          UU        3.136665e-18                   TRUE    8
#9  geneComp_00024284      MEF2C      MEF2C          UT          UU        3.379033e-17                   TRUE    9
#10 geneComp_00019938     HDAC11     HDAC11          UT          UU        5.600162e-14                  FALSE   10
#11 geneComp_00029131      PSMD1      PSMD1          UT          UU        4.049834e-12                  FALSE   11
#12 geneComp_00021197     KCNIP4     KCNIP4          UT          UU        6.781066e-11                   TRUE   12
#13 geneComp_00001363 AC008443.5 AC008443.5          UT          UU        2.996377e-10                  FALSE   13
#14 geneComp_00013561     BMPR1A     BMPR1A          UT          UU        3.475063e-09                   TRUE   14
#15 geneComp_00030026       RMI2       RMI2          UT          UU        4.760260e-09                   TRUE   15

#order the tables in same order (isoform_id)
SwitchList_neurons_late$isoformFeatures <- SwitchList_neurons_late$isoformFeatures[order(SwitchList_neurons_late$isoformFeatures$isoform_id),]
SwitchList_neurons_late$isoformSwitchAnalysis <- SwitchList_neurons_late$isoformSwitchAnalysis[order(SwitchList_neurons_late$isoformSwitchAnalysis$isoform_id),]
SwitchList_neurons_late$AlternativeSplicingAnalysis <- SwitchList_neurons_late$AlternativeSplicingAnalysis[order(SwitchList_neurons_late$AlternativeSplicingAnalysis$isoform_id),]

#keep only isoform with no dge
non_DE_IFN_late <-nonDE_list_late_IFN$SYMBOL
selected_isoforms3 <-SwitchList_neurons_late$isoformFeatures$gene_id %in% non_DE_IFN_late
SwitchList_neurons_late$isoformFeatures <- SwitchList_neurons_late$isoformFeatures[selected_isoforms3,]
SwitchList_neurons_late$orfAnalysis <- SwitchList_neurons_late$orfAnalysis[selected_isoforms3,]
SwitchList_neurons_late$isoformSwitchAnalysis <- SwitchList_neurons_late$isoformSwitchAnalysis[selected_isoforms3,]
SwitchList_neurons_late$AlternativeSplicingAnalysis <- SwitchList_neurons_late$AlternativeSplicingAnalysis[selected_isoforms3,]

#volcano plot
EnhancedVolcano(SwitchList_neurons_late$isoformFeatures,
                lab = SwitchList_neurons_late$isoformFeatures$gene_name,
                x = 'dIF',
                y = 'isoform_switch_q_value',
                title = 'neurons, IFN late treatment vs untreated',
                xlim= c(-1,1),
                ylim = c(-0.5,25),
                xlab="dIF",
                ylab="-log10(isoform_switch_q_value)", 
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = T,
                legendLabels =c("NS", "dIF", "q-value", "dIF and q-value")
)

#keep only the isoforms with padj <0.05
selected_isoforms <- SwitchList_neurons_late$isoformFeatures$isoform_switch_q_value <0.05
SwitchList_neurons_late$isoformFeatures <- SwitchList_neurons_late$isoformFeatures[selected_isoforms,]
SwitchList_neurons_late$orfAnalysis <- SwitchList_neurons_late$orfAnalysis[selected_isoforms,]
SwitchList_neurons_late$isoformSwitchAnalysis <- SwitchList_neurons_late$isoformSwitchAnalysis[selected_isoforms,]
SwitchList_neurons_late$AlternativeSplicingAnalysis <- SwitchList_neurons_late$AlternativeSplicingAnalysis[selected_isoforms,]

View(SwitchList_neurons_late)
View(SwitchList_neurons_late$isoformFeatures)
View(SwitchList_neurons_late$isoformSwitchAnalysis)
View(SwitchList_neurons_late$AlternativeSplicingAnalysis)

list_gene <- unique(SwitchList_neurons_late$isoformFeatures$gene_name)
length(list_gene) #127
length(SwitchList_neurons_late$isoformFeatures$gene_name) #158

SwitchList_neurons_late$isoformSwitchAnalysis <- SwitchList_neurons_late$isoformSwitchAnalysis[order(-abs(SwitchList_neurons_late$isoformFeatures$dIF)),]
SwitchList_neurons_late$AlternativeSplicingAnalysis <- SwitchList_neurons_late$AlternativeSplicingAnalysis[order(-abs(SwitchList_neurons_late$isoformFeatures$dIF)),]
SwitchList_neurons_late$isoformFeatures <- SwitchList_neurons_late$isoformFeatures[order(-abs(SwitchList_neurons_late$isoformFeatures$dIF)),]
write.csv(SwitchList_neurons_late$isoformFeatures, file= "SwitchList_neurons_late$isoformFeatures_IFN20220702.csv", row.names = T)
write.csv(SwitchList_neurons_late$AlternativeSplicingAnalysis, file= "SwitchList_neurons_late$AlternativeSplicingAnalysis_IFN20220702.csv", row.names = T)
write.csv(SwitchList_neurons_late$isoformSwitchAnalysis, file= "SwitchList_neurons_late$isoformSwitchAnalysis_IFN20220702.csv", row.names = T)
writeLines(SwitchList_neurons_late$isoformFeatures$gene_name, "DTUnonDE_IFN_late_genes_20220703.txt")

DTU_non_DE_IFN_Late <- SwitchList_neurons_late$isoformFeatures$gene_id
#  [1] "VPS33B-DT"   "TPBG"        "CNKSR1"      "RMI2"        "VPS33B-DT"   "RIOX2"       "GPR143"      "C1orf53"     "MMD2"       
#[10] "CNKSR1"      "UNC5C"       "RMI2"        "SLC10A3"     "CADPS2"      "C1orf53"     "VTN"         "RAD54L"      "TEKT3"      
#[19] "CASTOR1"     "GRIK1"       "UBAC2-AS1"   "B3GALNT1"    "ABHD3"       "LYSMD2"      "P3H2"        "MYO7A"       "MRPL45P2"   
#[28] "UCA1"        "C15orf39"    "HDAC11"      "CTRL"        "TMEM250"     "COL25A1"     "LANCL1"      "ARNTL2"      "FOXO4"      
#[37] "CYB5RL"      "LANCL1"      "STON1"       "PRR7-AS1"    "IDNK"        "POLD4"       "TBC1D31"     "RUNDC3A-AS1" "CDK2AP1"    
#[46] "LINC00910"   "LUZP2"       "DCLRE1A"     "LINC01586"   "VTN"         "DCLRE1A"     "TEDC1"       "FAM117A"     "CASC6"      
#[55] "MRPL45P2"    "CEP112"      "BMPR1A"      "LINC01114"   "FDX2"        "SPOCK3"      "EEF2KMT"     "GSN"         "PRKCQ"      
#[64] "SLC25A20"    "CCNA1"       "PROCR"       "SAMD3"       "GPSM3"       "CCNA1"       "KLHL32"      "CASP6"       "SAC3D1"     
#[73] "GCLM"        "C21orf62"    "ZCCHC7"      "CCNB2"       "FERMT3"      "PRANCR"      "CYSLTR2"     "COQ10A"      "MED29"      
#[82] "P4HA1"       "TFAP4"       "SPATA6L"     "ZNF580"      "BCL2L11"     "COX16"       "ZNF764"      "FILIP1L"     "ZNF175"     
#[91] "B3GALNT1"    "ZNF175"      "MRPL43"      "GPR61"       "C5"          "ZNF45"       "TULP3"       "TMEM150C"    "SAC3D1"     
#[100] "PMPCA"       "KCNB1"       "KYAT3"       "DPM3"        "PHF19"       "HCG17"       "POMGNT2"     "POMGNT2"     "HCG17"      
#[109] "INSYN1"      "SCNM1"       "PAK6"        "SEMA4G"      "CLTCL1"      "FECH"        "ZNF580"      "CNKSR2"      "MAML1"      
#[118] "PHF3"        "B4GALT4"     "DPM3"        "SIRT5"       "DOK1"        "ARHGAP4"     "ZNF764"      "CHCHD4"      "SEMA3B"     
#[127] "TENM3-AS1"   "ABHD18"      "PMPCA"       "PSMD1"       "C12orf75"    "MEF2C"       "LINC02134"   "PSMD1"       "KIAA1549"   
#[136] "KIAA1549"    "LUZP2"       "LINC01586"   "PIK3R4"      "ALKBH2"      "TP53INP1"    "SRPK1"       "CEP104"      "PACSIN1"    
#[145] "TOGARAM1"    "PHF19"       "HACL1"       "ACLY"        "NECAB3"      "TTTY14"      "POLD4"       "TTTY14"      "ARHGAP4"    
#[154] "TTTY14"      "BMPR1A"      "NECAB3"      "HDAC11"      "PMPCA"   

#only genes with dIF >0.3
DTUnonDE_IFN_late_top <-subset(SwitchList_neurons_late$isoformFeatures, abs(dIF) >0.3)
View(DTUnonDE_IFN_late_top)
unique(DTUnonDE_IFN_late_top$gene_name)

#switchplot for one given gene
gene_of_interest <- 'CADPS2'
switchPlot(SwitchList_neurons_late2, gene = gene_of_interest, condition1 = "UU", condition2 ="UT")

##overlap between DTU late  and SFARI
go.obj <-newGeneOverlap(list_SFARI, unique(SwitchList_neurons_late$isoformFeatures$gene_id), 20576 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
#listA size=895
#listB size=127
#Intersection size=7
#Overlapping p-value=0.32
#Jaccard Index=0.0

overlap_DTU_late_sfari <- getIntersection(go.obj)
writeLines(overlap_DTU_late_sfari, "overlap_IFN_DTU_late_sfari20220703.txt")
overlap_DTU_late_sfari

DTU_IFN_late_altern_splicing <- SwitchList_neurons_late$AlternativeSplicingAnalysis
DTU_IFN_late_features <-SwitchList_neurons_late$isoformFeatures
#subset: sfari genes
ASD_genes_IFN_late <- subset(DTU_IFN_late_features,  gene_id %in% c("CADPS2", "CLTCL1", "CNKSR2", "KCNB1", "MEF2C", "PHF3", "TBC1D31" ))
View(ASD_genes_IFN_late)

altern_late_IFN <-subset( DTU_IFN_late_altern_splicing, isoform_id %in% ASD_genes_IFN_late$isoform_id)
altern_late_IFN$gene_id <-ASD_genes_IFN_late$gene_id
altern_late_IFN$iso_biotype <-ASD_genes_IFN_late$iso_biotype
altern_late_IFN$dIF <-ASD_genes_IFN_late$dIF
altern_late_IFN$isoform_switch_q_value <-ASD_genes_IFN_late$isoform_switch_q_value
altern_late_IFN$switchConsequencesGene <-ASD_genes_IFN_late$switchConsequencesGene

altern_late_IFN <- altern_late_IFN[, c(26, 1, 27, 28, 29, 30, 2, 5, 8, 11, 14, 17, 20, 23)]

# global consequence analysis
extractConsequenceSummary(SwitchList_neurons_late)
extractConsequenceEnrichment(SwitchList_neurons_late)
extractConsequenceGenomeWide(SwitchList_neurons_late)

# global splicing analysis
extractSplicingSummary( SwitchList_neurons_late)
extractSplicingEnrichment(SwitchList_neurons_late)
extractSplicingGenomeWide(SwitchList_neurons_late)

#number of isoform switches with predicted functional consequences
extractSwitchSummary(SwitchList_neurons_late)
top_switches_late <- extractTopSwitches(SwitchList_neurons_late, filterForConsequences = TRUE, n = Inf, sortByQvals = FALSE)
genes_with_switches_late <- top_switches_late$gene_name
genes_with_switches_late

"RAD54L" %in% genes_with_switches_late

#GO enrichment
DTU_IFN_late_ENTREZID <-AnnotationDbi::select(org.Hs.eg.db, keys = SwitchList_neurons_late$isoformFeatures$gene_id, columns = "ENTREZID", keytype = "SYMBOL")
ego2_late_IFN_DTU <- clusterProfiler::enrichGO(gene          = DTU_IFN_late_ENTREZID$ENTREZID,
                                               OrgDb         = OrgDb,
                                               ont           = "BP",
                                               pAdjustMethod = "fdr",
                                               pvalueCutoff  = 0.05,
                                               qvalueCutoff  = 0.01, 
                                               readable      = TRUE)
#no enriched term found...


#neurons effect of early treatment
myDesign4 <-data.frame(
  sampleID = c("D30_TU_M1", "D30_TU_M2", "D30_TU_M3", "D30_UU_M1", "D30_UU_M2", "D30_UU_M3"),
  condition = c(rep("TU",3), rep("UU",3)),
  cell_line = c("M1", "M2", "M3", "M1", "M2", "M3")
)

myDesign4$condition <-factor(myDesign4$condition, levels = c("UU", "TU"))

neurons_early_SwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix         = myDesign4,
  isoformExonAnnoation = "~/Desktop/MRes/20210402/gencode.v34.annotation.gtf",
  isoformNtFasta       = "~/Desktop/MRes/20210402/gencode.v34.transcripts_norRNA.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

####
switchList_neurons_IFN_early <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = neurons_early_SwitchList,
  pathToOutput = "~/Desktop/MRes/20220306/neurons_early_treatment/",
  outputSequences      = TRUE, 
  prepareForWebServers = TRUE  
)

SwitchList_neurons_early2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = switchList_neurons_IFN_early,
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  codingCutoff = 0.725,
  removeNoncodinORFs = TRUE,
  pathToCPATresultFile      = "~/Desktop/MRes/20210516/neurons_early_treatment/result_neurons_early.txt",
  pathToPFAMresultFile      = "~/Desktop/MRes/20210516/neurons_early_treatment/Pfam_neurons_early.txt",
  pathToIUPred2AresultFile  = "~/Desktop/MRes/20210516/neurons_early_treatment/IUPred2A_neurons_early.result",
  pathToSignalPresultFile   = "~/Desktop/MRes/20210516/neurons_early_treatment/signalP_neurons_early.txt",
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'IDR_identified',
    'IDR_type',
    'signal_peptide_identified'),
  pathToOutput = "~/Desktop/MRes/20210516/neurons_early_treatment/",
  outputPlots               = TRUE,  # keeps the function from outputting the plots from this example
  quiet = FALSE
)

SwitchList_neurons_early <-SwitchList_neurons_early2

View(SwitchList_neurons_early$isoformFeatures)
View(SwitchList_neurons_early$isoformSwitchAnalysis)
View(SwitchList_neurons_early$AlternativeSplicingAnalysis)

#switchplot for one given gene
gene_with_DTU <-"TBC1D4"
switchPlot(SwitchList_neurons_early, gene = gene_with_DTU, condition1 = "UU", condition2 ="TU")

#order the tables in same order (isoform_id)
SwitchList_neurons_early$isoformFeatures <- SwitchList_neurons_early$isoformFeatures[order(SwitchList_neurons_early$isoformFeatures$isoform_id),]
SwitchList_neurons_early$isoformSwitchAnalysis <- SwitchList_neurons_early$isoformSwitchAnalysis[order(SwitchList_neurons_early$isoformSwitchAnalysis$isoform_id),]
SwitchList_neurons_early$AlternativeSplicingAnalysis <- SwitchList_neurons_early$AlternativeSplicingAnalysis[order(SwitchList_neurons_early$AlternativeSplicingAnalysis$isoform_id),]

#non dge (non DE list from DE IFN)
selected_isoforms6 <- SwitchList_neurons_early$isoformFeatures$gene_id %in% nonDE_list_early_IFN$SYMBOL
SwitchList_neurons_early$isoformFeatures <- SwitchList_neurons_early$isoformFeatures[selected_isoforms6,]
SwitchList_neurons_early$orfAnalysis <- SwitchList_neurons_early$orfAnalysis[selected_isoforms6,]
SwitchList_neurons_early$isoformSwitchAnalysis <- SwitchList_neurons_early$isoformSwitchAnalysis[selected_isoforms6,]
SwitchList_neurons_early$AlternativeSplicingAnalysis <- SwitchList_neurons_early$AlternativeSplicingAnalysis[selected_isoforms6,]

#volcano plot
EnhancedVolcano(SwitchList_neurons_early$isoformFeatures,
                lab = SwitchList_neurons_early$isoformFeatures$gene_name,
                x = 'dIF',
                y = 'isoform_switch_q_value',
                title = 'neurons, early IFN treatment vs no treatment',
                xlim= c(-0.8,0.8),
                ylim = c(-0.5,30),
                xlab="dIF",
                ylab="-log10(isoform_switch_q_value)", 
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = T,
                legendLabels =c("NS", "dIF", "q-value", "dIF and q-value")
)

#keep only the isoforms with padj <0.05
selected_isoforms4 <- SwitchList_neurons_early$isoformFeatures$isoform_switch_q_value <0.05
SwitchList_neurons_early$isoformFeatures <- SwitchList_neurons_early$isoformFeatures[selected_isoforms4,]
SwitchList_neurons_early$orfAnalysis <- SwitchList_neurons_early$orfAnalysis[selected_isoforms4,]
SwitchList_neurons_early$isoformSwitchAnalysis <- SwitchList_neurons_early$isoformSwitchAnalysis[selected_isoforms4,]
SwitchList_neurons_early$AlternativeSplicingAnalysis <- SwitchList_neurons_early$AlternativeSplicingAnalysis[selected_isoforms4,]

View(SwitchList_neurons_early$isoformFeatures)
View(SwitchList_neurons_early$isoformSwitchAnalysis)
View(SwitchList_neurons_early$AlternativeSplicingAnalysis)

SwitchList_neurons_early$isoformSwitchAnalysis <- SwitchList_neurons_early$isoformSwitchAnalysis[order(-abs(SwitchList_neurons_early$isoformFeatures$dIF)),]
SwitchList_neurons_early$AlternativeSplicingAnalysis <- SwitchList_neurons_early$AlternativeSplicingAnalysis[order(-abs(SwitchList_neurons_early$isoformFeatures$dIF)),]
SwitchList_neurons_early$isoformFeatures <- SwitchList_neurons_early$isoformFeatures[order(-abs(SwitchList_neurons_early$isoformFeatures$dIF)),]
write.csv(SwitchList_neurons_early$isoformFeatures, file= "SwitchList_neurons_early$isoformFeatures_IFN20220702.csv", row.names = T)
write.csv(SwitchList_neurons_early$AlternativeSplicingAnalysis, file= "SwitchList_neurons_early$AlternativeSplicingAnalysis_IFN20220702.csv", row.names = T)
write.csv(SwitchList_neurons_early$isoformSwitchAnalysis, file= "SwitchList_neurons_early$isoformSwitchAnalysis_IFN20220702.csv", row.names = T)
writeLines(SwitchList_neurons_early$isoformFeatures$gene_name, "DTUnonDE_IFN_early_genes_20220703.txt")

DTU_non_DE_IFN_early <- SwitchList_neurons_early$isoformFeatures$gene_id
#  [1] "CRYBG2"       "GPR143"       "SLC29A2"      "FAH"          "PDPN"         "COL18A1"      "CCND2-AS1"    "STON1"        "EFL1P1"      
#[10] "SLC29A2"      "CCNJL"        "SPTBN4"       "CNKSR1"       "CTRL"         "COL2A1"       "FAM13A"       "FBXO27"       "TNFRSF1A"    
#[19] "ETV4"         "TBC1D4"       "TENT5A"       "THNSL1"       "THNSL1"       "PLS3"         "NT5C3A"       "EZR"          "PLS3"        
#[28] "RNFT2"        "ATP23"        "MYO1C"        "PLS3"         "SLC30A6"      "GSN"          "C22orf23"     "RNFT2"        "DNAAF4-CCPG1"
#[37] "RASGEF1C"     "XBP1"         "DEPDC4"       "TMEM250"      "SPOCK3"       "MRPL23"       "STIL"         "SORBS3"       "RPL22L1"     
#[46] "TMEM225B"     "RAMP2"        "PPFIA4"       "CATSPERE"     "GMPPA"        "TIPARP"       "HACL1"        "PHLPP1"       "SLC22A5"     
#[55] "POLD4"        "GYG2P1"       "WWOX"         "WFDC2"        "LINC02263"    "C8orf33"      "FAM66B"       "STYX"         "SLC25A10"    
#[64] "LINC02263"    "ARHGEF25"     "CYB561D1"     "HSCB"         "SORBS3"       "FECH"         "PER3"         "TMEM164"      "TRAF3IP2"    
#[73] "LGALSL"       "NR4A1"        "ZEB1-AS1"     "RAB15"        "HMGN4"        "HMGN4"        "TULP3"        "NECTIN2"      "ZNF761"      
#[82] "SNRK"         "MAK16"        "SPCS1"        "MAML1"        "BRCA1"        "GMPPA"        "COL26A1"      "DAXX"         "NEK11"       
#[91] "ZNF596"       "CBFB"         "PIK3CD"       "DUS4L-BCAP29" "RAP1A"        "GPX3"         "SLC26A10"     "FAM13A"       "SAMD4A"      
#[100] "MEN1"         "SPRY1"        "PCDH19"       "C1S"          "ST3GAL4"      "MYLK-AS1"     "SNRK"         "GADD45G"      "ADAM22"      
#[109] "LPP"          "ZNF329"       "ALOX12-AS1"   "PHYH"         "PRDX4"        "PIK3R3"       "SHROOM4"      "TFAP4"        "ADAM22"      
#[118] "PRDX4"        "GADD45G"      "RHOC"         "TP53INP1"     "TP53INP1"     "HACL1"        "C8orf33" 
unique (DTU_non_DE_IFN_early ) #107

#number of isoform switches with predicted functional consequences
extractSwitchSummary(SwitchList_neurons_early)
top_switches_early <- extractTopSwitches(SwitchList_neurons_early, filterForConsequences = TRUE, n = Inf, sortByQvals = FALSE)
genes_with_switches_early <- top_switches_early$gene_name
unique(genes_with_switches_early)

#look for a gene in the list
list_gene <- unique(SwitchList_neurons_early$isoformFeatures$gene_name)
length(list_gene) #107
length(SwitchList_neurons_early$isoformFeatures$gene_name) #124

#genes DTU non DE with dIF>0.3
DTUnonDE_IFN_early_top <-subset(SwitchList_neurons_early$isoformFeatures, abs(dIF) >0.3)
View(DTUnonDE_IFN_early_top) #20
unique(DTUnonDE_IFN_early_top$gene_name) #19

View(SwitchList_neurons_early$AlternativeSplicingAnalysis)

# global consequence analysis
extractConsequenceSummary(SwitchList_neurons_early)
extractConsequenceEnrichment(SwitchList_neurons_early)
extractConsequenceGenomeWide(SwitchList_neurons_early)

# global splicing analysis
extractSplicingSummary(SwitchList_neurons_early)
extractSplicingEnrichment(SwitchList_neurons_early)
extractSplicingGenomeWide(SwitchList_neurons_early)

#number of isoform switches with predicted functional consequences
extractSwitchSummary(SwitchList_neurons_early)
top_switches_early <- extractTopSwitches(SwitchList_neurons_early, filterForConsequences = TRUE, n = Inf, sortByQvals = FALSE)
genes_with_switches_early <- top_switches_early$gene_name
genes_with_switches_early

"RAD54L" %in% genes_with_switches_early

#GO enrichment
DTU_IFN_early_ENTREZID <-AnnotationDbi::select(org.Hs.eg.db, keys = SwitchList_neurons_early$isoformFeatures$gene_id, columns = "ENTREZID", keytype = "SYMBOL")
ego2_early_IFN_DTU <- clusterProfiler::enrichGO(gene          = DTU_IFN_early_ENTREZID$ENTREZID,
                                                OrgDb         = OrgDb,
                                                ont           = "BP",
                                                pAdjustMethod = "fdr",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.01, 
                                                readable      = TRUE)
#no enriched term found...


#overlap between DTU early and SFARI
go.obj <-newGeneOverlap(list_SFARI, unique(SwitchList_neurons_early$isoformFeatures$gene_id), 20576 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
#listA size=895
#listB size=107
#Intersection size=2
#Overlapping p-value=0.95
#Jaccard Index=0.0

overlap_DTU_early_sfari <- getIntersection(go.obj)
writeLines(overlap_DTU_early_sfari, "overlap_IFN_DTU_early_sfari20220703.txt")
overlap_DTU_early_sfari

#only Sfari genes
DTU_IFN_early_features <-SwitchList_neurons_early$isoformFeatures
ASD_genes_IFN <- subset(DTU_IFN_early_features, gene_id %in% c("PCDH19", "WWOX" ))
ASD_genes_IFN2 <-ASD_genes_IFN[, c(3,4,8,9,27,28,37)]

DTU_IFN_early_altern_splicing <- SwitchList_neurons_early$AlternativeSplicingAnalysis
altern_early_IFN <- subset(DTU_IFN_early_altern_splicing, DTU_IFN_early_altern_splicing$isoform_id %in% ASD_genes_IFN2$isoform_id)
altern_early_IFN$gene_id <-ASD_genes_IFN2$gene_id
altern_early_IFN$iso_biotype <-ASD_genes_IFN2$iso_biotype
altern_early_IFN$dIF <-ASD_genes_IFN2$dIF
altern_early_IFN$isoform_switch_q_value <-ASD_genes_IFN2$isoform_switch_q_value
altern_early_IFN$switchConsequencesGene <-ASD_genes_IFN2$switchConsequencesGene
altern_early_IFN <- altern_early_IFN[, c(26, 1, 27, 28, 29, 30, 2, 5, 8, 11, 14, 17, 20, 23)]


#neurons effect of repeated treatment vs late
myDesign5 <-data.frame(
  sampleID = c("D30_TT_M1", "D30_TT_M2", "D30_TT_M3", "D30_UT_M1", "D30_UT_M2", "D30_UT_M3"),
  condition = c(rep("TT",3), rep("UT",3)),
  cell_line = c("M1", "M2", "M3", "M1", "M2", "M3")
)

myDesign5$condition <-factor(myDesign5$condition, levels = c("UT", "TT"))

neurons_rep_SwitchList <- importRdata(  
  isoformCountMatrix   = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix         = myDesign5,
  isoformExonAnnoation = "~/Desktop/MRes/20210402/gencode.v34.annotation.gtf",
  isoformNtFasta       = "~/Desktop/MRes/20210402/gencode.v34.transcripts_norRNA.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

#######
switchList3 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = neurons_rep_SwitchList,
  pathToOutput = "~/Desktop/MRes/20220306/neurons_rep_treatment/",
  outputSequences      = TRUE,  
  prepareForWebServers = TRUE
)

SwitchList_neurons_rep2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = switchList3,
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  codingCutoff = 0.725,
  removeNoncodinORFs = TRUE,
  pathToCPATresultFile      = "~/Desktop/MRes/20210516/neurons_repeated_treatment/result_neurons_repeated.txt",
  pathToPFAMresultFile      = "~/Desktop/MRes/20210516/neurons_repeated_treatment/Pfam_neurons_repeated.txt",
  pathToIUPred2AresultFile  = "~/Desktop/MRes/20210516/neurons_repeated_treatment/IUPred2A_neurons_repeated.result",
  pathToSignalPresultFile   = "~/Desktop/MRes/20210516/neurons_repeated_treatment/signalP_neurons_repeated.txt",
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'IDR_identified',
    'IDR_type',
    'signal_peptide_identified'),
  pathToOutput = "~/Desktop/MRes/20210516/neurons_repeated_treatment/",
  outputPlots               = TRUE,  # keeps the function from outputting the plots from this example
  quiet = FALSE
)

SwitchList_neurons_rep <-SwitchList_neurons_rep2

View(SwitchList_neurons_rep$isoformFeatures)   
View(SwitchList_neurons_rep$isoformSwitchAnalysis)
View(SwitchList_neurons_rep$AlternativeSplicingAnalysis)

#switchplot
gene_of_interest <- 'POLA2'
switchPlot(SwitchList_neurons_rep, gene = gene_of_interest, condition1 = "UT", condition2 ="TT")

#order the tables in same order (isoform_id) 
SwitchList_neurons_rep$isoformFeatures <- SwitchList_neurons_rep$isoformFeatures[order(SwitchList_neurons_rep$isoformFeatures$isoform_id),]
SwitchList_neurons_rep$isoformSwitchAnalysis <- SwitchList_neurons_rep$isoformSwitchAnalysis[order(SwitchList_neurons_rep$isoformSwitchAnalysis$isoform_id),]
SwitchList_neurons_rep$AlternativeSplicingAnalysis <- SwitchList_neurons_rep$AlternativeSplicingAnalysis[order(SwitchList_neurons_rep$AlternativeSplicingAnalysis$isoform_id),]

#non dge 
selected_isoforms9 <- SwitchList_neurons_rep$isoformFeatures$gene_id %in% nonDE_list_repvslate_IFN$SYMBOL
SwitchList_neurons_rep$isoformFeatures <- SwitchList_neurons_rep$isoformFeatures[selected_isoforms9,]
SwitchList_neurons_rep$orfAnalysis <- SwitchList_neurons_rep$orfAnalysis[selected_isoforms9,]
SwitchList_neurons_rep$isoformSwitchAnalysis <- SwitchList_neurons_rep$isoformSwitchAnalysis[selected_isoforms9,]
SwitchList_neurons_rep$AlternativeSplicingAnalysis <- SwitchList_neurons_rep$AlternativeSplicingAnalysis[selected_isoforms9,]

#volcano plot
EnhancedVolcano(SwitchList_neurons_rep$isoformFeatures,
                lab = SwitchList_neurons_rep$isoformFeatures$gene_name,
                x = 'dIF',
                y = 'isoform_switch_q_value',
                title = 'neurons, repeated IFN treatment vs late',
                xlim= c(-0.8,0.8),
                ylim = c(-0.5,40),
                xlab="dIF",
                ylab="-log10(isoform_switch_q_value)", 
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = T,
                legendLabels =c("NS", "dIF", "q-value", "dIF and q-value")
)

#keep only the isoforms with padj <0.05
selected_isoforms7 <- SwitchList_neurons_rep$isoformFeatures$isoform_switch_q_value <0.05
SwitchList_neurons_rep$isoformFeatures <- SwitchList_neurons_rep$isoformFeatures[selected_isoforms7,]
SwitchList_neurons_rep$orfAnalysis <- SwitchList_neurons_rep$orfAnalysis[selected_isoforms7,]
SwitchList_neurons_rep$isoformSwitchAnalysis <- SwitchList_neurons_rep$isoformSwitchAnalysis[selected_isoforms7,]
SwitchList_neurons_rep$AlternativeSplicingAnalysis <- SwitchList_neurons_rep$AlternativeSplicingAnalysis[selected_isoforms7,]

View(SwitchList_neurons_rep$isoformFeatures)
View(SwitchList_neurons_rep$isoformSwitchAnalysis)
View(SwitchList_neurons_rep$AlternativeSplicingAnalysis)


SwitchList_neurons_rep$isoformSwitchAnalysis <- SwitchList_neurons_rep$isoformSwitchAnalysis[order(-abs(SwitchList_neurons_rep$isoformFeatures$dIF)),]
SwitchList_neurons_rep$AlternativeSplicingAnalysis <- SwitchList_neurons_rep$AlternativeSplicingAnalysis[order(-abs(SwitchList_neurons_rep$isoformFeatures$dIF)),]
SwitchList_neurons_rep$isoformFeatures <- SwitchList_neurons_rep$isoformFeatures[order(-abs(SwitchList_neurons_rep$isoformFeatures$dIF)),]
write.csv(SwitchList_neurons_rep$isoformFeatures, file= "SwitchList_neurons_rep$isoformFeatures_IFN20220702.csv", row.names = T)
write.csv(SwitchList_neurons_rep$AlternativeSplicingAnalysis, file= "SwitchList_neurons_rep$AlternativeSplicingAnalysis_IFN20220702.csv", row.names = T)
write.csv(SwitchList_neurons_rep$isoformSwitchAnalysis, file= "SwitchList_neurons_rep$isoformSwitchAnalysis_IFN20220702.csv", row.names = T)
writeLines(SwitchList_neurons_rep$isoformFeatures$gene_name, "DTUnonDE_IFN_rep_genes_20220703.txt")

DTU_non_DE_IFN_rep <- SwitchList_neurons_rep$isoformFeatures$gene_id
#[1] "ADCYAP1"     "GBP4"        "SNRNP25"     "GJA1"        "HTRA1"       "HTRA1"       "ORAI1"       "ORAI1"       "SKAP2"      
#[10] "CLMN"        "PCDHGA4"     "PCDHGA4"     "SFXN3"       "SKAP2"       "SERPINH1"    "SLC25A5-AS1" "UBAC2-AS1"   "C7orf31"    
#[19] "SEMA4B"      "SFXN3"       "SNHG30"      "ANO4"        "SEMA4B"      "B3GALNT1"    "PCK2"        "KCNMB2-AS1"  "STIL"       
#[28] "GNG5"        "MAPKAPK3"    "WNT5B"       "CLIC1"       "PCYT1B"      "SEZ6"        "COQ10A"      "B4GALT4"     "KCNH3"      
#[37] "HDAC11"      "PRRG1"       "SERPINH1"    "FBLN5"       "GNG5"        "KCNH1"       "INTS2"       "HDHD5"       "H4C14"      
#[46] "H4C14"       "CCNB1"       "UBN1"        "NOCT"        "HSCB"        "NMRAL1"      "CYTOR"       "USP44"       "SYPL1"      
#[55] "ZFAT"        "CDK2AP1"     "IQCJ-SCHIP1" "LINC01114"   "NKRF"        "SYPL1"       "ATG4C"       "ACTN1"       "KLHL32"     
#[64] "PIP4K2C"     "HELQ"        "SLC39A11"    "NTAN1"       "SNHG3"       "FOXO3"       "RBM11"       "FOXO3"       "CNTNAP5"    
#[73] "HMGCL"       "SNCA"        "YARS2"       "YARS2"       "LIFR-AS1"    "PEX12"       "FNTB"        "DUS2"        "ST3GAL4"    
#[82] "KLF11"       "STX10"       "CRISPLD1"    "ASB13"       "ECD"         "PPP4R4"      "ARL8B"       "AGTRAP"      "TIPIN"      
#[91] "MXI1"        "MAGI2"       "DNMT1"       "SNHG3"       "NKAIN4"      "PCDHB9"      "NDUFA7"      "SLC7A5"      "GMDS"       
#[100] "HMGCL"       "WAS"         "NDUFA7"      "LRRC8B"      "SLC7A5"      "TRIM66"      "KLF11"       "KHK"         "PRDM10"     
#[109] "DMAC2"       "EXOC3L1"     "ZNF213"      "DNAJC3-DT"   "GCSH"        "GCH1"        "PLPP1"       "TMEM216"     "MYD88"      
#[118] "LINC00910"   "CPT1B"       "VASP"        "PSMB10"      "PLD2"        "CGRRF1"      "KCNAB1"      "MAN2B1"      "TSPAN4"     
#[127] "KCNH3"       "CCDC47"      "ZNF689"      "MAGI2"       "ARHGAP4"     "PAQR8"       "RGS8"        "DENND3"      "FAM185A"    
#[136] "RMND5A"      "REEP6"       "EXOC3L1"     "NOL12"       "CCDC47"      "DENND3"      "GCK"         "PPP2CB"      "ENKD1"      
#[145] "POLA2"       "ADAMTS6"     "PPP2CB"      "MYO15B"      "SNX21"       "DNAJC21"     "WDSUB1"      "DNAJC21"     "DENND3"     
#[154] "MXI1"        "ST3GAL4"     "PRDM10"      "DNMT1"

#genes DTU non DE with dIF >0.3
DTUnonDE_IFN_rep_top <-subset(SwitchList_neurons_rep$isoformFeatures, abs(dIF) >0.3)
View(DTUnonDE_IFN_rep_top) #25
unique(DTUnonDE_IFN_rep_top$gene_name) #19

DTU_IFN_rep_features <-SwitchList_neurons_rep$isoformFeatures

#number of isoform switches with predicted functional consequences
extractSwitchSummary(SwitchList_neurons_rep)
top_switches_rep <- extractTopSwitches(SwitchList_neurons_rep, filterForConsequences = TRUE, n = Inf, sortByQvals = FALSE)
genes_with_switches_rep <- top_switches_rep$gene_name
unique(genes_with_switches_rep)

# global consequence analysis
extractConsequenceSummary(SwitchList_neurons_rep)
extractConsequenceEnrichment(SwitchList_neurons_rep)
extractConsequenceGenomeWide(SwitchList_neurons_rep)

# global splicing analysis
extractSplicingSummary(SwitchList_neurons_rep)
extractSplicingEnrichment(SwitchList_neurons_rep)
extractSplicingGenomeWide(SwitchList_neurons_rep)

#number of isoform switches with predicted functional consequences
extractSwitchSummary(SwitchList_neurons_rep)
top_switches_rep <- extractTopSwitches(SwitchList_neurons_rep, filterForConsequences = TRUE, n = Inf, sortByQvals = FALSE)
genes_with_switches_rep <- top_switches_rep$gene_name
genes_with_switches_rep

"RAD54L" %in% genes_with_switches_rep

#GO enrichment
DTU_IFN_rep_ENTREZID <-AnnotationDbi::select(org.Hs.eg.db, keys = SwitchList_neurons_rep$isoformFeatures$gene_id, columns = "ENTREZID", keytype = "SYMBOL")
ego2_rep_IFN_DTU <- clusterProfiler::enrichGO(gene          = DTU_IFN_rep_ENTREZID$ENTREZID,
                                              OrgDb         = OrgDb,
                                              ont           = "BP",
                                              pAdjustMethod = "fdr",
                                              pvalueCutoff  = 0.05,
                                              qvalueCutoff  = 0.01, 
                                              readable      = TRUE)

#no enriched term found...

#overlap between DTU rep and SFARI
go.obj <-newGeneOverlap(list_SFARI, unique(SwitchList_neurons_rep$isoformFeatures$gene_id), 20576 )
go.obj <-testGeneOverlap(go.obj)
go.obj

#GeneOverlap object:
#listA size=895
#listB size=128
#Intersection size=4
#Overlapping p-value=0.81
#Jaccard Index=0.0

overlap_DTU_rep_sfari <- getIntersection(go.obj)
writeLines(overlap_DTU_rep_sfari, "overlap_IFN_DTU_rep_sfari20220703.txt")
overlap_DTU_rep_sfari

ASD_genes_IFN <- subset(DTU_IFN_rep_features, gene_id %in% c("CNTNAP5", "TSPAN4", "SLC7A5", "POLA2" ))
ASD_genes_IFN2 <-ASD_genes_IFN[, c(3,4,8,9,27,28,37)]

DTU_IFN_rep_altern_splicing <- SwitchList_neurons_rep$AlternativeSplicingAnalysis
altern_rep_IFN <- subset(DTU_IFN_rep_altern_splicing, DTU_IFN_rep_altern_splicing$isoform_id %in% ASD_genes_IFN2$isoform_id)
altern_rep_IFN $gene_id <-ASD_genes_IFN2$gene_id
altern_rep_IFN$iso_biotype <-ASD_genes_IFN2$iso_biotype
altern_rep_IFN$dIF <-ASD_genes_IFN2$dIF
altern_rep_IFN$isoform_switch_q_value <-ASD_genes_IFN2$isoform_switch_q_value
altern_rep_IFN$switchConsequencesGene <-ASD_genes_IFN2$switchConsequencesGene
altern_rep_IFN <- altern_rep_IFN[, c(26, 1, 27, 28, 29, 30, 2, 5, 8, 11, 14, 17, 20, 23)]


####
#DTU non DE neurons whatever IFN treatment
list_DTUnonDGE <- c(SwitchList_neurons_rep$isoformFeatures$gene_name,
                    SwitchList_neurons_early$isoformFeatures$gene_name, 
                    SwitchList_neurons_late$isoformFeatures$gene_name)

length(list_DTUnonDGE)
list_DTUnonDGE <-unique(list_DTUnonDGE)
#335

list_DTUnonDGE
writeLines(list_DTUnonDGE, "DTUnonDE_IFN_neurons_20220703.txt")

list_DTUnonDGE_ENTREZID <-AnnotationDbi::select(org.Hs.eg.db, keys = list_DTUnonDGE, columns = "ENTREZID", keytype = "SYMBOL")
writeLines(list_DTUnonDGE_ENTREZID$ENTREZID, "DTUnonDE_IFN_neurons_ENTREZID20220703.txt")

