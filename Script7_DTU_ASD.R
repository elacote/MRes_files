#BiocManager::install("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)

#BiocManager::install("clusterProfiler")
library("clusterProfiler")
library(enrichplot)

#install.packages("readxl")
library(readxl)

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#install.packages("gplots")
library(gplots)

#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

#BiocManager::install("GeneOverlap")
library("GeneOverlap")


########################
### Import quantifications
kallistoQuant_griesi <- importIsoformExpression(
  parentDir = "~/Desktop/MRes/20220412/kal_results/")

#design matrix
#NPCs
design_NPC_griesi <-data.frame(
  sampleID = colnames(kallistoQuant_griesi$abundance)[2:30],
  condition = c(rep("ASD",16), rep("CTL",13))
)
#CTL = reference
design_NPC_griesi$condition <-factor(design_NPC_griesi$condition, levels = c("CTL", "ASD"))

### Create switchAnalyzeRlist
NPC_griesi_SwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant_griesi$counts,
  isoformRepExpression = kallistoQuant_griesi$abundance,
  designMatrix         = design_NPC_griesi,
  isoformExonAnnoation = "~/Desktop/MRes/20210402/gencode.v34.annotation.gtf",
  isoformNtFasta       = "~/Desktop/MRes/20210402/gencode.v34.transcripts_norRNA.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

################
switchList_NPC_griesi <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = NPC_griesi_SwitchList ,
  pathToOutput = "~/Desktop/MRes/20220625/NPCs_griesi/",
  outputSequences      = TRUE,
  prepareForWebServers = TRUE
)

NPCSwitchListAnalyzed1 <- analyzeCPAT(
  switchAnalyzeRlist   = switchList_NPC_griesi,
  pathToCPATresultFile = "~/Desktop/MRes/20220625/NPCs_griesi/result-11.txt",
  codingCutoff = 0.725,
  removeNoncodinORFs   = TRUE   
)

NPCSwitchListAnalyzed1  <- analyzePFAM(
  switchAnalyzeRlist   = NPCSwitchListAnalyzed1,
  pathToPFAMresultFile = "~/Desktop/MRes/20220625/NPCs_griesi/NPC_griesi_Pfam.txt",
  showProgress=TRUE
)

#NPCSwitchListAnalyzed1  <- analyzeSignalP(switchAnalyzeRlist       = NPCSwitchListAnalyzed1,pathToSignalPresultFile  = "~/Desktop/MRes/20220625/NPCs_griesi/prediction_results.txt")
#  No signal peptides were found

#NPCSwitchListAnalyzed1 <- analyzeNetSurfP2(switchAnalyzeRlist        = NPCSwitchListAnalyzed1,pathToNetSurfP2resultFile  = "~/Desktop/MRes/20220625/NPCs_griesi/griesi_Net.csv")
#The 'pathToIUPred2AresultFile' does not appear to contain results for the isoforms stored in the switchAnalyzeRlist...


# analyse alternative splicing and intro retention
NPCSwitchListAnalyzed1<- analyzeAlternativeSplicing(NPCSwitchListAnalyzed1)
NPCSwitchListAnalyzed1 <- analyzeIntronRetention(NPCSwitchListAnalyzed1)

#volcano plot
EnhancedVolcano(NPCSwitchListAnalyzed1$isoformFeatures,
                lab = NPCSwitchListAnalyzed1$isoformFeatures$gene_name,
                x = 'dIF',
                y = 'isoform_switch_q_value',
                title = 'NPCs, ASD vs CTL',
                xlim= c(-1,1),
                ylim = c(-0.5,4),
                xlab="dIF",
                ylab="-log10(isoform_switch_q_value)", 
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = T,
                legendLabels =c("NS", "dIF", "q-value", "dIF and q-value")
)

#create a copy to modify
NPCSwitchListAnalyzed2 <-NPCSwitchListAnalyzed1

#keep only the isoforms with padj <0.05
selected_isoforms <- NPCSwitchListAnalyzed1$isoformFeatures$isoform_switch_q_value <0.05
NPCSwitchListAnalyzed2$isoformFeatures <- NPCSwitchListAnalyzed1$isoformFeatures[selected_isoforms,]
NPCSwitchListAnalyzed2$orfAnalysis <- NPCSwitchListAnalyzed1$orfAnalysis[selected_isoforms,]
NPCSwitchListAnalyzed2$isoformSwitchAnalysis <- NPCSwitchListAnalyzed1$isoformSwitchAnalysis[selected_isoforms,]
NPCSwitchListAnalyzed2$AlternativeSplicingAnalysis <- NPCSwitchListAnalyzed1$AlternativeSplicingAnalysis[selected_isoforms,]

#order by abs(dIF)
NPCSwitchListAnalyzed2$isoformSwitchAnalysis <- NPCSwitchListAnalyzed2$isoformSwitchAnalysis[order(-abs(NPCSwitchListAnalyzed2$isoformFeatures$dIF)),]
NPCSwitchListAnalyzed2$AlternativeSplicingAnalysis <- NPCSwitchListAnalyzed2$AlternativeSplicingAnalysis[order(-abs(NPCSwitchListAnalyzed2$isoformFeatures$dIF)),]
NPCSwitchListAnalyzed2$isoformFeatures <- NPCSwitchListAnalyzed2$isoformFeatures[order(-abs(NPCSwitchListAnalyzed2$isoformFeatures$dIF)),]


NPC_griesi_DTU <-NPCSwitchListAnalyzed2$isoformFeatures[, c(3,4,8,9,27,28,32)]
NPC_griesi_DTU <-NPCSwitchListAnalyzed2$isoformFeatures
list_NPC_griesi_DTU <-NPC_griesi_DTU$gene_id

#overlap with IFN DE
overlap_DTU_griesi_DE_IFN_NPC <-NULL
for (gene in list_NPC_griesi_DTU) {
  if (gene %in% df_DE_list_NPC_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_NPC<-c(overlap_DTU_griesi_DE_IFN_NPC, gene) 
}

overlap_DTU_griesi_DE_IFN_NPC #NULL

overlap_DTU_griesi_DE_IFN_late <-NULL
for (gene in list_NPC_griesi_DTU) {
  if (gene %in% df_DE_list_late_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_late <-c(overlap_DTU_griesi_DE_IFN_late, gene) 
}

overlap_DTU_griesi_DE_IFN_late #NULL

overlap_DTU_griesi_DE_IFN_early <-NULL
for (gene in list_NPC_griesi_DTU) {
  if (gene %in% df_DE_list_early_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_early <-c(overlap_DTU_griesi_DE_IFN_early, gene) 
}
overlap_DTU_griesi_DE_IFN_early

overlap_DTU_griesi_DE_IFN_rep <-NULL
for (gene in list_NPC_griesi_DTU) {
  if (gene %in% df_DE_list_repvslate_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_rep <-c(overlap_DTU_griesi_DE_IFN_rep, gene) 
}

overlap_DTU_griesi_DE_IFN_rep #NULL

Sfari_genes <- read.csv("~/Desktop/MRes/20220611/SFARI.csv", header =TRUE, stringsAsFactors = FALSE)
Sfari_high_score <- subset(Sfari_genes, gene.score %in% c(1,2))

list_SFARI <- Sfari_high_score$gene.symbol

go.obj <-newGeneOverlap(list_SFARI,list_NPC_griesi_DTU , 20068 )
go.obj <-testGeneOverlap(go.obj)
go.obj

View(NPCSwitchListAnalyzed2$isoformFeatures)
View(NPCSwitchListAnalyzed2$isoformSwitchAnalysis)
View(NPCSwitchListAnalyzed2$AlternativeSplicingAnalysis)
write.csv(NPCSwitchListAnalyzed2$isoformFeatures, file= "GriesiNPCSwitchListAnalyzed2$isoformFeatures_IFN20220704.csv", row.names = T)
write.csv(NPCSwitchListAnalyzed2$isoformSwitchAnalysis, file= "GriesiNPCSwitchListAnalyzed2$isoformSwitchAnalysis_IFN20220704.csv", row.names = T)
write.csv(NPCSwitchListAnalyzed2$AlternativeSplicingAnalysis, file= "GRiesiNPCSwitchListAnalyzed2$AlternativeSplicingAnalysis_IFN20220704.csv", row.names = T)
writeLines(NPCSwitchListAnalyzed2$isoformFeatures$gene_name, "DTU_griesi_NPC_genessignif_20220704.txt")

griesi_genes_NPC <-NPCSwitchListAnalyzed2$AlternativeSplicingAnalysis[, c(1,2,5,8,11,14,17,20,23)]
griesi_genes_NPC$gene_id <-NPC_griesi_DTU$gene_id
griesi_genes_NPC<-griesi_genes_NPC[, c(1,10,2,3,4,5,6,7,8,9)]
griesi_genes_NPC <-griesi_genes_NPC[order(griesi_genes_NPC$gene_id),]

NPC_griesi_DTU <- NPC_griesi_DTU[order(NPC_griesi_DTU$gene_id),]
griesi_genes_NPC$iso_biotype <-NPC_griesi_DTU$iso_biotype
griesi_genes_NPC$dIF <-NPC_griesi_DTU$dIF
griesi_genes_NPC$isoform_switch_q_value <-NPC_griesi_DTU$isoform_switch_q_value
griesi_genes_NPC$switchConsequencesGene <-NPC_griesi_DTU$switchConsequencesGene
griesi_genes_NPC<- griesi_genes_NPC[, c(1, 2, 11, 12, 13, 3, 4, 5, 6, 7, 8, 9, 10)]


# visualise consequences
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

NPCSwitchListAnalyzed1 <- analyzeSwitchConsequences(
  NPCSwitchListAnalyzed1,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0,
  showProgress=TRUE
)

extractTopSwitches(
  NPCSwitchListAnalyzed2, 
  filterForConsequences = FALSE, 
  n = 10, 
  sortByQvals = TRUE
)

switchingIso <- extractTopSwitches( 
  NPCSwitchListAnalyzed1,
  filterForConsequences = FALSE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

#switchplots
switchPlot(NPCSwitchListAnalyzed1, gene = 'TMEM54')
switchPlot(NPCSwitchListAnalyzed1, gene = 'RPS4X')

extractSplicingSummary(NPCSwitchListAnalyzed2)

splicing_enrichment <-extractSplicingEnrichment(NPCSwitchListAnalyzed1, 
                                                minEventsForPlotting =1)
splicing_enrichment
#No features left for plotting after filtering with via "minEventsForPlotting" argument.

extractSplicingGenomeWide(NPCSwitchListAnalyzed1)

######
#neurons
design_neurons_griesi <-data.frame(
  sampleID = colnames(kallistoQuant_griesi$abundance)[31:49],
  condition = c(rep("ASD",9), rep("CTL",10))
)

design_neurons_griesi$condition <-factor(design_neurons_griesi$condition, levels = c("CTL", "ASD"))

neurons_griesi_SwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant_griesi$counts,
  isoformRepExpression = kallistoQuant_griesi$abundance,
  designMatrix         = design_neurons_griesi,
  isoformExonAnnoation = "~/Desktop/MRes/20210402/gencode.v34.annotation.gtf",
  isoformNtFasta       = "~/Desktop/MRes/20210402/gencode.v34.transcripts_norRNA.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

switchList_neurons <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = neurons_griesi_SwitchList,
  pathToOutput = "~/Desktop/MRes/20220625/neurons_griesi/",
  outputSequences      = TRUE,
  prepareForWebServers = TRUE
)

View(switchList_neurons$isoformFeatures)

neuronsSwitchListAnalyzed1 <- analyzeCPC2(
  switchAnalyzeRlist   = switchList_neurons,
  pathToCPC2resultFile = "~/Desktop/MRes/20220625/neurons_griesi/result_cpc2.txt",
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

neuronsSwitchListAnalyzed1  <- analyzePFAM(
  switchAnalyzeRlist   = neuronsSwitchListAnalyzed1 ,
  pathToPFAMresultFile = "~/Desktop/MRes/20220625/neurons_griesi/Pfam_neurons_griesi.txt",
  showProgress=TRUE
)

#neuronsSwitchListAnalyzed1  <- analyzeSignalP(switchAnalyzeRlist       = neuronsSwitchListAnalyzed1 ,pathToSignalPresultFile  = "~/Desktop/MRes/20220625/neurons_griesi/signalP_neurons_griesi.txt")
#  No signial peptides were found

#neuronsSwitchListAnalyzed1 <- analyzeIUPred2A(switchAnalyzeRlist        = neuronsSwitchListAnalyzed1,pathToIUPred2AresultFile = "~/Desktop/MRes/20220625/neurons_griesi/IUPred2A_neurons_griesi.txt")
#The 'pathToIUPred2AresultFile' does not appear to contain results for the isoforms stored in the switchAnalyzeRlist...


#assess alternative splicing
neuronsSwitchListAnalyzed1 <- analyzeIntronRetention(neuronsSwitchListAnalyzed1 )
neuronsSwitchListAnalyzed1 <-analyzeAlternativeSplicing(neuronsSwitchListAnalyzed1)

#consequences
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

neuronsSwitchListAnalyzed1 <- analyzeSwitchConsequences(
  neuronsSwitchListAnalyzed1,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0,
  showProgress=TRUE
)

extractTopSwitches(
  neuronsSwitchListAnalyzed1, 
  filterForConsequences = FALSE, 
  n = 10, 
  sortByQvals = TRUE
)

switchingIso <- extractTopSwitches( 
  neuronsSwitchListAnalyzed1, 
  filterForConsequences = FALSE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

#switchplots
switchPlot(neuronsSwitchListAnalyzed1, gene = 'HLA-A')
switchPlot(neuronsSwitchListAnalyzed1, gene = 'TMEM119')

#volcano plot
EnhancedVolcano(neuronsSwitchListAnalyzed1$isoformFeatures,
                lab = neuronsSwitchListAnalyzed1$isoformFeatures$gene_name,
                x = 'dIF',
                y = 'isoform_switch_q_value',
                title = 'neuronal cells, ASD vs CTL',
                xlim= c(-0.5,0.5),
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

#copy to modify
neuronsSwitchListAnalyzed2 <-neuronsSwitchListAnalyzed1

#keep only the isoforms with padj <0.05
selected_isoforms <- neuronsSwitchListAnalyzed1$isoformFeatures$isoform_switch_q_value <0.05
neuronsSwitchListAnalyzed2$isoformFeatures <- neuronsSwitchListAnalyzed1$isoformFeatures[selected_isoforms,]
neuronsSwitchListAnalyzed2$orfAnalysis <- neuronsSwitchListAnalyzed1$orfAnalysis[selected_isoforms,]
neuronsSwitchListAnalyzed2$isoformSwitchAnalysis <- neuronsSwitchListAnalyzed1$isoformSwitchAnalysis[selected_isoforms,]
neuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis <- neuronsSwitchListAnalyzed1$AlternativeSplicingAnalysis[selected_isoforms,]

#order by dIF
neuronsSwitchListAnalyzed2$isoformSwitchAnalysis <- neuronsSwitchListAnalyzed2$isoformSwitchAnalysis[order(-abs(neuronsSwitchListAnalyzed2$isoformFeatures$dIF)),]
neuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis <- neuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis[order(-abs(neuronsSwitchListAnalyzed2$isoformFeatures$dIF)),]
neuronsSwitchListAnalyzed2$isoformFeatures <- neuronsSwitchListAnalyzed2$isoformFeatures[order(-abs(neuronsSwitchListAnalyzed2$isoformFeatures$dIF)),]

View(neuronsSwitchListAnalyzed2$isoformFeatures)
View(neuronsSwitchListAnalyzed2$isoformSwitchAnalysis)
View(neuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis)
write.csv(neuronsSwitchListAnalyzed2$isoformFeatures, file= "GriesineuronsSwitchListAnalyzed2$isoformFeatures_IFN20220704.csv", row.names = T)
write.csv(neuronsSwitchListAnalyzed2$isoformSwitchAnalysis, file= "GriesineuronsSwitchListAnalyzed2$isoformSwitchAnalysis_IFN20220704.csv", row.names = T)
write.csv(neuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis, file= "GRiesineuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis_IFN20220704.csv", row.names = T)
writeLines(neuronsSwitchListAnalyzed2$isoformFeatures$gene_name, "DTU_griesi_neurons_genessignif_20220704.txt")

neurons_griesi_DTU <-neuronsSwitchListAnalyzed2$isoformFeatures[, c(3,4,8,9,27,28,32)]
View(neurons_griesi_DTU)
list_neurons_griesi_DTU <-neurons_griesi_DTU$gene_id

#overlap with IFN DE
overlap_DTU_griesi_DE_IFN_NPC <-NULL
for (gene in list_neurons_griesi_DTU) {
  if (gene %in% df_DE_list_NPC_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_NPC<-c(overlap_DTU_griesi_DE_IFN_NPC, gene) 
}
overlap_DTU_griesi_DE_IFN_NPC #HLA-A"

overlap_DTU_griesi_DE_IFN_late <-NULL
for (gene in list_neurons_griesi_DTU) {
  if (gene %in% df_DE_list_late_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_late <-c(overlap_DTU_griesi_DE_IFN_late, gene)
}

overlap_DTU_griesi_DE_IFN_late #"HLA-A"

overlap_DTU_griesi_DE_IFN_early <-NULL
for (gene in list_neurons_griesi_DTU) {
  if (gene %in% df_DE_list_early_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_early <-c(overlap_DTU_griesi_DE_IFN_early, gene) 
}
overlap_DTU_griesi_DE_IFN_early #HLA-A

overlap_DTU_griesi_DE_IFN_rep <-NULL
for (gene in list_neurons_griesi_DTU) {
  if (gene %in% df_DE_list_repvslate_IFN$SYMBOL)
    overlap_DTU_griesi_DE_IFN_rep <-c(overlap_DTU_griesi_DE_IFN_rep, gene) 
}

overlap_DTU_griesi_DE_IFN_rep #HLA-A

#overlap Sfari
go.obj <-newGeneOverlap(list_SFARI,list_neurons_griesi_DTU, 20068 )
go.obj <-testGeneOverlap(go.obj)
go.obj
getIntersection((go.obj))

griesi_genes_neurons <-neuronsSwitchListAnalyzed2$AlternativeSplicingAnalysis[, c(1,2,5,8,11,14,17,20,23)]
griesi_genes_neurons$gene_id <-neurons_griesi_DTU$gene_id
griesi_genes_neurons<-griesi_genes_neurons[, c(1,10,2,3,4,5,6,7,8,9)]
View(griesi_genes_neurons)
griesi_genes_neurons <-griesi_genes_neurons[order(griesi_genes_neurons$gene_id),]
neurons_griesi_DTU <- neurons_griesi_DTU[order(neurons_griesi_DTU$gene_id),]
griesi_genes_neurons$iso_biotype <-neurons_griesi_DTU$iso_biotype
griesi_genes_neurons$dIF <-neurons_griesi_DTU$dIF
griesi_genes_neurons$isoform_switch_q_value <-neurons_griesi_DTU$isoform_switch_q_value
griesi_genes_neurons$switchConsequencesGene <-neurons_griesi_DTU$switchConsequencesGene
griesi_genes_neurons <- griesi_genes_neurons[, c(1, 2, 11, 12, 13, 3, 4, 5, 6, 7, 8, 9, 10)]


extractSplicingSummary(neuronsSwitchListAnalyzed2)

splicing_enrichment_neurons <-extractSplicingEnrichment( neuronsSwitchListAnalyzed2, 
                                                         minEventsForPlotting =1)
splicing_enrichment
#No features left for ploting after filtering with via "minEventsForPlotting" argument.
extractSplicingEnrichmentComparison(neuronsSwitchListAnalyzed2)
#Cannot do a contrast of different comparisons since only one comparison is analyzed in the switchAnalyzeRlist.
extractSplicingGenomeWide( neuronsSwitchListAnalyzed2)
