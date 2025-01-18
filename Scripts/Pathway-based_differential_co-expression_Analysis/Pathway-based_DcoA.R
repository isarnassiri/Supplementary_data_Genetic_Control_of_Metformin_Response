
############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
library(data.table)
library(qusage)
library(clusterProfiler)
library(simplifyEnrichment)
library(cowplot)
library(DGCA, quietly = TRUE) 
library(ggpubr)
library(dplyr)
library(biomaRt)

#----------------- Description -----------------
# A coExQTL's functional impact would be vary depending on the specific context and the biological significance of the gene within the pathway. 
# We used a pathway-based differential co-expression analysis to determine the permissible range of expression variation within a pathway that would classify a coExQTL as functionally disruptive. 
# Pathway-based differential co-expression analysis was performed on a subset of genes that have allele-specific co-expression relationships and enriched with curated gene sets from online pathway databases (FDR < 0.01). 
# The analysis uncovers the average shift in correlation between gene expression among two genotype classes, as well as its statistical significance.

#----------------- Output -----------------
# Results of pathway enrichment analysis and Pathway-based differential co-expression analysis as text files.

#----------------- Examples -----------------
treatments = c('IFN')
selected_genes = c('ERAP2')
selected_SNPs = c('rs2910789')

# Execute the rest of the code line by line.

##########################################################################################################################################################


#----- Step-01: pathway enrichment analysis

setwd('~')
IFN = fread('All_RELICATED_array_RNAseq_IFN.txt', stringsAsFactors = F, header = T)
LPS24 = fread('All_RELICATED_array_RNAseq_LPS24.txt', stringsAsFactors = F, header = T)
UT = fread('All_RELICATED_array_RNAseq_UT.txt', stringsAsFactors = F, header = T)

IFN$state = rep('IFN', dim(IFN)[1])
LPS24$state = rep('LPS24', dim(LPS24)[1])
UT$state = rep('UT', dim(UT)[1])

list = list(ifn=IFN, lps=LPS24, ut=UT)

merged = do.call(rbind, list)
merged = as.data.frame(merged)
merged = merged[!duplicated(merged$ind_RNAseq2), ]

#----- read in the reference gene sets
gmtfile <- c("~/hc2.all.v2024.1.Hs.symbols.gmt")  #h.all.v2024.1.Hs.symbols.gmt
c5 <- read.gmt(gmtfile)
#-----

#--------------------------- query ---------------------------
# treatments = c('IFN')
# selected_genes = c('ERAP2')
# selected_SNPs = c('rs2910789')
#---------------------------

#-- enrichment of all allele-specific co-expression relationships 
for (treatment in treatments){
  
  selected_gene = selected_genes[grep(treatment, treatments)]
  print(treatment)
  print(selected_gene)
  
  #--
  selected = merged[which( merged$state == treatment ), ] # & merged$replication == 'Yes' ; merged[which(merged$Gene2_RNAseq == selected_gene & merged$state == treatment ), ]
  print(dim(selected))

  fwrite(selected, paste0(treatment, '_RNAseq_Array.txt'), quote = F, row.names = F, sep = '\t')
  
  unique(selected$var_id_array)
  summary(abs(selected$zScoreDiff_RNAseq))
  
  egmt <- enricher(unique(c(selected$Gene1_RNAseq, selected$Gene2_RNAseq)), TERM2GENE=c5)
  fwrite(egmt@result, paste(treatment, '_AllcoExQTL_PathwayEnrichment_c2_HALLMARK.txt', sep = '_'), quote = F, row.names = F, sep = '\t')
 
}

#----- Step-02: pathway score analysis
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
head(listAttributes(ensembl))
filters = listFilters(ensembl)
GTF <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(GTF) = c('symbol', 'id')

for (treatment in treatments) {

  Gene2 = selected_genes[grep(treatment, treatments)]
  SNP = selected_SNPs[grep(treatment, treatments)]
  
  print(treatment)
  print(Gene2)
  print(SNP)

  Enriched_pathways = fread(paste0('~/', treatment, '__AllcoExQTL_PathwayEnrichment_c2_HALLMARK.txt'), stringsAsFactors = F)
  Enriched_pathways = as.data.frame(Enriched_pathways)
  Enriched_pathways = Enriched_pathways[order(Enriched_pathways$p.adjust, decreasing = F),]
  Enriched_pathways = Enriched_pathways[which(Enriched_pathways$p.adjust < 1e-2),]

  print(dim(Enriched_pathways))
 
  pathway_genes = strsplit(Enriched_pathways$geneID, "/")
  names(pathway_genes) = Enriched_pathways$ID
  
  pathways = stack(pathway_genes) 
  pathways$ind = as.character(pathways$ind) 
  str(pathways)
  
  #-------------------
  expression = fread(paste0('~/expression_',treatment,'.txt'), stringsAsFactors = F)
  expression = as.data.frame(expression)
  row.names(expression) = expression$id
  expression = expression[,-1]
  dim(expression)
  
  PROBEID = data.frame(id = row.names(expression), stringsAsFactors = F)
  SYMBOL = right_join(GTF, PROBEID, by ='id')
  SYMBOL = SYMBOL[-which(SYMBOL$symbol ==''),]
  SYMBOL = SYMBOL[!is.na(SYMBOL$symbol),]
  SYMBOL = SYMBOL[!duplicated(SYMBOL$symbol),]
  SYMBOL = SYMBOL[!duplicated(SYMBOL$id),]
  
  expression = expression[which(row.names(expression) %in% SYMBOL$id),]
  SYMBOL = SYMBOL[which(SYMBOL$id %in% row.names(expression)),]
  
  expression = expression[!duplicated(row.names(expression)),]
  dim(expression)
  length(SYMBOL$symbol)
  
  SYMBOL = SYMBOL[match(row.names(expression), SYMBOL$id),]
  identical(row.names(expression), SYMBOL$id)
  
  row.names(expression) = SYMBOL$symbol
  
  TP_GENOTYPE = fread(paste0('~/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2_USED_for_BOXPLOT.txt'),skip = paste0(SNP, '_'), nrows = 1, stringsAsFactors = F, header = F)
  TP_GENOTYPE$V1 = gsub('_', '', TP_GENOTYPE$V1)
  
  GENOTYPE_sample_names = fread(paste0('~/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),nrows = 1, stringsAsFactors = F)
  colnames(TP_GENOTYPE) = colnames(GENOTYPE_sample_names)
  
  expression=as.data.frame(expression)
  TP_GENOTYPE=as.data.frame(TP_GENOTYPE)
  expression = expression[,which(colnames(expression) %in% colnames(TP_GENOTYPE))]
  TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression))]
  
  dim(TP_GENOTYPE)
  dim(expression) 
  
  #---------- Ref homozygote
  expr0=expression[, which(colnames(expression) %in% colnames(TP_GENOTYPE)[which(TP_GENOTYPE[1,]==0)])]
  expr1=expression[, which(colnames(expression) %in% colnames(TP_GENOTYPE)[which(TP_GENOTYPE[1,]==1)])]
  expr2=expression[, which(colnames(expression) %in% colnames(TP_GENOTYPE)[which(TP_GENOTYPE[1,]==2)])]
  #-------------------
  
  set.seed(1234)
  design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr0)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr0)[2])))
  colnames(design_matrix) <- c("G2", "G0")
  
  #-- G2_G0 
  INPUT = cbind(expr2, expr0)
  
  moduleDC_res_G2_G0 = moduleDC(inputMat = INPUT, design = design_matrix,
                                compare = c("G2", "G0"), genes = pathways$values,
                                labels = pathways$ind, nPerm = 1000, number_DC_genes = 10,
                                dCorAvgMethod = "median", corrType = "spearman")
  
  head(moduleDC_res_G2_G0)
  
  design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr1)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr1)[2])))
  colnames(design_matrix) <- c("G2", "G1")
  
  #-- G2_G1
  INPUT = cbind(expr2, expr1)
  
  moduleDC_res_G2_G1 = moduleDC(inputMat = INPUT, design = design_matrix,
                                compare = c("G2", "G1"), genes = pathways$values,
                                labels = pathways$ind, nPerm = 1000, number_DC_genes = 10,
                                dCorAvgMethod = "median", corrType = "spearman")
  
  head(moduleDC_res_G2_G1)
  
  fwrite(moduleDC_res_G2_G0, paste0('~/ModuleDC_', treatment,'_G2_G0.txt'), quote = F, row.names = F, sep = '\t')
  fwrite(moduleDC_res_G2_G1, paste0('~/ModuleDC_', treatment,'_G2_G1.txt'), quote = F, row.names = F, sep = '\t')
  #-------------------
}
