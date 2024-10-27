######################## First do pathway enrichment analysis

#-- intersection of RNAseq and array
library(data.table)
treat = c('IFN', 'LPS24', 'UT')
i=1

setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/gQTL_coExQTL/')

for(i in 1:length(treat))
{
  
  RNAseq = fread(paste0(treat[i], '_coExQTL_gQTL.txt'), stringsAsFactors = F, header = T)
  array = fread(paste0(treat[i], '_coExQTL_array_gQTL.txt'), stringsAsFactors = F, header = T)
  
  array$Gene1 = gsub('\\..*' , '', array$Gene1)
  array$Gene2 = gsub('\\..*' , '', array$Gene2)
  
  RNAseq$index = paste(RNAseq$Gene1, RNAseq$Gene2, RNAseq$var_id, sep='_')
  array$index = paste(array$Gene1, array$Gene2, array$var_id, sep='_')
  array$index2 = paste(array$Gene2, array$Gene1, array$var_id, sep='_')
  
  RNAseq = RNAseq[order(RNAseq$pValDiff_adj, decreasing = F), ]
  RNAseq = RNAseq[RNAseq$pValDiff_adj <= 5e-2, ]
  array = array[array$pValDiff_adj <= 5e-2, ] 
  
  dim(RNAseq)
  dim(array)
  
  length(unique(RNAseq$Gene1))
  length(unique(RNAseq$Gene2))
  length(unique(RNAseq$var_id))
  range(RNAseq$pValDiff_adj)
  
  length(unique(array$Gene1))
  length(unique(array$Gene2))
  length(unique(array$var_id))
  range(array$pValDiff_adj)
  
  MERGED = merge(RNAseq, array, by = 'index')
  MERGED = MERGED[!duplicated(MERGED$index), ]
  
  length(unique(MERGED$Gene1.x))
  length(unique(MERGED$Gene2.x))
  length(unique(MERGED$var_id.x))
  range(MERGED$pValDiff_adj.x)
  
  print(treat[i])
  # print(length(unique(RNAseq$Gene1)))
  # print(length(unique(RNAseq$var_id)))
  print(table((MERGED$Classes.x == MERGED$Classes.y)))
  df = data.frame(table((MERGED$Classes.x == MERGED$Classes.y)))
  print(dim(MERGED))
  print(df[2, 2]/dim(MERGED)[1])

  #---
  head(MERGED)
  MERGED$replication = NA
  MERGED$replication[which(MERGED$Classes.x == MERGED$Classes.y)] = 'Yes'
  MERGED$replication[-which(MERGED$Classes.x == MERGED$Classes.y)] = 'No'
  table(MERGED$replication)
  colnames(MERGED) = gsub('.x', '_RNAseq', colnames(MERGED))
  colnames(MERGED) = gsub('.y', '_array', colnames(MERGED))
  
  table(MERGED$replication)
  
  #MERGED = MERGED[, -grep('ind', colnames(MERGED))]
  
  fwrite(MERGED, paste0('All_RELICATED_array_RNAseq_', treat[i], '.txt'), quote = F, row.names = F, sep = '\t')
}

#--------- comparison and make summary
library(data.table)
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

#=========================================================== enricher + generating supplementary file 7 ===============================
#----- clusterProfiler
library(qusage)
gmtfile <- c("/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/hc2.all.v2024.1.Hs.symbols.gmt")  #h.all.v2024.1.Hs.symbols.gmt
c5 <- read.gmt(gmtfile)
#-----

library(data.table)
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/')

treatments = c('IFN', 'UT', 'LPS24')
selected_genes = c('ERAP2', 'RPS26', 'OXR1')

# treatment = 'IFN'
# selected_gene = 'ERAP2' # "rs2910789", 1576
# 
# treatment = 'UT'
# selected_gene = 'RPS26' # "rs7305461", 496
# 
# treatment = 'LPS24'
# selected_gene = 'OXR1' # "rs3110426", 1057   

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
  
  library(clusterProfiler)
  egmt <- enricher(unique(c(selected$Gene1_RNAseq, selected$Gene2_RNAseq)), TERM2GENE=c5)
  fwrite(egmt@result, paste(treatment, '_AllcoExQTL_PathwayEnrichment_c2_HALLMARK.txt', sep = '_'), quote = F, row.names = F, sep = '\t')
 
}

#=========================================================== pathway score
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
head(listAttributes(ensembl))
filters = listFilters(ensembl)
GTF <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(GTF) = c('symbol', 'id')

library(data.table)

treatments = c('IFN', 'UT', 'LPS24')
selected_genes = c('ERAP2', 'RPS26', 'OXR1')
selected_SNPs = c('rs2910789', 'rs7305461', 'rs3110426')

# treatment = c('UT' )
# Gene2 = c('RPS26' )
# SNP = c('rs7305461' )

# treatment = c('LPS24' )
# Gene2 = c('OXR1' )
# SNP = c('rs3110426' ) 

treatment = c('IFN' )
Gene2 = c('ERAP2' )
SNP = c('rs2910789' ) 

for (treatment in treatments) {

  Gene2 = selected_genes[grep(treatment, treatments)]
  SNP = selected_SNPs[grep(treatment, treatments)]
  
  print(treatment)
  print(Gene2)
  print(SNP)

  Enriched_pathways = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/', treatment, '__AllcoExQTL_PathwayEnrichment_c2_HALLMARK.txt'), stringsAsFactors = F)
  Enriched_pathways = as.data.frame(Enriched_pathways)
  Enriched_pathways = Enriched_pathways[order(Enriched_pathways$p.adjust, decreasing = F),]
  Enriched_pathways = Enriched_pathways[which(Enriched_pathways$p.adjust < 1e-2),]
  
  if(treatment == 'IFN')
  {
    # Enriched_pathways = Enriched_pathways[which(Enriched_pathways$p.adjust < 1e-10),]
    
    # Enriched_pathways = Enriched_pathways[grep('HALLMARK', Enriched_pathways$ID),]
    # 
    # toMatch <- c("MYC", "INTERFERON", "OXIDATIVE", "MITOCHOND", "RAS")
    # matches <- unique(grep(paste(toMatch,collapse="|"), Enriched_pathways$ID, value=TRUE))
    # 
    # Enriched_pathways = Enriched_pathways[which(Enriched_pathways$ID %in% matches),]
  }
  
  print(dim(Enriched_pathways))
 
  pathway_genes = strsplit(Enriched_pathways$geneID, "/")
  names(pathway_genes) = Enriched_pathways$ID
  
  pathways = stack(pathway_genes) 
  pathways$ind = as.character(pathways$ind) 
  str(pathways)
  
  #-------------------
  setwd('/Users/isarnassiri/Documents/RESULTS_USED')
  
  library(data.table)
  expression = fread(paste0('Expression/Input_files_gene/expression_',treatment,'.txt'), stringsAsFactors = F)
  expression = as.data.frame(expression)
  row.names(expression) = expression$id
  expression = expression[,-1]
  dim(expression)
  
  library(dplyr)
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
  
  TP_GENOTYPE = fread(paste0('Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2_USED_for_BOXPLOT.txt'),skip = paste0(SNP, '_'), nrows = 1, stringsAsFactors = F, header = F)
  TP_GENOTYPE$V1 = gsub('_', '', TP_GENOTYPE$V1)
  
  GENOTYPE_sample_names = fread(paste0('Genotype//Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),nrows = 1, stringsAsFactors = F)
  colnames(TP_GENOTYPE) = colnames(GENOTYPE_sample_names)
  
  expression=as.data.frame(expression)
  TP_GENOTYPE=as.data.frame(TP_GENOTYPE)
  expression = expression[,which(colnames(expression) %in% colnames(TP_GENOTYPE))]
  TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression))]
  
  dim(TP_GENOTYPE)
  dim(expression) 
  
  #---------- Ref homozygote
  library(ggpubr)
  expr0=expression[, which(colnames(expression) %in% colnames(TP_GENOTYPE)[which(TP_GENOTYPE[1,]==0)])]
  expr1=expression[, which(colnames(expression) %in% colnames(TP_GENOTYPE)[which(TP_GENOTYPE[1,]==1)])]
  expr2=expression[, which(colnames(expression) %in% colnames(TP_GENOTYPE)[which(TP_GENOTYPE[1,]==2)])]
  #-------------------
  
  set.seed(1234)
  design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr0)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr0)[2])))
  colnames(design_matrix) <- c("G2", "G0")
  
  #-- G2_G0 
  INPUT = cbind(expr2, expr0)
  
  library(DGCA, quietly = TRUE) 
  moduleDC_res_G2_G0 = moduleDC(inputMat = INPUT, design = design_matrix,
                          compare = c("G2", "G0"), genes = pathways$values,
                          labels = pathways$ind, nPerm = 1000, number_DC_genes = 10,
                          dCorAvgMethod = "median")
  
  head(moduleDC_res_G2_G0)
  
  design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr1)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr1)[2])))
  colnames(design_matrix) <- c("G2", "G1")
  
  #-- G2_G1
  INPUT = cbind(expr2, expr1)
  
  moduleDC_res_G2_G1 = moduleDC(inputMat = INPUT, design = design_matrix,
                          compare = c("G2", "G1"), genes = pathways$values,
                          labels = pathways$ind, nPerm = 1000, number_DC_genes = 10,
                          dCorAvgMethod = "median")
  
  head(moduleDC_res_G2_G1)
  
  fwrite(moduleDC_res_G2_G0, paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/ModuleDC_', treatment,'_G2_G0.txt'), quote = F, row.names = F, sep = '\t')
  fwrite(moduleDC_res_G2_G1, paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/ModuleDC_', treatment,'_G2_G1.txt'), quote = F, row.names = F, sep = '\t')
  #-------------------
  

  # Module-based differential correlation
  # set.seed(1234)
  # design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr0)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr0)[2])))
  # colnames(design_matrix) <- c("G2", "G0")
  # dim(design_matrix)
  # 
  # #-- G2_G0 
  # INPUT = cbind(expr2, expr0)
  # 
  # pathway_genes = pathways[pathways$ind == "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "values"]
  # INPUT_pathway = INPUT[pathway_genes, ]
  # dim(INPUT_pathway)
  # 
  # moduleDC_res = ddcorAll(inputMat = INPUT_pathway, design = design_matrix,
  #                         compare = c("G2", "G0"), nPerm = 50, getDCorAvg = TRUE, dCorAvgType = "gene_average", dCorAvgMethod = "median")
  # 
  # View(moduleDC_res[["avg_dcor"]])
  
}

#=========================================================== enricher + generating supplementary file 7 + visualization ===============================

setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/gQTL_coExQTL/')

#--------- read in profiles
library(data.table)
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

#--
# save merged as a supp

#----- clusterProfiler
library(qusage)
gmtfile <- c("/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/hc2.all.v2024.1.Hs.symbols.gmt")  #h.all.v2024.1.Hs.symbols.gmt
c5 <- read.gmt(gmtfile)
#-----

library(data.table)
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/')

treatments = c('IFN', 'UT', 'LPS24')
selected_genes = c('ERAP2', 'RPS26', 'OXR1')

treatment = 'IFN'

#-- enrichment of all allele-specific co-expression relationships 
for (treatment in treatments){

  selected_gene = selected_genes[grep(treatment, treatments)]
  print(treatment)
  print(selected_gene)
  
  selected = merged[which( merged$state == treatment ), ] # & merged$replication == 'Yes' ; merged[which(merged$Gene2_RNAseq == selected_gene & merged$state == treatment ), ]
  print(dim(selected))

  library(clusterProfiler)
  library(simplifyEnrichment)
  egmt <- enricher(unlist(unique(c(selected$Gene1_RNAseq, selected$Gene2_RNAseq))), TERM2GENE=c5, pvalueCutoff = 0.01)

  #--
  moduleDC_res_G2_G0 = fread( paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/USED/ModuleDC_', treatment,'_G2_G0.txt'), stringsAsFactors = F)
  moduleDC_res_G2_G1 = fread( paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/PathwayScore/USED/ModuleDC_', treatment,'_G2_G1.txt'), stringsAsFactors = F)
  
  moduleDC_res_G2_G0 = moduleDC_res_G2_G0[which(moduleDC_res_G2_G0$pVal < 0.05),]
  moduleDC_res_G2_G1 = moduleDC_res_G2_G1[which(moduleDC_res_G2_G1$pVal < 0.05),]
  
  moduleDC = rbind(moduleDC_res_G2_G0, moduleDC_res_G2_G1)
  moduleDC = moduleDC[order(moduleDC$pVal, decreasing = F),]
  moduleDC = moduleDC[!duplicated(moduleDC$Module),]
  #--
  
  if(treatment == 'UT')
  {
    
    egmt_subset = subset_enrichResult(egmt, moduleDC$Module)
    View(egmt_subset@result)
    
    terms = grep(paste(c("HALLMARK_PROTEIN_SECRETION", "REACTOME_SARS_COV_1_INFECTION", "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"), collapse="|"), egmt_subset@result$ID[egmt_subset@result$p.adjust < 0.01], ignore.case=TRUE, value=TRUE)
    terms

    egmt_subset = subset_enrichResult(egmt_subset, terms)
    
    egmt_subset@result$Description = gsub('_',' ', gsub('^.*?_', '', egmt_subset@result$Description))
    
    library(cowplot)
    setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/FIGURES_NEW/')
    pdf(paste(treatment, selected_gene, 'cnetplot_gene.pdf', sep = '_'), width = 10, height = 8, useDingbats = F)
    print(cnetplot(egmt_subset, showCategory = 4, circular = TRUE, colorEdge = TRUE, categorySize="pvalue", color_category='firebrick', color_gene='steelblue'))
    dev.off()
    
    # library(cowplot)
    # setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/FIGURES_NEW/')
    # pdf(paste(treatment, selected_gene, 'cnetplot_category.pdf', sep = '_'), width = 10, height = 10, useDingbats = F)
    # cnetplot(egmt_subset, node_label = 'category', showCategory = 4, circular = TRUE, colorEdge = TRUE, categorySize="pvalue", color_category='firebrick', color_gene='steelblue')
    # dev.off()
    
  }
  
  if(treatment == 'LPS24')
  {
    
    egmt_subset = subset_enrichResult(egmt, moduleDC$Module)
    View(egmt_subset@result)

    terms = grep(paste(c("HALLMARK_TNFA_SIGNALING_VIA_NFKB" , "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_HYPOXIA"), collapse="|"), egmt_subset@result$ID[egmt_subset@result$p.adjust < 0.01], ignore.case=TRUE, value=TRUE)
    terms
    
    egmt_subset = subset_enrichResult(egmt_subset, terms)
    
    egmt_subset@result$Description = gsub('_',' ', gsub('^.*?_', '', egmt_subset@result$Description))
    
    # pathway_genes = strsplit(egmt_subset@result$geneID, "/")
    # selected[selected$Gene1_RNAseq %in% unlist(pathway_genes, use.names=FALSE),]
    
    library(cowplot)
    setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/FIGURES_NEW/')
    pdf(paste(treatment, selected_gene, 'cnetplot_gene.pdf', sep = '_'), width = 10, height = 8, useDingbats = F)
    print(cnetplot(egmt_subset, node_label = 'gene', showCategory = 4, circular = F, colorEdge = TRUE, categorySize="pvalue", color_category='firebrick', color_gene='steelblue'))
    dev.off()
    
    pdf(paste(treatment, selected_gene, 'cnetplot_category.pdf', sep = '_'), width = 10, height = 10, useDingbats = F)
    print(cnetplot(egmt_subset, node_label = 'category', showCategory = 4, circular = F, colorEdge = TRUE, categorySize="pvalue", color_category='firebrick', color_gene='steelblue'))
    dev.off()
    
  }

  if(treatment == 'IFN')
  {
    
    egmt_subset = egmt #subset_enrichResult(egmt, moduleDC$Module)
    View(egmt_subset@result)
    
    # These three have sig p-value for pathway score
    # c("BROWNE_HCMV_INFECTION_14HR_DN", "KAUFFMANN_MELANOMA_RELAPSE_UP", "REACTOME_INTERLEUKIN_6_SIGNALING")
    
    terms = grep(paste(c("BROWNE_HCMV_INFECTION_14HR_DN", "KAUFFMANN_MELANOMA_RELAPSE_UP", "REACTOME_INTERLEUKIN_6_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "REACTOME_TRANSCRIPTIONAL_ACTIVATION_OF_MITOCHONDRIAL_BIOGENESIS"), collapse="|"), egmt_subset@result$ID[egmt_subset@result$p.adjust < 0.01], ignore.case=TRUE, value=T)
    terms
    
    egmt_subset = subset_enrichResult(egmt, terms)
    egmt_subset@result$Description = gsub('^.*?_', '', egmt_subset@result$Description)

    library(cowplot)
    setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/FIGURES_NEW/')
    pdf(paste(treatment, selected_gene, 'cnetplot_gene.pdf', sep = '_'), width = 11, height = 8, useDingbats = F)
    print(cnetplot(egmt_subset, node_label = 'gene', showCategory = 5, circular = F, colorEdge = TRUE, categorySize="pvalue", color_category='firebrick', color_gene='steelblue'))
    dev.off()
    
    pdf(paste(treatment, selected_gene, 'cnetplot_category.pdf', sep = '_'), width = 10, height = 10, useDingbats = F)
    print(cnetplot(egmt_subset, node_label = 'category', showCategory = 5, circular = F, colorEdge = TRUE, categorySize="pvalue", color_category='firebrick', color_gene='steelblue'))
    dev.off()
    
  }
  
}
