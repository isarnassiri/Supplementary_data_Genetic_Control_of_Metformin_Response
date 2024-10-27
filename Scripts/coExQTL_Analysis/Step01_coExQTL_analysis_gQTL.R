# module purge
# module add R/4.1.0-foss-2021a
# module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

.libPaths("/well/parkkinen/users/gbf362/R/4.1/skylake/")

treatment="UT"

setwd('/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/')
dir.create('coExQTL/gQTL_coExQTL/', recursive = T)

#--- gene symbol to Ensemble transcript IDs
library(rtracklayer)
gtf_gene <- rtracklayer::import('GENOMEANNOT/GRCh38.gtf')
gtf_gene <- as.data.frame(gtf_gene)

#-- for gQTL
GTF <- gtf_gene[which( gtf_gene$type == 'gene'),]
GTF = GTF[!is.na(gtf_gene$gene_id),]
GTF=GTF[,c("gene_id", "gene_name")]
#--

GTF=as.data.frame(GTF)
colnames(GTF) = c('id','symbol')

library(data.table)
library(dplyr)

expression = fread(paste0('Input_files_gene/expression_',treatment,'.txt'), stringsAsFactors = F)
expression = as.data.frame(expression)
row.names(expression) = expression$id
expression = expression[,-1] 

expression = expression[which(row.names(expression) %in% GTF$id),]
GTF_subset = GTF[which(GTF$id %in% row.names(expression)),]

dim(GTF_subset)
dim(expression)

GTF_subset = GTF_subset[match(row.names(expression), GTF_subset$id),]
identical(GTF_subset$id, row.names(expression))

row.names(expression) = make.names(GTF_subset$symbol,unique=T)

# if you want to read an SNP you should make the file tab separated
TP_GENOTYPE = fread(paste0('Input_files_Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),stringsAsFactors = F)
TP_GENOTYPE = as.data.frame(TP_GENOTYPE)
TP_GENOTYPE = TP_GENOTYPE[!duplicated(TP_GENOTYPE$id),]
row.names(TP_GENOTYPE) = TP_GENOTYPE$id
TP_GENOTYPE = TP_GENOTYPE[,-1] 

expression = expression[,which(colnames(expression) %in% colnames(TP_GENOTYPE))]
TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression))]

dim(expression)
dim(TP_GENOTYPE)

MAF_imputed_Allsamples = fread('MAF_imputed_Allsamples_TYPED.txt', stringsAsFactors = F, header = T)
MAF_imputed_Allsamples = as.data.frame(MAF_imputed_Allsamples)
colnames(MAF_imputed_Allsamples)[which(colnames(MAF_imputed_Allsamples) == 'ID')] = 'var_id'

headerC = fread('gQTL_conditional_Header.txt', stringsAsFactors = F, header = F)
cis = fread(paste0('gQTL/QTLtools_outputs/conditional_pass/',treatment,'/gQTL_conditional_pass_1.txt'), stringsAsFactors = F)
cis = as.data.frame(cis)

colnames(cis) = as.character(headerC)
cis = cis[-which(cis$var_id == '.'),]

dim(cis)
cis = merge(cis, MAF_imputed_Allsamples, by = 'var_id')
dim(cis)

# bwd_pval: The nominal backward p-value of the association between the most significant variant and the phenotype.
library(qvalue)
cis$FDR <- qvalue(p = cis$bwd_pval)$qvalues
range(cis$FDR)

cis = cis[which(cis$FDR < 1e-5),]
cis = cis[which(cis$maf > 0.039),]
cis = cis[which(cis$`2` > 5),]
cis_sub <- cis[!(is.na(cis$phe_id) | cis$phe_id == "" ),]
cis_sub = cis_sub[order(as.numeric(cis_sub$FDR), decreasing = F),]
# cis_sub <- cis_sub[!duplicated(cis_sub$s),] 

cis_sub <- cis_sub[which(cis_sub$bwd_best_hit == 1),] 
dim(cis_sub)
#--

RESULT = data.frame()
i=1;
library(DGCA)
# 
for(i in 1:dim(cis_sub)[1]) {
  tryCatch({
  
  QTL_temp = cis_sub[i,] # cis[which(cis$var_id == 'rs2910792' & cis$phe_id == 'ERAP2'),]
  Gene1_expression = expression
  TP_GENOTYPE_sub = TP_GENOTYPE[which(row.names(TP_GENOTYPE) == QTL_temp$var_id),] 
  
  #---------- Ref homozygote
  library(ggpubr)
  library(data.table)
  expr0=Gene1_expression[, which(colnames(Gene1_expression) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==0)])]
  expr1=Gene1_expression[, which(colnames(Gene1_expression) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==1)])]
  expr2=Gene1_expression[, which(colnames(Gene1_expression) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==2)])]
  #------------------------------------
  
  set.seed(1234)
  design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr0)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr0)[2])))
  colnames(design_matrix) <- c("G2", "G0")
  INPUT = cbind(expr2, expr0)
  ddcor_G2_G0 = ddcorAll(inputMat = INPUT, design = design_matrix,  compare = c("G2","G0"),  adjust = "fdr", nPerm = 0, corrType = "pearson", splitSet = QTL_temp$phe_id, sigThresh = 0.05, sortBy = "pValDiff_adj", verbose = TRUE)
  
  design_matrix <- cbind(rep(c(1,0), c(dim(expr2)[2],dim(expr1)[2])), rep(c(0,1), c(dim(expr2)[2],dim(expr1)[2])))
  colnames(design_matrix) <- c("G2", "G1")
  INPUT = cbind(expr2, expr1)
  ddcor_G2_G1 = ddcorAll(inputMat = INPUT, design = design_matrix,  compare = c("G2","G1"),  
             adjust = "fdr", nPerm = 0, corrType = "pearson", splitSet = QTL_temp$phe_id, sigThresh = 0.05, sortBy = "pValDiff_adj", verbose = TRUE)
  
  ddcor_G2_G0$index = 'G2_G0'
  ddcor_G2_G1$index = 'G2_G1'
  
  ddcor_G2_G0[,12] = QTL_temp$var_id 
  ddcor_G2_G0[,13] = QTL_temp$REF 
  ddcor_G2_G0[,14] = QTL_temp$ALT 
  colnames(ddcor_G2_G0) = c('Gene1','Gene2','GM_cor','GM_pVal','GR_cor','GR_pVal','zScoreDiff','pValDiff','pValDiff_adj','Classes','index','var_id','REF','ALT')
  
  ddcor_G2_G1[,12] = QTL_temp$var_id 
  ddcor_G2_G1[,13] = QTL_temp$REF 
  ddcor_G2_G1[,14] = QTL_temp$ALT 
  colnames(ddcor_G2_G1) = c('Gene1','Gene2','GM_cor','GM_pVal','GR_cor','GR_pVal','zScoreDiff','pValDiff','pValDiff_adj','Classes','index','var_id','REF','ALT')
  
  ddcor_G2_G1 = ddcor_G2_G1[which(ddcor_G2_G1$Classes != 'NonSig'),]
  ddcor_G2_G0 = ddcor_G2_G0[which(ddcor_G2_G0$Classes != 'NonSig'),]
  
  dim(ddcor_G2_G1)
  dim(ddcor_G2_G0)
  
  if(i==1){RESULT = rbind(ddcor_G2_G0, ddcor_G2_G1)}else{RESULT = rbind(RESULT, rbind(ddcor_G2_G0, ddcor_G2_G1))}
  
  print(unique(RESULT$Gene2))
  
  # if(dim(RESULT)[1] != 0 )
  # {
    #--- save results
    write.table(RESULT, paste0('coExQTL/gQTL_coExQTL/', treatment,'_coExQTL_gQTL.txt'), quote = F, row.names = F, sep = '\t')
  # } 

  print(i)
  print(dim(RESULT))
  
  }, error=function(e){})
}
