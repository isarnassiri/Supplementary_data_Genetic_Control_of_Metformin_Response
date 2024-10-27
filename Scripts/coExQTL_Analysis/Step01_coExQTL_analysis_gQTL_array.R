module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
#srun -p short --pty bash

.libPaths("/well/parkkinen/users/gbf362/R/4.1/skylake/")

#-------------------
setwd('/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/')

treatment="UT"

library(data.table)
if(treatment=='UT'){TP_EXPESSION = fread(paste0('Input_files_gene/expression_CD14_array_FF14.bed.gz'), stringsAsFactors = F, header = T)
}else{TP_EXPESSION = fread(paste0('Input_files_gene/expression_',treatment,'_array_FF14.bed.gz'), stringsAsFactors = F, header = T)
}

TP_EXPESSION = as.data.frame(TP_EXPESSION)
row.names(TP_EXPESSION) = make.names(TP_EXPESSION$pid, unique=T)
TP_EXPESSION = TP_EXPESSION[,-c(1:6)]
dim(TP_EXPESSION)

GENOTYPE = fread(paste0('Input_files_Genotype/Genotypes_metaData/', treatment,'_array_genotype_012.txt'), stringsAsFactors = F)
GENOTYPE = as.data.frame(GENOTYPE)

row.names(GENOTYPE) = make.names(GENOTYPE$ID, unique = T)
GENOTYPE = GENOTYPE[,-1]

TP_EXPESSION = TP_EXPESSION[,which(colnames(TP_EXPESSION) %in% colnames(GENOTYPE))]
GENOTYPE = GENOTYPE[,which(colnames(GENOTYPE) %in% colnames(TP_EXPESSION))]
dim(GENOTYPE)
dim(TP_EXPESSION)

#---
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

cis_sub <- cis_sub[which(cis_sub$bwd_best_hit == 1),] 
dim(cis_sub)
#--

RESULT = data.frame()
i=1;
library(DGCA)

for(i in 1:dim(cis_sub)[1]) {
 tryCatch({
 
 QTL_temp = cis_sub[i,] # cis[which(cis$var_id == 'rs2910792' & cis$phe_id == 'ERAP2'),]
 GENOTYPE_sub = GENOTYPE[which(row.names(GENOTYPE) == QTL_temp$var_id),] 
 
 #---------- Ref homozygote
 library(ggpubr)
 library(data.table)
 expr0=TP_EXPESSION[, which(colnames(TP_EXPESSION) %in% colnames(GENOTYPE_sub)[which(GENOTYPE_sub[1,]==0)])]
 expr1=TP_EXPESSION[, which(colnames(TP_EXPESSION) %in% colnames(GENOTYPE_sub)[which(GENOTYPE_sub[1,]==1)])]
 expr2=TP_EXPESSION[, which(colnames(TP_EXPESSION) %in% colnames(GENOTYPE_sub)[which(GENOTYPE_sub[1,]==2)])]
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
 # Note that you should also set adjust = “perm”, because otherwise a different p-value adjustment method will be used, and the time spent generating the permutation samples will have been wasted.
 # https://htmlpreview.github.io/?https://github.com/andymckenzie/DGCA/blob/master/inst/doc/DGCA.html
 
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
 write.table(RESULT, paste0('coExQTL/gQTL_coExQTL/', treatment,'_coExQTL_array_gQTL.txt'), quote = F, row.names = F, sep = '\t')
 # } 
 
 print(i)
 print(dim(RESULT))
 
 }, error=function(e){})
}

#======================== prepration of array expression profiles

# library(data.table)
# setwd('/t1-data/data/fairfaxlab/transcriptomics/RNAseq/eQTL/monocytes/analysis/QTL_analysis/Fastqtl/RESULTS_USED/Kaur_results/Scipaper_2014_expression_profile/')
# expression_all = fread('Fairfax_2014.HumanHT-12_V4_norm_exprs.tsv', stringsAsFactors = F)
# expression_all = as.data.frame(expression_all)
# 
# probeID=expression_all$phenotype_id
# library("illuminaHumanv4.db")  
# gene = data.frame(pid=unlist(mget(x = probeID,envir = illuminaHumanv4SYMBOL))) 
# 
# expression_all = cbind(gene, expression_all)
# expression_all = expression_all[,-2]
# expression_all[1:10,1:10]
# fwrite(expression_all,'/t1-data/data/fairfaxlab/transcriptomics/RNAseq/eQTL/monocytes/analysis/QTL_analysis/co-expression_QTL/co-expression_array/Fairfax_2014_norm_GeneName.txt', quote = F, row.names = F, sep = '\t')
# #---------------------------
# 
# library(data.table)
# library(rtracklayer)
# gtf <- rtracklayer::import('/t1-data/user/nassisar/HISAT/GENOMEANNOT/GRCh38.gtf')
# gtf_df <- as.data.frame(gtf)
# gtf_df <- gtf_df[which(gtf_df$type == 'gene'),]
# gtf_df = gtf_df[,c('seqnames', 'start', 'end', 'gene_name', 'gene_id', 'strand')]
# dim(gtf_df)
# colnames(gtf_df) = c('Chr', 'start', 'end', 'pid', 'gid', 'strand')
# 
# #------------- expression array
# setwd('/t1-data/data/fairfaxlab/transcriptomics/RNAseq/eQTL/monocytes/analysis/QTL_analysis/co-expression_QTL/co-expression_array/')
# 
# expression_all = fread('Fairfax_2014_norm_GeneName.txt', stringsAsFactors = F)
# expression_all = as.data.frame(expression_all)
# dim(expression_all)
# 
# treatment = c('IFN', 'LPS24', 'CD14')
# for(i in 1:3)
# {
#   Expression = expression_all[,grep(treatment[i], colnames(expression_all))]
#   colnames(Expression) = gsub(paste0(treatment[i],'_'),'',colnames(Expression))
#   Expression = Expression[,order(as.numeric(colnames(Expression)),decreasing = F)]
#   Expression = data.frame(id = expression_all$pid,Expression, stringsAsFactors = F)
#   dim(Expression)
#   colnames(Expression) = gsub('X','',colnames(Expression))
#  print( dim(Expression))
#   fwrite(Expression,paste0('expression_',treatment[i],'.txt'), quote = F, row.names = F, sep = '\t')
# }
