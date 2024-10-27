
#-------------------
umask 002
module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

srun -p short --pty bash
#-

.libPaths("/well/parkkinen/users/gbf362/R/4.1/skylake/")

treatment="IFN"
setwd('/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/')

library(data.table)
library(SAVER)
library(dplyr)
expression = fread(paste0('Input_files_transcript/expression_',treatment,'.txt'), stringsAsFactors = F)
expression = as.data.frame(expression)
row.names(expression) = expression$id
expression = expression[,-1]

MAF_imputed_Allsamples = fread('MAF_imputed_Allsamples_TYPED.txt', stringsAsFactors = F, header = T)
MAF_imputed_Allsamples = as.data.frame(MAF_imputed_Allsamples)
colnames(MAF_imputed_Allsamples)[which(colnames(MAF_imputed_Allsamples) == 'ID')] = 'var_id'

headerC = fread('header_conditional_tQTL.txt', stringsAsFactors = F, header = F)
cis = fread(paste0('tQTL/QTLtools_outputs-100kb-window/',treatment,'/tQTL_conditional_pass_1.txt'), stringsAsFactors = F)
cis = as.data.frame(cis)

colnames(cis) = as.character(headerC)
cis = cis[-which(cis$var_id == '.'),]
cis = merge(cis, MAF_imputed_Allsamples, by = 'var_id')
dim(cis)

# bwd_pval: The nominal backward p-value of the association between the most significant variant and the phenotype.
library(qvalue)
cis$FDR <- qvalue(p = cis$bwd_pval)$qvalues
range(cis$FDR)

cis = cis[which(cis$FDR < 1e-5),]
cis = cis[which(cis$maf > 0.039),]
cis = cis[which(cis$`2` > 5),]
cis <- cis[!(is.na(cis$phe_id) | cis$phe_id == "" ),]

cis = cis[order(as.numeric(cis$FDR), decreasing = F),]
cis_sub = cis[!duplicated(cis$phe_id),]
#------------

set.seed(123)
expression_saver <- saver(expression[which(row.names(expression) %in% cis_sub$phe_id),], ncores = 6, estimates.only = TRUE)
expression_saver = data.frame(GeneID=row.names(expression_saver), expression_saver)

fwrite(expression_saver, paste0('Input_files_transcript/expression_',treatment,'_SAVER.txt'), row.names = F, quote = F, sep = '\t')
