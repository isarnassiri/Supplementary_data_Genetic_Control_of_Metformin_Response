
#======================================= visualization of interaction - per genotype and per all
library(data.table)
treatment='LPS24'

#--- gene symbol to Ensemble transcript IDs
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
head(listAttributes(ensembl))
filters = listFilters(ensembl)
GTF <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(GTF) = c('symbol', 'id')

setwd('/Users/isarnassiri/Documents/RESULTS_USED/')

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

identical(row.names(expression), SYMBOL$id)

row.names(expression) = SYMBOL$symbol

TP_GENOTYPE = fread(paste0('Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),stringsAsFactors = F)
TP_GENOTYPE = as.data.frame(TP_GENOTYPE)
TP_GENOTYPE = TP_GENOTYPE[!duplicated(TP_GENOTYPE$id),]
row.names(TP_GENOTYPE) = TP_GENOTYPE$id
TP_GENOTYPE = TP_GENOTYPE[,-1]

expression = expression[,which(colnames(expression) %in% colnames(TP_GENOTYPE))]
TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression))]
identical(colnames(TP_GENOTYPE), colnames(expression))

METHYLATION = fread(paste0('Expression/Input_files_methylation/Methylation_',treatment,'_hg38_QTLtools.bed.gz'), stringsAsFactors = F, header = T)
METHYLATION = as.data.frame(METHYLATION)
METHYLATION = METHYLATION[,-c(1:3)]
colnames(METHYLATION)[1] = 'id'
colnames(METHYLATION)[-1] = sub('X','',colnames(METHYLATION)[-1])
row.names(METHYLATION) = METHYLATION[,1]
METHYLATION = METHYLATION[,which(colnames(METHYLATION) %in% colnames(expression))]
expression = expression[,which(colnames(expression) %in%  colnames(METHYLATION))]

dim(METHYLATION)

identical(colnames(METHYLATION), colnames(expression))
expression_methylation = rbind(METHYLATION, expression)

expression_methylation = expression_methylation[,which(colnames(expression_methylation) %in% colnames(TP_GENOTYPE))]
TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression_methylation))]
identical(colnames(expression_methylation), colnames(TP_GENOTYPE))

#-------------------

library(data.table)
query = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/gmQTL_coExQTL/',treatment,'_coExQTL_gmQTL.txt'), stringsAsFactors = F)
query = as.data.frame(query)

query = query[order(query$pValDiff, decreasing = F),]

library(qvalue)
query$FDR <- qvalue(query$pValDiff, lambda = 0)$qvalues

query = query[which(query$Gene2_name == 'OXR1' & query$var_id == 'rs3110426'),]
dim(query)

query = query[1:10,]

# query = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/eQTL_mQTL_pairs_analysis/eQTL_mQTL_corr_CS_coloc_selected.txt', stringsAsFactors = F, header = T)
# query = as.data.frame(query)
# query = query[which(query$gene_name == 'DTX4' & query$cpg_ID == 'cg07745373' & query$SNP_ID == 'rs7934971'),]
# query

#-------------------
dim(query)
dim(TP_GENOTYPE)
dim(expression_methylation)
y=1

for (y in 1:dim(query)[1]) {
  
  library(stringr)
  INPUT_genes = data.frame(V1=query$Gene1[y], V2=query$Gene2_name[y], V3 = query$var_id[y], REF=query$REF[y], ALT =query$ALT[y], stringsAsFactors = F)
 
  Gene1=INPUT_genes$V1

  Gene2=INPUT_genes$V2

  SNP_ID = INPUT_genes$V3
  
  REF=query$REF[y]
  ALT =query$ALT[y]
  
  TP_EXPESSION_Gene1 = expression_methylation[which(row.names(expression_methylation) %in% Gene1),]
  dim(TP_EXPESSION_Gene1) 
  
  TP_EXPESSION_Gene2 = expression_methylation[which(row.names(expression_methylation) %in% Gene2),]
  dim(TP_EXPESSION_Gene2) 
  
  TP_GENOTYPE_sub = TP_GENOTYPE[which(row.names(TP_GENOTYPE) == SNP_ID),]
  dim(TP_GENOTYPE_sub)
  
  TP_EXPESSION_Gene1 = data.table(TP_EXPESSION_Gene1)
  TP_EXPESSION_Gene2 = data.table(TP_EXPESSION_Gene2)
  TP_GENOTYPE_sub = data.table(TP_GENOTYPE_sub)
  
  #=============== visualization type 1
  #---------- Ref homozygote
  library(ggpubr)
  library(data.table)
  expr0=data.frame(Gene1=as.numeric(select(TP_EXPESSION_Gene1, which(colnames(TP_EXPESSION_Gene1) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==0)]))), 
                   Gene2=as.numeric(select(TP_EXPESSION_Gene2, which(colnames(TP_EXPESSION_Gene2) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==0)]))) )
  expr1=data.frame(Gene1=as.numeric(select(TP_EXPESSION_Gene1, which(colnames(TP_EXPESSION_Gene1) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==1)]))), 
                   Gene2=as.numeric(select(TP_EXPESSION_Gene2, which(colnames(TP_EXPESSION_Gene2) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==1)]))) )
  expr2=data.frame(Gene1=as.numeric(select(TP_EXPESSION_Gene1, which(colnames(TP_EXPESSION_Gene1) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==2)]))), 
                   Gene2=as.numeric(select(TP_EXPESSION_Gene2, which(colnames(TP_EXPESSION_Gene2) %in% colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,]==2)]))) )
  
  expr0$Genotype = rep(paste0('0: ' ,REF,'_', REF), dim(expr0)[1])
  expr1$Genotype = rep(paste0('1: ' ,REF,'_', ALT), dim(expr1)[1])
  expr2$Genotype = rep(paste0('2: ' ,ALT,'_', ALT), dim(expr2)[1])
  
  input.plot = do.call(rbind,list(expr0,expr1,expr2))
  
  input.plot$Genotype = factor(input.plot$Genotype, levels = unique(input.plot$Genotype))
  
  library(ggpubr)
  library(ggplot2)
  
  theme_set(theme_bw())# pre-set the bw theme.
  g <- ggplot(input.plot, aes( Gene1, Gene2)) + labs(subtitle="Allele-specific correlation for expresseion values per genotype", title="", x = paste0(Gene1, ' Methylation'), y = paste0(Gene2, ' Expression'), size="Expression", color = SNP_ID)
  
  p = g + geom_jitter(aes(col=Genotype, size=Gene1)) +
    geom_smooth(aes(col=Genotype), method="lm", se=F) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + geom_rug(aes(color = Genotype)) +
    theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.subtitle = element_text(size = 25 ),legend.title=element_text(size=25), legend.text=element_text(size=20))#+ facet_wrap(~Genotype) , face="bold"
  
  scientific_10 <- function(x, ...) {
    parse(text = gsub("e", "%*%10^", scales::label_scientific(...)(x)))
  }
  
  FDR = query$FDR
  
  if(length(FDR) == 0)
  {
    dat_text <- data.frame(
      label = paste0('Correlation-FDR > 0.05')
    )
  }else{
    
    FDR = FDR[order(FDR, decreasing = F)]
    
    dat_text <- data.frame(
      label = paste0('Correlation-FDR: ', scientific_10(FDR[1] , digits = 3) )
    )
  }
  
  print(dat_text)
  print(FDR)
  
  p <- p +
    geom_text(
      size    = 8,
      data    = dat_text,
      mapping = aes(x = Inf, y = Inf, label = label),
      hjust   = 1.05,
      vjust   = 1.5,
      parse = TRUE 
    )
  
  print(p)
  
  if(!dir.exists(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gmQTL_coExQTL_gQTL/')))
  {
    dir.create(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gmQTL_coExQTL_gQTL/'), recursive = T)
  }
  
  library(grDevices)
  library(grid)
  pdf(file=paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gmQTL_coExQTL_gQTL/', Gene1,'_', Gene2,'_', SNP_ID,'_', treatment, "_gQTL_gcoExQTL.pdf"), width = 12, height = 10, useDingbats = F)
  print(p)
  dev.off()
  
}







 

