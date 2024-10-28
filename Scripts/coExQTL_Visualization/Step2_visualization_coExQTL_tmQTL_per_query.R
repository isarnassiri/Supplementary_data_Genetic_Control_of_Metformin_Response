
#======================================= visualization of interaction - per genotype and per all
library(data.table)
treatment='UT'
  
setwd('/Users/isarnassiri/Documents/RESULTS_USED/')

library(dplyr)
expression=fread(paste0('Expression/Input_files_transcript/expression_',treatment,'_withtQTL_SAVER.txt'),stringsAsFactors=F)
expression=as.data.frame(expression)
row.names(expression)=expression$GeneID
expression=expression[,-1]
colnames(expression)=gsub('X','',colnames(expression))

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
query_all = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/tmQTL_coExQTL/',treatment,'_coExQTL_tmQTL.txt'), stringsAsFactors = F, header = T)
query_all = as.data.frame(query_all)
query_all = query_all[order(query_all$pValDiff_adj, decreasing = F),]
dim(query_all)

library(qvalue)
query_all$FDR <- qvalue(query_all$pValDiff, lambda = 0)$qvalues
range(query_all$FDR)

length(unique(query_all$Gene1))
query = query_all[!duplicated(query_all$Gene2_name),]
#query = query_all[which(query_all$Gene2_name == 'ELP5' & query_all$Gene2 == 'ENST00000574841'),] # 'NFKB1'
query = query_all[which(query_all$Gene2_name == 'NFKB1'),] # 'USMG5'

dim(query)
dim(TP_GENOTYPE)
dim(expression_methylation)

y=1
for (y in 1:dim(query)[1]) {
  
  Gene1=query$Gene1[y]
  Gene1_id=query$Gene1[y]
  
  Gene2=query$Gene2[y]
  Gene2_id=query$Gene2_name[y]
  
  SNP_ID = query$var_id[y]
  
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
  
  TP_EXPESSION_Gene1 = as.data.frame(TP_EXPESSION_Gene1)
  TP_EXPESSION_Gene2 = as.data.frame(TP_EXPESSION_Gene2)
  TP_GENOTYPE_sub = as.data.frame(TP_GENOTYPE_sub)
  #---
  
  #--- remove outliers - Methyaltion
  Q <- quantile(as.vector(t(TP_EXPESSION_Gene1[1,])), probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR( as.vector(t(TP_EXPESSION_Gene1[1,])) )
  up <-  Q[2]+1.5*iqr # Upper Range
  low<- Q[1]-1.5*iqr  # Lower Range
  
  if(length(which(as.vector(t(TP_EXPESSION_Gene1[1,])) < up & as.vector(t(TP_EXPESSION_Gene1[1,])) > low)) != 0)
  {
    TP_EXPESSION_Gene1 = TP_EXPESSION_Gene1[,which(as.vector(t(TP_EXPESSION_Gene1[1,])) < up & as.vector(t(TP_EXPESSION_Gene1[1,])) > low)]
  }
  
  dim(TP_EXPESSION_Gene1)
  #---
  
  #--- remove outliers - Expression
  Q <- quantile(as.vector(t(TP_EXPESSION_Gene2[1,])), probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR( as.vector(t(TP_EXPESSION_Gene2[1,])) )
  up <-  Q[2]+1.5*iqr # Upper Range
  low<- Q[1]-1.5*iqr  # Lower Range

  if(length(which(as.vector(t(TP_EXPESSION_Gene2[1,])) < up & as.vector(t(TP_EXPESSION_Gene2[1,])) > low)) != 0 & !is.null(dim(TP_EXPESSION_Gene2[,which(as.vector(t(TP_EXPESSION_Gene2[1,])) < up & as.vector(t(TP_EXPESSION_Gene2[1,])) > low)])) )
  {
    TP_EXPESSION_Gene2 = TP_EXPESSION_Gene2[,which(as.vector(t(TP_EXPESSION_Gene2[1,])) < up & as.vector(t(TP_EXPESSION_Gene2[1,])) > low)]
  }

  dim(TP_EXPESSION_Gene2)
  #---

  TP_EXPESSION_Gene1 = TP_EXPESSION_Gene1[,which(colnames(TP_EXPESSION_Gene1) %in% colnames(TP_EXPESSION_Gene2))]
  TP_EXPESSION_Gene2 = TP_EXPESSION_Gene2[,which(colnames(TP_EXPESSION_Gene2) %in% colnames(TP_EXPESSION_Gene1))]
  dim(TP_EXPESSION_Gene1)
  dim(TP_EXPESSION_Gene2)

  TP_GENOTYPE_sub = TP_GENOTYPE_sub[,which(colnames(TP_GENOTYPE_sub) %in% colnames(TP_EXPESSION_Gene2))]
  dim(TP_GENOTYPE_sub)
  
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
  
  #--- remove outliers
  # Q <- quantile(c(input.plot$Gene1, input.plot$Gene2), probs=c(.25, .75), na.rm = FALSE)
  # iqr <- IQR(c(input.plot$Gene1, input.plot$Gene2))
  # up <-  as.double(Q[2]+1.5*iqr) # Upper Range
  # low <- as.double(Q[1]-1.5*iqr)  # Lower Range
  # 
  # input.plot = input.plot[which(input.plot$Gene1 < up & input.plot$Gene2 < up & input.plot$Gene1 > low & input.plot$Gene2 > low),]
  #---
  
  library(ggpubr)
  library(ggplot2)
  
  theme_set(theme_bw())# pre-set the bw theme.
  g <- ggplot(input.plot, aes( Gene1, Gene2)) + labs(subtitle="Allele-specific correlation for expresseion values per genotype", title="", x = paste0(Gene1, ' Methylation'), y = paste0(Gene2_id, ' (', Gene2, ')', ' Expression'), size="Expression", color = SNP_ID)
  
  p = g + geom_jitter(aes(col=Genotype, size=Gene1)) +
    geom_smooth(aes(col=Genotype), method="lm", se=F) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + geom_rug(aes(color = Genotype)) +
    theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.subtitle = element_text(size = 25 ),legend.title=element_text(size=25), legend.text=element_text(size=20))#+ facet_wrap(~Genotype) , face="bold"
  
  scientific_10 <- function(x, ...) {
    parse(text = gsub("e", "%*%10^", scales::label_scientific(...)(x)))
  }
  
  FDR = query$FDR[y]
  
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
  
  if(!dir.exists(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tmQTL_coExQTL/')))
  {
    dir.create(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tmQTL_coExQTL/'), recursive = T)
  }
  
  library(grDevices)
  library(grid)
  pdf(file=paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tmQTL_coExQTL/', Gene1,'_', Gene2,'_', Gene2_id, '_', SNP_ID,'_', treatment,"_tmcoExQTL.pdf"), width = 12, height = 10, useDingbats = F)
  print(p)
  dev.off()
   
}

