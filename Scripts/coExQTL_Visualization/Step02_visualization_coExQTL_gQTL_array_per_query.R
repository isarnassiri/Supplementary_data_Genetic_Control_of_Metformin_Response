library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
head(listAttributes(ensembl))
filters = listFilters(ensembl)
GTF <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(GTF) = c('symbol', 'id')

library(data.table)
SNPs = fread('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'

ListSNP = fread('/Users/isarnassiri/Documents/RESULTS_USED/QTLtools_input_files_for_trans_analysis/QTLtools/Genotypes_metaData/list_SNPs.txt',header = T, stringsAsFactors = F)
dim(ListSNP)

#--------------------------- query ---------------------------
setwd('/Users/isarnassiri/Documents/RESULTS_USED')
treat = c('IFN', 'LPS24', 'UT')
i=1

input1 = data.frame(var_id='rs7305461', Gene2='RPS26', Gene1='KPNA2', stringsAsFactors = F) 
input2 = data.frame(var_id='rs3110426', Gene2='OXR1', Gene1='PTGS2', stringsAsFactors = F) 
input3 = data.frame(var_id='rs2910789', Gene2='ERAP2', Gene1='EPM2AIP1', stringsAsFactors = F) 

for(i in 1:3)
{
  for(inp in 1:3)
  {
    print(inp)
    input = get(paste0('input', inp))
    visualization()
  }
}

#-------------------

visualization = function()
{
  
  #--
  Gene1=input$Gene1
  Gene2=input$Gene2 
  SNP=input$var_id 
  
  REF=SNPs$REF[which(SNPs$SNP_ID == SNP)]
  ALT=SNPs$ALT[which(SNPs$SNP_ID == SNP)]
  #--
  
  print(treat[i])
  treatment = treat[i]
  
  gQTL_coExQTL = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/gQTL_coExQTL/', treat[i], '_coExQTL_array_gQTL.txt'), stringsAsFactors = F, header = T)
  gQTL_coExQTL = as.data.frame(gQTL_coExQTL)
  
  # as I do not use permutation, so I need to perform multiple testing correction
  library(qvalue)
  gQTL_coExQTL$FDR <- qvalue(gQTL_coExQTL$pValDiff, lambda = 0)$qvalues
  
  setwd('/Users/isarnassiri/Documents/RESULTS_USED/QTLtools_input_files_for_trans_analysis/QTLtools/')
  library(data.table)
  
  if(treatment=='UT'){TP_EXPESSION = fread(paste0('expression_CD14_array_FF14.bed.gz'), stringsAsFactors = F, header = T)
  }else{TP_EXPESSION = fread(paste0('expression_',treatment,'_array_FF14.bed.gz'), stringsAsFactors = F, header = T)
  }
  
  TP_EXPESSION = as.data.frame(TP_EXPESSION)
  row.names(TP_EXPESSION) = make.names(TP_EXPESSION$pid, unique=T)
  TP_EXPESSION = TP_EXPESSION[,-c(1:6)]
  dim(TP_EXPESSION)
  
  #NOTE: first row of GENOTYPE is an index for columns, so I use which(ListSNP$ID %in% SNP)+1
  GENOTYPE = fread(paste0('Genotypes_metaData/',treatment,'/', treatment,'_metaData_012.012'), select=(which(ListSNP$ID %in% SNP)+1), stringsAsFactors = F)
  TP_GENOTYPE = t(GENOTYPE)
  dim(TP_GENOTYPE)
  
  GENOTYPE_sample_names = fread(paste0('Genotypes_metaData/',treatment,'/',treatment,'_metaData_012.012.indv'), stringsAsFactors = F,header = F)
  dim(GENOTYPE_sample_names)
  
  colnames(TP_GENOTYPE) = GENOTYPE_sample_names$V1
  TP_GENOTYPE = data.frame(ID = SNP, TP_GENOTYPE)
  row.names(TP_GENOTYPE) = NULL
  colnames(TP_GENOTYPE) = gsub('X','',colnames(TP_GENOTYPE))
  colnames(TP_GENOTYPE)
  
  TP_EXPESSION = as.data.frame(TP_EXPESSION)
  TP_GENOTYPE = as.data.frame(TP_GENOTYPE)
  
  TP_EXPESSION = TP_EXPESSION[,which(colnames(TP_EXPESSION) %in% colnames(TP_GENOTYPE))]
  TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(TP_EXPESSION))]
  dim(TP_GENOTYPE)
  dim(TP_EXPESSION)
  
  library(stringr)
  
  # Gene1_id=unique(GTF[which(GTF$symbol == Gene1),'id'])
  # Gene2_id=unique(GTF[which(GTF$symbol == Gene2),'id'])
  
  TP_EXPESSION_Gene1 = TP_EXPESSION[which(row.names(TP_EXPESSION) %in% Gene1),]
  dim(TP_EXPESSION_Gene1) 
  
  TP_EXPESSION_Gene2 = TP_EXPESSION[which(row.names(TP_EXPESSION) %in% Gene2),]
  dim(TP_EXPESSION_Gene2) 
  
  TP_GENOTYPE_sub = TP_GENOTYPE#[which(row.names(TP_GENOTYPE) == SNP),]
  dim(TP_GENOTYPE_sub)
  
  TP_EXPESSION_Gene1 = as.data.frame(TP_EXPESSION_Gene1)
  TP_EXPESSION_Gene2 = as.data.frame(TP_EXPESSION_Gene2)
  TP_GENOTYPE_sub = as.data.frame(TP_GENOTYPE_sub)
  TP_GENOTYPE_sub
  
  #=============== visualization type 1
  #---------- Ref homozygote
  
  expr0=data.frame(Gene1=as.numeric(TP_EXPESSION_Gene1[,which(TP_GENOTYPE_sub[1,]==0)]), 
                   Gene2=as.numeric(TP_EXPESSION_Gene2[,which(TP_GENOTYPE_sub[1,]==0)]) )
  expr1=data.frame(Gene1=as.numeric(TP_EXPESSION_Gene1[,which(TP_GENOTYPE_sub[1,]==1)]), 
                   Gene2=as.numeric(TP_EXPESSION_Gene2[,which(TP_GENOTYPE_sub[1,]==1)]) )
  expr2=data.frame(Gene1=as.numeric(TP_EXPESSION_Gene1[,which(TP_GENOTYPE_sub[1,]==2)]), 
                   Gene2=as.numeric(TP_EXPESSION_Gene2[,which(TP_GENOTYPE_sub[1,]==2)]) )
  
  expr0$Genotype = rep(paste0('0: ' ,REF,'_', REF), dim(expr0)[1])
  expr1$Genotype = rep(paste0('1: ' ,REF,'_', ALT), dim(expr1)[1])
  expr2$Genotype = rep(paste0('2: ' ,ALT,'_', ALT), dim(expr2)[1])
  
  input.plot = do.call(rbind,list(expr0,expr1,expr2))
  
  input.plot$Genotype = factor(input.plot$Genotype, levels = unique(input.plot$Genotype))
  
  library(ggpubr)
  library(ggplot2)
  
  theme_set(theme_bw())# pre-set the bw theme.
  g <- ggplot(input.plot, aes( Gene1, Gene2)) + labs(subtitle="Allele-specific correlation for expresseion values per genotype", title="", x = paste0(Gene1, ' Expression'), y = paste0(Gene2, ' Expression'), size="Expression", color = SNP)
  
  p = g + geom_jitter(aes(col=Genotype, size=Gene1)) +
    geom_smooth(aes(col=Genotype), method="lm", se=F) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + geom_rug(aes(color = Genotype)) +
    theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.subtitle = element_text(size = 25 ),legend.title=element_text(size=25), legend.text=element_text(size=20))#+ facet_wrap(~Genotype) , face="bold"
  
  scientific_10 <- function(x, ...) {
    parse(text = gsub("e", "%*%10^", scales::label_scientific(...)(x)))
  }
  
  FDR = gQTL_coExQTL[which(gQTL_coExQTL$Gene1 == Gene1 & gQTL_coExQTL$Gene2 == Gene2 & gQTL_coExQTL$var_id == SNP), c('FDR', 'index')]
  
  if(length(FDR$FDR) == 0)
  {
    dat_text <- data.frame(
      label = paste0('DC-FDR > 0.05')
    )
  }else{
    
    FDR = FDR[order(FDR$FDR, decreasing = F),]
    
    dat_text <- data.frame(
      label = paste0('DC-FDR: ', scientific_10(FDR$FDR[1] , digits = 3) )
    )
  }
  
  print(dat_text)
  print(input)
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
  
  if(!dir.exists(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gQTL_coExQTL_array/')))
  {
    dir.create(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gQTL_coExQTL_array/'), recursive = T)
  }
  
  library(grDevices)
  library(grid)
  pdf(file=paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gQTL_coExQTL_array/', Gene1,'_', Gene2,'_', SNP,'_', treatment,"_array_gcoExQTL.pdf"), width = 12, height = 10, useDingbats = F)
  print(p)
  dev.off()
  
}
