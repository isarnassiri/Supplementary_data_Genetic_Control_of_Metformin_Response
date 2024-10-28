######################## visualization of tQTL ########################

#============================== co-expression tQTL ==============================

library(data.table)
library(dplyr)
SNPs = fread('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'

#--- gene symbol to Ensemble transcript IDs
library(rtracklayer)
GTF <- rtracklayer::import('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/tQTL/QTLtools_inputs/GRCh38.gtf')
GTF <- as.data.frame(GTF)

# if you want to read an SNP you should make the file tab separated
TP_GENOTYPE = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),stringsAsFactors = F)
TP_GENOTYPE = as.data.frame(TP_GENOTYPE)
TP_GENOTYPE = TP_GENOTYPE[!duplicated(TP_GENOTYPE$id),]
row.names(TP_GENOTYPE) = TP_GENOTYPE$id
TP_GENOTYPE = TP_GENOTYPE[,-1] 


#============================== input
i=3
treat = c('IFN', 'LPS24', 'UT')

listFiles = list.files(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/Selected_USED/',treat[i],'/'), pattern = 'good')

listFiles = gsub('_RNAseq.*', '', listFiles)
listFiles = gsub('__', '_', listFiles)

library(stringr)
listFiles = data.frame(str_split_fixed(listFiles, "_", str_count(listFiles[1], "_")+1 ))
colnames(listFiles) = c("Gene1", "Gene2", "Gene1_name", "Gene2_name", "var_id","treat")

merged_subset = listFiles

#==============================

q=1

for (q in 1:dim(merged_subset)[1]) {
  
caseStudy = FALSE

if(!caseStudy)
{
  Gene1=merged_subset[q,"Gene1"]
  Gene2=merged_subset[q,"Gene2"]
  Gene1_name=merged_subset[q,"Gene1_name"]
  Gene2_name=merged_subset[q,"Gene2_name"]
  SNP=merged_subset[q,"var_id"]
  
  REF=SNPs$REF[which(SNPs$SNP_ID == SNP)]
  ALT=SNPs$ALT[which(SNPs$SNP_ID == SNP)]

}

if(caseStudy)
{
  #ENST00000529254__ENST00000468864_ARL2_TRAPPC10_rs2243999_IFN_RNAseq_coExQTL_good
  
  Gene1='ENST00000529254'
  Gene2='ENST00000468864'
  Gene1_name='ARL2'
  Gene2_name='TRAPPC10'
  SNP='rs2243999'

  REF=SNPs$REF[which(SNPs$SNP_ID == SNP)]
  ALT=SNPs$ALT[which(SNPs$SNP_ID == SNP)]
}

for(i in 1:length(treat))
{
#-------------------

treatment=treat[i]
setwd('/Users/isarnassiri/Documents/RESULTS_USED/')

print(treat[i])

tQTL_coExQTL = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/tQTL_coExQTL_subset_of_TP/', treat[i], '_coExQTL_tQTL_R2.txt'), stringsAsFactors = F, header = T)
tQTL_coExQTL = as.data.frame(tQTL_coExQTL)
tQTL_coExQTL = tQTL_coExQTL[-which(is.na(tQTL_coExQTL$pValDiff_adj_R2)), ]

# as I do not use permutation, so I need to perform multiple testing correction
library(qvalue)
tQTL_coExQTL$FDR <- qvalue(tQTL_coExQTL$pValDiff_adj_R2)$qvalues

if(dim(tQTL_coExQTL[which(tQTL_coExQTL$Gene1 == Gene1 & tQTL_coExQTL$Gene2 == Gene2 & tQTL_coExQTL$var_id == SNP), ])[1] != 0)
{
  expression = fread(paste0('Expression/Input_files_transcript/expression_',treatment,'_withtQTL_SAVER.txt'), stringsAsFactors = F)
  expression = as.data.frame(expression)
  row.names(expression) = expression$GeneID
  expression = expression[,-1] 
  colnames(expression) = gsub('X', '', colnames(expression))
  
  FDR = tQTL_coExQTL[which(tQTL_coExQTL$Gene1 == Gene1 & tQTL_coExQTL$Gene2 == Gene2 & tQTL_coExQTL$var_id == SNP), c('FDR', 'index')]

}else{
  expression = fread(paste0('Expression/Input_files_transcript/expression_',treatment,'_AllTP_SAVER.txt'), stringsAsFactors = F)
  expression = as.data.frame(expression)
  row.names(expression) = expression$GeneID
  expression = expression[,-1] 
  colnames(expression) = gsub('X', '', colnames(expression))
  
  FDR = data.frame(FDR = 1, index = 'NA')
}

expression = expression[which(row.names(expression) %in% GTF$transcript_id),]
GTF_subset = GTF[which(GTF$transcript_id %in% row.names(expression)),]
GTF_subset = GTF_subset[!duplicated(GTF_subset$transcript_id),]

dim(GTF_subset)
dim(expression)

GTF_subset = GTF_subset[match(row.names(expression), GTF_subset$transcript_id),]
identical(GTF_subset$transcript_id, row.names(expression))

expression = expression[,which(colnames(expression) %in% colnames(TP_GENOTYPE))]
TP_GENOTYPE_sub = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression))]

dim(expression)
dim(TP_GENOTYPE_sub)

library(stringr)
TP_EXPESSION_Gene1 = expression[which(row.names(expression) %in% Gene1),]
dim(TP_EXPESSION_Gene1) 

TP_EXPESSION_Gene2 = expression[which(row.names(expression) %in% Gene2),]
dim(TP_EXPESSION_Gene2) 

TP_GENOTYPE_sub = TP_GENOTYPE_sub[which(row.names(TP_GENOTYPE_sub) == SNP),]
dim(TP_GENOTYPE_sub)

TP_EXPESSION_Gene1 = as.data.frame(TP_EXPESSION_Gene1)
TP_EXPESSION_Gene2 = as.data.frame(TP_EXPESSION_Gene2)

TP_GENOTYPE_sub = as.data.frame(TP_GENOTYPE_sub)

#=============== visualization type 1
#--------------- Ref homozygote

expr0=data.frame(Gene1=as.numeric(TP_EXPESSION_Gene1[,which(TP_GENOTYPE_sub[1,]==0)]), 
 Gene2=as.numeric(TP_EXPESSION_Gene2[,which(TP_GENOTYPE_sub[1,]==0)]) )
expr1=data.frame(Gene1=as.numeric(TP_EXPESSION_Gene1[,which(TP_GENOTYPE_sub[1,]==1)]), 
 Gene2=as.numeric(TP_EXPESSION_Gene2[,which(TP_GENOTYPE_sub[1,]==1)]) )
expr2=data.frame(Gene1=as.numeric(TP_EXPESSION_Gene1[,which(TP_GENOTYPE_sub[1,]==2)]), 
 Gene2=as.numeric(TP_EXPESSION_Gene2[,which(TP_GENOTYPE_sub[1,]==2)]) )

expr0$Genotype = rep(paste0('0: ' ,REF,'/', REF), dim(expr0)[1])
expr1$Genotype = rep(paste0('1: ' ,REF,'/', ALT), dim(expr1)[1])
expr2$Genotype = rep(paste0('2: ' ,ALT,'/', ALT), dim(expr2)[1])

# In a VCF (Variant Call Format) file, the numbers 0, 1, and 2 refer to the number of alternate alleles present at a specific genomic position for a given individual.

input.plot = do.call(rbind,list(expr0,expr1,expr2))

input.plot$Genotype = factor(input.plot$Genotype, levels = unique(input.plot$Genotype))

#--- remove outliers
Q <- quantile(c(input.plot$Gene1, input.plot$Gene2), probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(c(input.plot$Gene1, input.plot$Gene2))
up <-  as.double(Q[2]+1.5*iqr) # Upper Range
low <- as.double(Q[1]-1.5*iqr)  # Lower Range

input.plot = input.plot[which(input.plot$Gene1 < up & input.plot$Gene2 < up & input.plot$Gene1 > low & input.plot$Gene2 > low),]
#---

library(ggpubr)
library(ggplot2)

theme_set(theme_bw())# pre-set the bw theme.
g <- ggplot(input.plot, aes( Gene1, Gene2)) + labs(subtitle="Allele-specific correlation for expresseion values per genotype", title="", x = paste0(Gene1_name, ' (', Gene1, ')'), y = paste0(Gene2_name, ' (', Gene2, ')'), size="Expression", color = SNP)

p = g + geom_jitter(aes(col=Genotype, size=Gene1)) +
geom_smooth(aes(col=Genotype), method="lm", se=F) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + geom_rug(aes(color = Genotype)) +
theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.subtitle = element_text(size = 25 ),legend.title=element_text(size=25), legend.text=element_text(size=20))#+ facet_wrap(~Genotype) , face="bold"

scientific_10 <- function(x, ...) {
  parse(text = gsub("e", "%*%10^", scales::label_scientific(...)(x)))
}

# if(FDR$index == "G2_G1")
# {
#   Geno = paste0(ALT,'/', ALT, ' and ', REF,'/', ALT)
# }
# 
# if(FDR$index == "G2_G0")
# {
#   Geno = paste0(ALT,'/', ALT, ' and ', REF,'/', REF)
# }

if(dim(FDR)[1] == 1)
{
  if(FDR$FDR != 1)
  {
    dat_text <- data.frame(
      label = paste0('DC-FDR: ', scientific_10(FDR$FDR , digits = 3) )
    )
  }
  
  if(FDR$FDR == 1)
  {
    dat_text <- data.frame(
      label = paste0('DC-FDR > 0.05')
    )
  }
}else{
  dat_text <- data.frame(
    label = paste0('DC-FDR: ', scientific_10(FDR$FDR[1] , digits = 3) )
  )
}

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

if(!dir.exists(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tQTL_coExQTL/')))
{
  dir.create(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tQTL_coExQTL/'), recursive = T)
}

library(grDevices)
library(grid)
pdf(file=paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tQTL_coExQTL/', Gene1,'_','_', Gene2,'_', Gene1_name,'_',Gene2_name, '_', SNP,'_', treatment,"_RNAseq_coExQTL.pdf"), width = 12, height = 10, useDingbats = F)
print(p)
dev.off()

}

}


#------------------------------- visualization of isoforms inputs -------------------------------
library(rtracklayer)
GTF_TP <- rtracklayer::import('/Users/isarnassiri/Documents/RESULTS_USED/references/gencode.v38.annotation_revised.gtf')
GTF_TP <- as.data.frame(GTF_TP)
colnames(GTF_TP)

query='ENST00000337003'
GTF_TP_subset = GTF_TP[grep(query, GTF_TP$transcript_id),]
GTF_TP_subset = GTF_TP_subset[which(GTF_TP_subset$type %in% c('exon', 'transcript')),]

GTF_TP_subset = GTF_TP_subset[, match(c("seqnames", "gene_name", "type", "start", "end", "strand", "width", "gene_id", "transcript_id", "exon_number"), colnames(GTF_TP_subset))]

# GTF_TP_subset$gene_id = paste0('gene_id ', GTF_TP_subset$gene_id, ';')
# GTF_TP_subset$transcript_id = paste0('transcript_id ', GTF_TP_subset$transcript_id, ';')

rtracklayer::export(GTF_TP_subset, paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/',query ,'.gtf'))

#------------------------------- Enrichment analysis
Gene1='ENST00000529254'
Gene2='ENST00000468864'
Gene1_name='ARL2'
Gene2_name='TRAPPC10'
SNP='rs2243999'

REF=SNPs$REF[which(SNPs$SNP_ID == SNP)]
ALT=SNPs$ALT[which(SNPs$SNP_ID == SNP)]

i=1
treat = c('IFN', 'LPS24', 'UT')
print(treat[i])

tQTL_coExQTL = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/tQTL_coExQTL_subset_of_TP/', treat[i], '_coExQTL_tQTL_R2.txt'), stringsAsFactors = F, header = T)
tQTL_coExQTL = as.data.frame(tQTL_coExQTL)
tQTL_coExQTL = tQTL_coExQTL[-which(is.na(tQTL_coExQTL$pValDiff_adj_R2)), ]
# as I do not use permutation, so I need to perform multiple testing correction
library(qvalue)
tQTL_coExQTL$FDR <- qvalue(tQTL_coExQTL$pValDiff_adj_R2)$qvalues
tQTL_coExQTL = tQTL_coExQTL[which(tQTL_coExQTL$FDR < 0.01),]
tQTL_coExQTL_subset = tQTL_coExQTL[which(tQTL_coExQTL$Gene2 == Gene2 & tQTL_coExQTL$var_id == SNP),]
dim(tQTL_coExQTL_subset)


interaction = get('Query_interaction')
interaction_canonical = interaction[which(interaction$index %in% CS_DCEG[[1]]),]
dim(interaction_canonical)
interaction_canonical = interaction_canonical[which(interaction_canonical$gene1symbol %in% GTF$symbol),]
interaction_canonical = interaction_canonical[which(interaction_canonical$gene2symbol %in% GTF$symbol),]
interaction_canonical = interaction_canonical[order(interaction_canonical$pValDiff_adj, decreasing = F),]

#------- array profile [you can omit this part for dot plot]
array_interaction = fread(paste0('/t1-data/data/fairfaxlab/transcriptomics/RNAseq/eQTL/monocytes/analysis/QTL_analysis/co-expression_QTL/co-expression_array/All_DCA_',treatment,'.txt'),stringsAsFactors = F)
array_interaction$index = paste(array_interaction$Gene1,array_interaction$Gene2,sep='_')
interaction_canonical = interaction_canonical[which(interaction_canonical$index %in% array_interaction$index),]

selected_gene = 'ERAP2'
selected = interaction_canonical[which(interaction_canonical$gene2symbol == selected_gene),]
selected = unique(as.vector(c(selected$gene1symbol, selected_gene)))
length(selected)

#----- clusterProfiler
library(qusage)
gmtfile <- c("/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/h.all.v7.1.symbols.gmt") #c2.cp.reactome.v7.0.symbols.gmt
c5 <- read.gmt(gmtfile)

library(clusterProfiler)
egmt <- enricher(c(selected, selected_gene), TERM2GENE=c5)

expression_sub = expression[which(row.names(expression) %in% selected),]
dim(expression_sub)

library(cowplot)
setwd('/t1-data/data/fairfaxlab/transcriptomics/RNAseq/eQTL/monocytes/analysis/QTL_analysis/co-expression_QTL/BoxPlots/')
pdf('cnetplot_egmt.pdf', width = 8, height = 10, useDingbats = F)
cnetplot(egmt, foldChange=rowMeans(expression_sub), showCategory =4)
dev.off()

egmt@result$ID = gsub('HALLMARK_', '', egmt@result$ID)
egmt@result$Description = gsub('HALLMARK_', '', egmt@result$Description)

library(cowplot)
pdf('dotplot_egmt.pdf', width = 8, height = 10, useDingbats = F)
dotplot(egmt)
dev.off()



