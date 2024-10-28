library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
head(listAttributes(ensembl))
filters = listFilters(ensembl)
GTF <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(GTF) = c('symbol', 'id')

library(data.table)
SNPs = fread('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'


#--------------------------- query ---------------------------
setwd('/Users/isarnassiri/Documents/RESULTS_USED')
treat = c('IFN', 'LPS24', 'UT')

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
  
gQTL_coExQTL = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/gQTL_coExQTL/', treat[i], '_coExQTL_gQTL.txt'), stringsAsFactors = F, header = T)
gQTL_coExQTL = as.data.frame(gQTL_coExQTL)

# as I do not use permutation, so I need to perform multiple testing correction
library(qvalue)
gQTL_coExQTL$FDR <- qvalue(gQTL_coExQTL$pValDiff, lambda = 0)$qvalues

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
 
library(stringr)

Gene1_id=unique(GTF[which(GTF$symbol == Gene1),'id'])
Gene2_id=unique(GTF[which(GTF$symbol == Gene2),'id'])

TP_EXPESSION_Gene1 = expression[which(row.names(expression) %in% Gene1),]
dim(TP_EXPESSION_Gene1) 

TP_EXPESSION_Gene2 = expression[which(row.names(expression) %in% Gene2),]
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

if(!dir.exists(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gQTL_coExQTL/')))
{
  dir.create(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gQTL_coExQTL/'), recursive = T)
}

library(grDevices)
library(grid)
pdf(file=paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/gQTL_coExQTL/', Gene1,'_', Gene2,'_', SNP,'_', treatment,"_RNAseq_gcoExQTL.pdf"), width = 12, height = 10, useDingbats = F)
print(p)
dev.off()

}




#------------------------------- Enrichment analysis

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
gmtfile <- c("/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures/h.all.v7.1.symbols.gmt") #c2.cp.reactome.v7.0.symbols.gmt
c5 <- read.gmt(gmtfile)
#c5$ont = gsub('REACTOME_','',c5$ont)
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



