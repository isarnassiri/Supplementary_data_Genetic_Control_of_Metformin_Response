
############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
library(dplyr)
library(data.table)
library(grDevices)
library(grid)
library(ggpubr)
library(ggplot2)
library(stringr)
library(qvalue)
library(biomaRt)

#----------------- Description -----------------
# The expression of a gene pair is coloured by the genotype of the SNP.
# Each data point on the horizontal and vertical axis indicates values for a single sample. Regression lines are shown for categories of genotypes. 

#----------------- Output -----------------
# An image file in pdf format

#----------------- Examples -----------------
treat = c('IFN')
input1 = data.frame(var_id='rs7305461', Gene2='RPS26', Gene1='KPNA2', stringsAsFactors = F) 
input = get(input1)
visualization()

#- Note: You need to call the visualization() function first (Select lines 47-147 and hit the Run button/Enter).

##########################################################################################################################################################

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
GTF <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
colnames(GTF) = c('symbol', 'id')

SNPs = fread('~/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'

#--------------------------- query ---------------------------
treat = c('IFN')
input1 = data.frame(var_id='rs7305461', Gene2='RPS26', Gene1='KPNA2', stringsAsFactors = F) 
input = get(input1)
visualization()
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
gQTL_coExQTL$FDR <- qvalue(gQTL_coExQTL$pValDiff, lambda = 0)$qvalues

expression = fread(paste0('Expression/Input_files_gene/expression_',treatment,'.txt'), stringsAsFactors = F)
expression = as.data.frame(expression)
row.names(expression) = expression$id
expression = expression[,-1]

PROBEID = data.frame(id = row.names(expression), stringsAsFactors = F)
SYMBOL = right_join(GTF, PROBEID, by ='id')
SYMBOL = SYMBOL[-which(SYMBOL$symbol ==''),]
SYMBOL = SYMBOL[!is.na(SYMBOL$symbol),]
SYMBOL = SYMBOL[!duplicated(SYMBOL$symbol),]
SYMBOL = SYMBOL[!duplicated(SYMBOL$id),]

expression = expression[which(row.names(expression) %in% SYMBOL$id),]
SYMBOL = SYMBOL[which(SYMBOL$id %in% row.names(expression)),]

expression = expression[!duplicated(row.names(expression)),]
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

Gene1_id=unique(GTF[which(GTF$symbol == Gene1),'id'])
Gene2_id=unique(GTF[which(GTF$symbol == Gene2),'id'])

TP_EXPESSION_Gene1 = expression[which(row.names(expression) %in% Gene1),]
TP_EXPESSION_Gene2 = expression[which(row.names(expression) %in% Gene2),]
TP_GENOTYPE_sub = TP_GENOTYPE#[which(row.names(TP_GENOTYPE) == SNP),]

TP_EXPESSION_Gene1 = as.data.frame(TP_EXPESSION_Gene1)
TP_EXPESSION_Gene2 = as.data.frame(TP_EXPESSION_Gene2)
TP_GENOTYPE_sub = as.data.frame(TP_GENOTYPE_sub)

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
theme_set(theme_bw())

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

if(!dir.exists(paste0('~/BoxPlots/gQTL_coExQTL/')))
{
  dir.create(paste0('~/BoxPlots/gQTL_coExQTL/'), recursive = T)
}

pdf(file=paste0('~/BoxPlots/gQTL_coExQTL/', Gene1,'_', Gene2,'_', SNP,'_', treatment,"_RNAseq_gcoExQTL.pdf"), width = 12, height = 10, useDingbats = F)
print(p)
dev.off()

}
