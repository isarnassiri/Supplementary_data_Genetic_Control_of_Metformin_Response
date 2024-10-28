######################## case study ########################
# REF='G'
# ALT='A'
# 
# Gene2='ENST00000398979'
# Gene1='ENST00000515837'
# SNP='rs2860519'
# treatment='UT'

library(data.table)
treat = c('IFN', 'LPS24', 'UT')

i=3
j=2
for(i in 1:length(treat))
{
#-------------------
treatment=treat[i]
setwd('/Users/isarnassiri/Documents/RESULTS_USED/')

#--- gene symbol to Ensemble transcript IDs
library(rtracklayer)
GTF <- rtracklayer::import('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/tQTL/QTLtools_inputs/GRCh38.gtf')
GTF <- as.data.frame(GTF)

library(data.table)
library(dplyr)
expression = fread(paste0('Expression/Input_files_transcript/expression_',treatment,'_withtQTL_SAVER.txt'), stringsAsFactors = F)
expression = as.data.frame(expression)
row.names(expression) = expression$GeneID
expression = expression[,-1] 
colnames(expression) = gsub('X', '', colnames(expression))

expression = expression[which(row.names(expression) %in% GTF$transcript_id),]
GTF_subset = GTF[which(GTF$transcript_id %in% row.names(expression)),]
GTF_subset = GTF_subset[!duplicated(GTF_subset$transcript_id),]

dim(GTF_subset)
dim(expression)

GTF_subset = GTF_subset[match(row.names(expression), GTF_subset$transcript_id),]
identical(GTF_subset$transcript_id, row.names(expression))

# if you want to read an SNP you should make the file tab separated
TP_GENOTYPE = fread(paste0('Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),stringsAsFactors = F)
TP_GENOTYPE = as.data.frame(TP_GENOTYPE)
TP_GENOTYPE = TP_GENOTYPE[!duplicated(TP_GENOTYPE$id),]
row.names(TP_GENOTYPE) = TP_GENOTYPE$id
TP_GENOTYPE = TP_GENOTYPE[,-1] 

expression = expression[,which(colnames(expression) %in% colnames(TP_GENOTYPE))]
TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression))]

dim(expression)
dim(TP_GENOTYPE)

tQTL_coExQTL = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Input_files/coExQTL_results_rerun_using_conditional_peaks/tQTL_coExQTL_subset_of_TP/',treat[i],'_coExQTL_tQTL_R2.txt'), stringsAsFactors = F, header = T)
tQTL_coExQTL = as.data.frame(tQTL_coExQTL)

tQTL_coExQTL = tQTL_coExQTL[which(tQTL_coExQTL$Gene2_name %in% c('NT5C3B', 'LGALS2', 'ELP5')),]
tQTL_coExQTL = tQTL_coExQTL[which(tQTL_coExQTL$Gene2 %in% 'ENST00000574841'),]

table(is.na(tQTL_coExQTL$pValDiff_adj_R2))

#tQTL_coExQTL = tQTL_coExQTL[-which(is.na(tQTL_coExQTL$pValDiff_adj_R2)), ]

# as I do not use permutation I need to perform multiple testing correction
library(qvalue)
tQTL_coExQTL$FDR <- qvalue(tQTL_coExQTL$pValDiff_adj_R2, lambda = 0)$qvalues
range(tQTL_coExQTL$pValDiff_adj_R2)

#tQTL_coExQTL = tQTL_coExQTL[which(tQTL_coExQTL$FDR < 1e-5), ]
tQTL_coExQTL = tQTL_coExQTL[order(tQTL_coExQTL$FDR, decreasing = F), ]

print(dim(tQTL_coExQTL))

#-------------------
for(j in 1:dim(tQTL_coExQTL)[1])
{

REF=tQTL_coExQTL$REF[j]
ALT=tQTL_coExQTL$ALT[j]
Gene2=tQTL_coExQTL$Gene2[j]
Gene1=tQTL_coExQTL$Gene1[j]
SNP=tQTL_coExQTL$var_id[j]
Gene1_name=tQTL_coExQTL$Gene1_name[j]
Gene2_name=tQTL_coExQTL$Gene2_name[j]

library(stringr)
TP_EXPESSION_Gene1 = expression[which(row.names(expression) %in% Gene1),]
dim(TP_EXPESSION_Gene1) 

TP_EXPESSION_Gene2 = expression[which(row.names(expression) %in% Gene2),]
dim(TP_EXPESSION_Gene2) 

TP_GENOTYPE_sub = TP_GENOTYPE[which(row.names(TP_GENOTYPE) == SNP),]
dim(TP_GENOTYPE_sub)

TP_EXPESSION_Gene1 = as.data.frame(TP_EXPESSION_Gene1)
TP_EXPESSION_Gene2 = as.data.frame(TP_EXPESSION_Gene2)

TP_GENOTYPE_sub = as.data.frame(TP_GENOTYPE_sub)
TP_GENOTYPE_sub

TP_EXPESSION_Gene1 = as.data.frame(TP_EXPESSION_Gene1)
TP_EXPESSION_Gene2 = as.data.frame(TP_EXPESSION_Gene2)

#=============== visualization type 1
#--------------- Ref homozygote

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

if(length(tQTL_coExQTL$FDR[j]) == 1)
{
  dat_text <- data.frame(
      label = paste0('DC-FDR: ', scientific_10(tQTL_coExQTL$FDR[j], digits = 3) )
    )
}else{
  dat_text <- data.frame(
    label = paste0('DC-FDR > 0.05')
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

if(!dir.exists(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/tQTL_coExQTL_',treatment,'/')))
{
  dir.create(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/tQTL_coExQTL_',treatment,'/'), recursive = T)
}

library(grDevices)
library(grid)
pdf(file=paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/tQTL_coExQTL_',treatment,'/', Gene1,'_','_', Gene2,'_', Gene1_name,'_',Gene2_name, '_', SNP,'_', treatment,"_RNAseq_coExQTL.pdf"), width = 15, height = 10, useDingbats = F)
print(p)
dev.off()

}

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



