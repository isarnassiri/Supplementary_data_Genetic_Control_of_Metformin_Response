library(dplyr)
library( ggpubr)
library(data.table)

REMOVE_OUTLIERS = TRUE
PC = TRUE

SNPs = fread('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'

incorporate.PC = data.frame(LPS24=9,UT=11)

# destination folder of plots
folder = paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/metQTL_BoxPlots/')  

if(!dir.exists(folder))
{
  dir.create(folder, recursive = T)
}

#---------- single input

#-- nominal
lps = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_outputs/LPS24/nominal_pass/mQTL_nominal_pass_1.txt', stringsAsFactors = F)
ut = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_outputs/UT/nominal_pass/mQTL_nominal_pass_1.txt', stringsAsFactors = F)

header_nominal = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/QTLtools_outputs/gQTL_nominal_Header.txt', stringsAsFactors = F)
colnames(lps) = colnames(header_nominal)
colnames(ut) = colnames(header_nominal)
 
#--------------------------- query ---------------------------
treatment='UT' 
input1 = data.frame(SNP_ID='rs2914937', gene_id='cg22687766', gene_name = 'CD55', stringsAsFactors = F) 
input2 = data.frame(SNP_ID='rs6591507', gene_id='cg07745373', gene_name = 'DTX4', stringsAsFactors = F) 
input3 = data.frame(SNP_ID='rs61822619', gene_id='cg03948781', gene_name = 'DSTYK', stringsAsFactors = F) 
input4 = data.frame(SNP_ID='rs13331559', gene_id='cg04976154', gene_name = 'TELO2', stringsAsFactors = F) 

i=1
for(i in 1:1)
{
  print(i)
  input = get(paste0('input', i))
}

visualization()

visualization = function()
{
#----------
library(dplyr)
input = left_join(input, SNPs, by = 'SNP_ID')
t=1

# for (t in 1:dim(input)[1]) {
# 
#   tryCatch({

SNP = input$SNP_ID[t]
CPGSite = input$gene_id[t] 
GENE = input$gene_name[t] 

Allele_0 = paste0(input$REF[t], input$REF[t])
Allele_1 = paste0(input$REF[t], input$ALT[t])
Allele_2 = paste0(input$ALT[t], input$ALT[t])

#===== some of the SNP ids are subset of others and the fread retrive the first one (e.g. rs570631764 and rs570631764000); therefore I add _ to the SNP id to distinguish them.
TP_GENOTYPE = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2_USED_for_BOXPLOT.txt'),skip = paste0(SNP, '_'), nrows = 1, stringsAsFactors = F)
TP_GENOTYPE$V1 = gsub('_', '', TP_GENOTYPE$V1)

GENOTYPE_sample_names = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),nrows = 1, stringsAsFactors = F)
colnames(TP_GENOTYPE) = colnames(GENOTYPE_sample_names)

#------------------------------------------- UT
treatment='UT'
library(data.table)
TP_EXPESSION = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_inputs/expression_',treatment,'_Methylation.bed.gz'),skip = CPGSite, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_inputs/expression_',treatment,'_Methylation.bed.gz'),nrows = 1,stringsAsFactors = F) # includes IFB1
colnames(TP_EXPESSION) = colnames(EXPESSION_sample_names)

TP_EXPESSION = TP_EXPESSION[,which(colnames(TP_EXPESSION) %in% colnames(TP_GENOTYPE)),with=F]
TP_GENOTYPE_treatment = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(TP_EXPESSION)),with=F]
dim(TP_GENOTYPE_treatment)
dim(TP_EXPESSION)

identical(colnames(TP_EXPESSION) , colnames(TP_GENOTYPE_treatment))

dim(TP_EXPESSION)
dim(TP_GENOTYPE_treatment)

INPUT = rbind(as.matrix(TP_GENOTYPE_treatment), TP_EXPESSION)
INPUT = t(INPUT)
dim(INPUT)
colnames(INPUT) = c('genotype', paste0(GENE, '-', CPGSite) )
INPUT = as.data.frame(INPUT)

INPUT[,1] = as.numeric(as.character(INPUT[,1]))
table(INPUT[,1])
INPUT[,2] = as.numeric(as.character(INPUT[,2]))

#INPUT[,2] = INPUT[,2]+10
summary(INPUT[,2])
#------------- regress out

message("Calculating PC")
nPC = incorporate.PC[,(colnames(incorporate.PC) == treatment)]

if(t==1 & PC)
{
expression.data = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_inputs/expression_',treatment,'_Methylation.bed.gz'),stringsAsFactors = F) # includes IFB1
dim(expression.data)
expression.data = expression.data[,which(colnames(expression.data) %in% colnames(TP_GENOTYPE_treatment)),with=F]

identical(colnames(expression.data), colnames(TP_GENOTYPE_treatment))

variance<-apply(expression.data,1,var)

test.cov<-prcomp(t(expression.data[!variance==0,]),center=T,scale.=T)
test.cov<-data.frame(test.cov$x[1:length(expression.data[!variance==0,]),1:length(expression.data[!variance==0,])])
test.cov_UT<-test.cov[,1:nPC]
}

data <- data.frame(probe=INPUT[,2], genotype=INPUT[,1], test.cov_UT)

pcFit <- lm(data$probe ~ data$genotype + as.matrix(test.cov_UT))
summary(pcFit)
remodelled <- data$probe - rowSums(sapply(1:nPC, function(i) pcFit$coefficients[i+2]*test.cov_UT[,i]))
summary(remodelled)

# cor(remodelledd10, remodelledm10)
# summary(INPUT[,2]/10)
# summary(remodelled/10)

INPUT[,2] = remodelled
#-------------

INPUT$genotype[which(INPUT$genotype == 0)] = Allele_0
INPUT$genotype[which(INPUT$genotype == 1)] = Allele_1
INPUT$genotype[which(INPUT$genotype == 2)] = Allele_2

INPUT$genotype = factor(INPUT$genotype, levels = c(Allele_0,Allele_1,Allele_2))

x <- paste("INPUT_",treatment, sep="") # INPUT_UT
eval(call("<-", as.name(x), INPUT))

#------------------------------------------- LPS24
treatment='LPS24'
library(data.table)
TP_EXPESSION = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_inputs/expression_',treatment,'_Methylation.bed.gz'),skip = CPGSite, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_inputs/expression_',treatment,'_Methylation.bed.gz'),nrows = 1,stringsAsFactors = F) # includes IFB1
colnames(TP_EXPESSION) = colnames(EXPESSION_sample_names)

TP_EXPESSION = TP_EXPESSION[,which(colnames(TP_EXPESSION) %in% colnames(TP_GENOTYPE)),with=F]
TP_GENOTYPE_treatment = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(TP_EXPESSION)),with=F]
dim(TP_GENOTYPE_treatment)
dim(TP_EXPESSION)

identical(colnames(TP_EXPESSION) , colnames(TP_GENOTYPE_treatment))

dim(TP_EXPESSION)
dim(TP_GENOTYPE_treatment)

INPUT = rbind(as.matrix(TP_GENOTYPE_treatment), TP_EXPESSION)
INPUT = t(INPUT)
dim(INPUT)
colnames(INPUT) = c('genotype', paste0(GENE, '-', CPGSite) )
INPUT = as.data.frame(INPUT)

INPUT[,1] = as.numeric(as.character(INPUT[,1]))
table(INPUT[,1])
INPUT[,2] = as.numeric(as.character(INPUT[,2]))
summary((INPUT[,2]))

table(INPUT[,2] < 1)
#INPUT[,2] = INPUT[,2]+10
summary(INPUT[,2])
#------------- regress out

message("Calculating PC")
nPC = incorporate.PC[,(colnames(incorporate.PC) == treatment)]

if(t==1 & PC)
{
  expression.data = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_inputs/expression_',treatment,'_Methylation.bed.gz'),stringsAsFactors = F) # includes IFB1
  dim(expression.data)
  expression.data = expression.data[,which(colnames(expression.data) %in% colnames(TP_GENOTYPE_treatment)),with=F]
  
  identical(colnames(expression.data), colnames(TP_GENOTYPE_treatment))
  
  variance<-apply(expression.data,1,var)
  
  test.cov<-prcomp(t(expression.data[!variance==0,]),center=T,scale.=T)
  test.cov<-data.frame(test.cov$x[1:length(expression.data[!variance==0,]),1:length(expression.data[!variance==0,])])
  test.cov_LPS<-test.cov[,1:nPC]
}
data <- data.frame(probe=INPUT[,2], genotype=INPUT[,1], test.cov_LPS)

pcFit <- lm(data$probe ~ data$genotype + as.matrix(test.cov_LPS))
summary(pcFit)
remodelled <- data$probe - rowSums(sapply(1:nPC, function(i) pcFit$coefficients[i+2]*test.cov_LPS[,i]))
summary(INPUT[,2])
summary(remodelled)

INPUT[,2] = remodelled
summary(remodelled)
summary(INPUT[,2])
 
#-------------
INPUT$genotype[which(INPUT$genotype == 0)] = Allele_0
INPUT$genotype[which(INPUT$genotype == 1)] = Allele_1
INPUT$genotype[which(INPUT$genotype == 2)] = Allele_2

INPUT$genotype = factor(INPUT$genotype, levels = c(Allele_0,Allele_1,Allele_2))

x <- paste("INPUT_",treatment, sep="") # INPUT_UT
eval(call("<-", as.name(x), INPUT))

#----------------- visualization

#-----------
setwd(folder)
#-----------

#-----------
INPUT_LPS24$treatment = rep('LPS', dim(INPUT_LPS24)[1])
INPUT_UT$treatment = rep('Naive', dim(INPUT_UT)[1])

data<-rbind(INPUT_LPS24, INPUT_UT)
colnames(data)<-c("snp","gene.x","treatment") #,"stats"

data$treatment <- factor(data$treatment, levels = c('Naive','LPS'))
gene.name.1<-input$gene_id[t]
str(data)

set.seed(123)

library(ggplot2)
graph<-ggplot(data, aes((snp), gene.x)) + geom_jitter(colour="grey35", position=position_jitter(width=0.175), alpha=0.7, size=2.8) + geom_boxplot(colour="grey15",fill="grey85",alpha=0.45,outlier.size=0) + facet_wrap(~treatment, nrow=1) + theme_bw() + ylab(paste0(gene.name.1, ' Methylation')) + xlab(input$SNP_ID)+ ggtitle( paste0('mQTL: ', input$SNP_ID[t], ' & ', input$gene_name[t], ' (', input$gene_id[t], ')') )
graph <- graph + theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.title = element_text(size = 25 ), strip.text = element_text(size = 25, color = c("black"))) 

scientific_10 <- function(x, ...) {
  parse(text = gsub("e", "%*%10^", scales::label_scientific(...)(x)))
}

UT_FDR = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_outputs/UT_FDR_nominal_metQTL.txt', nrows=1, skip=which(ut$phe_id == input$gene_id & ut$var_id == input$SNP_ID))
LPS_FDR = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/mQTL/QTLtools_outputs/LPS24_FDR_nominal_metQTL.txt', nrows=1, skip=which(lps$phe_id == input$gene_id & lps$var_id == input$SNP_ID))

dat_text <- data.frame(
  label = c( paste0('P: ', scientific_10(UT_FDR$V3, digits = 3)), paste0('P: ', scientific_10(LPS_FDR$V3, digits = 3)) ),
  treatment = factor(c("Naive", "LPS"), levels = c('Naive','LPS'))
)

if(length(grep("charact", dat_text$label)) > 0)
{
  dat_text$label[grep("charact", dat_text$label)] = "P: 1 %*% 10^-1"
}

graph <- graph +
  geom_text(
    size    = 8,
    data    = dat_text,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.05,
    vjust   = 1.5,
    parse = TRUE 
  )

print(graph)

g <- ggplot_gtable(ggplot_build(graph))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c('#66CC99', '#FF9966')
k <- 1

for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

library(grid)
grid.draw(g)

pdf(paste(input$SNP_ID[t], input$gene_id[t], input$gene_name[t],'.pdf',sep = '_'),width = 8, height = 8, useDingbats = F)
print(grid.draw(g))
dev.off()

#   }, error=function(e){})
# }
getwd()  

}
