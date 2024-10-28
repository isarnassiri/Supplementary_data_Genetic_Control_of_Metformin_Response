
############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
library(dplyr)
library(ggpubr)
library(data.table)
library(biomaRt) 
library(ggplot2)
library(grid)

#----------------- Description -----------------
# Each box indicates a particular SNP allele, and the y-axis displays the expression levels of individuals with that allele. 
# The Fragments Per Kilobase Million (FPKM) metrix displays the expression levels of transcripts. 
# The 25th and 75th percentiles are represented in the box on the bottom and top. The box has a line that indicates the median of expression. 
# The whiskers cover both the minimum and maximum values, with the exception of outliers.

#----------------- Output -----------------
# An image file in pdf format

#----------------- Examples -----------------
treatment='LPS24'
input1 = data.frame(SNP_ID='rs2760527',gene_name='RGS1',gene_id='ENST00000367459', stringsAsFactors = F)
input = get('input1')
visualization()
## Note: You need to run the visualization() function first (Select lines 111-448 and hit the Run button/Enter).

##########################################################################################################################################################


#--- To allow comparison with the output of the regression model the optimal number of P.C. used to regress out expression changes attributable to the effect of the non-genetic covariates. 
PC = TRUE;
REMOVE_OUTLIERS = TRUE
incorporate.PC = data.frame(IFN=22,LPS24=29,UT=30)
#---

SNPs = fread('~/Genotype/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'

incorporate.PC = data.frame(IFN=6,LPS24=5,UT=7)

# destination folder of plots
folder = paste0('~/tQTL_BoxPlots/')  

if(!dir.exists(folder))
{
  dir.create(folder, recursive = T)
}

#---------- single input

#-- nominal
ifn = fread('~/IFN/nominal_pass/tQTL_nominal_pass_1.txt', stringsAsFactors = F)
lps = fread('~/LPS24/nominal_pass/tQTL_nominal_pass_1.txt', stringsAsFactors = F)
ut = fread('~/UT/nominal_pass/tQTL_nominal_pass_1.txt', stringsAsFactors = F)

header_nominal = fread('~/tQTL/header_nominal_tQTL.txt', stringsAsFactors = F)
colnames(ifn) = colnames(header_nominal)
colnames(lps) = colnames(header_nominal)
colnames(ut) = colnames(header_nominal)

library(qvalue)
ifn$FDR <- qvalue(ifn$nom_pval)$qvalues
lps$FDR <- qvalue(lps$nom_pval)$qvalues
ut$FDR <- qvalue(ut$nom_pval)$qvalues

#--------------------------- query ---------------------------
# treatment='LPS24' 
# input1 = data.frame(SNP_ID='rs2760527',gene_name='RGS1',gene_id='ENST00000367459', stringsAsFactors = F)
# input = get('input1')
# visualization()
## Note: You need to run the visualization() function first (Select lines 69-400 and hit the Run button/Enter).
#---------------------------

if(treatment == 'IFN')
{
  assign('query', ifn)
}

if(treatment == 'UT')
{
  assign('query', ut)
}

if(treatment == 'LPS24')
{
  assign('query', lps)
}

visualization = function()
{
#----------
TP_id = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),select='id',stringsAsFactors = F) # includes IFB1
input = input[which(input$gene_id %in% TP_id$id),]
dim(input)

#----------
input = input[order(input$gene_name),]
library(dplyr)
input = left_join(input, SNPs, by = 'SNP_ID')
t=1

# for (t in 1:dim(input)[1]) {
# 
#   tryCatch({

SNP = input$SNP_ID[t]
TRANSCRIPT = input$gene_id[t]
GENE = input$gene_name[t]

print(GENE)

Allele_0 = paste0(input$REF[t], input$REF[t])
Allele_1 = paste0(input$REF[t], input$ALT[t])
Allele_2 = paste0(input$ALT[t], input$ALT[t])

TP_GENOTYPE = fread(paste0('~/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2_USED_for_BOXPLOT.txt'),skip = paste0(SNP, '_'), nrows = 1, stringsAsFactors = F)
TP_GENOTYPE$V1 = gsub('_', '', TP_GENOTYPE$V1)

GENOTYPE_sample_names = fread(paste0('~/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),nrows = 1, stringsAsFactors = F)
colnames(TP_GENOTYPE) = colnames(GENOTYPE_sample_names)

#------------------------------------------- UT
treatment='UT'
library(data.table)
TP_EXPESSION = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),skip = TRANSCRIPT, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),nrows = 1,stringsAsFactors = F) # includes IFB1
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
colnames(INPUT) = c('genotype', paste0(GENE, '-', TRANSCRIPT) )
INPUT = as.data.frame(INPUT)

INPUT = INPUT[-1,]

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
expression.data = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),stringsAsFactors = F) # includes IFB1
dim(expression.data)
expression.data = expression.data[,which(colnames(expression.data) %in% colnames(TP_GENOTYPE_treatment)),with=F]

identical(colnames(expression.data), colnames(TP_GENOTYPE_treatment))

variance<-apply(expression.data[,-1],1,var)

test.cov<-prcomp(t(expression.data[!variance==0,-1]),center=T,scale.=T)
test.cov<-data.frame(test.cov$x[1:length(expression.data[!variance==0,-1]),1:length(expression.data[!variance==0,-1])])
test.cov_UT<-test.cov[,1:nPC]
}

data <- data.frame(probe=INPUT[,2], genotype=INPUT[,1], test.cov_UT)

pcFit <- lm(data$probe ~ data$genotype + as.matrix(test.cov_UT))
summary(pcFit)
remodelled <- data$probe - rowSums(sapply(1:nPC, function(i) pcFit$coefficients[i+2]*test.cov_UT[,i]))
summary(remodelled)

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
TP_EXPESSION = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),skip = TRANSCRIPT, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),nrows = 1,stringsAsFactors = F) # includes IFB1
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
colnames(INPUT) = c('genotype', paste0(GENE, '-', TRANSCRIPT) )
INPUT = as.data.frame(INPUT)

INPUT = INPUT[-1,]

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
  expression.data = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),stringsAsFactors = F) # includes IFB1
  dim(expression.data)
  expression.data = expression.data[,which(colnames(expression.data) %in% colnames(TP_GENOTYPE_treatment)),with=F]
  
  identical(colnames(expression.data), colnames(TP_GENOTYPE_treatment))
  
  variance<-apply(expression.data[,-1],1,var)
  
  test.cov<-prcomp(t(expression.data[!variance==0,-1]),center=T,scale.=T)
  test.cov<-data.frame(test.cov$x[1:length(expression.data[!variance==0,-1]),1:length(expression.data[!variance==0,-1])])
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

#------------------------------------------- IFNg
treatment='IFN'
library(data.table)
TP_EXPESSION = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),skip = TRANSCRIPT, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),nrows = 1,stringsAsFactors = F) # includes IFB1
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
colnames(INPUT) = c('genotype', paste0(GENE, '-', TRANSCRIPT) )
INPUT = as.data.frame(INPUT)

INPUT = INPUT[-1,]

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
  expression.data = fread(paste0('~/Expression/Input_files_transcript/expression_',treatment,'.txt'),stringsAsFactors = F) # includes IFB1
  dim(expression.data)
  expression.data = expression.data[,which(colnames(expression.data) %in% colnames(TP_GENOTYPE_treatment)),with=F]
  
  identical(colnames(expression.data), colnames(TP_GENOTYPE_treatment))
  
  variance<-apply(expression.data[,-1],1,var)
  
  test.cov<-prcomp(t(expression.data[!variance==0,-1]),center=T,scale.=T)
  test.cov<-data.frame(test.cov$x[1:length(expression.data[!variance==0,-1]),1:length(expression.data[!variance==0,-1])])
  test.cov_IFN<-test.cov[,1:nPC]
}
data <- data.frame(probe=INPUT[,2], genotype=INPUT[,1], test.cov_IFN)

pcFit <- lm(data$probe ~ data$genotype + as.matrix(test.cov_IFN))
summary(pcFit)
remodelled <- data$probe - rowSums(sapply(1:nPC, function(i) pcFit$coefficients[i+2]*test.cov_IFN[,i]))
summary(INPUT[,2])
summary(remodelled)

INPUT[,2] = remodelled

#-------------

INPUT$genotype[which(INPUT$genotype == 0)] = Allele_0
INPUT$genotype[which(INPUT$genotype == 1)] = Allele_1
INPUT$genotype[which(INPUT$genotype == 2)] = Allele_2

INPUT$genotype = factor(INPUT$genotype, levels = c(Allele_0,Allele_1,Allele_2))

x <- paste("INPUT_",treatment, sep="") # INPUT_UT
eval(call("<-", as.name(x), INPUT))

# fwrite(INPUT, '/Users/isar.nassiri/Desktop/OAS1_Isoforms_expression_genotypes_IFNg_with_regressout.txt', quote = F, sep = '\t', row.names = T)

#----------------- visualization

#-----------
setwd(folder)
#-----------

#-----------
INPUT_LPS24$treatment = rep('LPS', dim(INPUT_LPS24)[1])
INPUT_IFN$treatment = rep('IFNg', dim(INPUT_IFN)[1])
INPUT_UT$treatment = rep('Naive', dim(INPUT_UT)[1])

if(REMOVE_OUTLIERS)
{

#--- remove outliers
Q <- quantile(as.vector(t(INPUT_LPS24[,2])), probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(as.vector(t(INPUT_LPS24[,2])))
up <-  as.double(Q[2]+1.5*iqr) # Upper Range
low <- as.double(Q[1]-1.5*iqr)  # Lower Range
INPUT_LPS24_subset = INPUT_LPS24[which(INPUT_LPS24[,2] < up &  INPUT_LPS24[,2]  > low ),]

Q <- quantile(as.vector(t(INPUT_UT[,2])), probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(as.vector(t(INPUT_UT[,2])))
up <-  as.double(Q[2]+1.5*iqr) # Upper Range
low <- as.double(Q[1]-1.5*iqr)  # Lower Range
INPUT_UT_subset = INPUT_UT[which(INPUT_UT[,2] < up &  INPUT_UT[,2]  > low ),]

Q <- quantile(as.vector(t(INPUT_IFN[,2])), probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(as.vector(t(INPUT_IFN[,2])))
up <-  as.double(Q[2]+1.5*iqr) # Upper Range
low <- as.double(Q[1]-1.5*iqr)  # Lower Range
INPUT_IFN_subset = INPUT_IFN[which(INPUT_IFN[,2] < up &  INPUT_IFN[,2]  > low ),]
#---

data<-rbind(INPUT_LPS24_subset, INPUT_UT_subset, INPUT_IFN_subset)
}else{
  data<-rbind(INPUT_LPS24, INPUT_UT, INPUT_IFN)
}

colnames(data)<-c("snp","gene.x","treatment") #,"stats"

data$treatment <- factor(data$treatment, levels = c('Naive','LPS','IFNg'))
gene.name.1<-input$gene_name[t]
str(data)

set.seed(123)

library(ggplot2)
graph<-ggplot(data, aes((snp), gene.x)) + geom_jitter(colour="grey35", position=position_jitter(width=0.175), alpha=0.7, size=2.8) + geom_boxplot(colour="grey15",fill="grey85",alpha=0.45,outlier.size=0) + facet_wrap(~treatment, nrow=1) + theme_bw() + ylab(paste0(gene.name.1, ' Expression')) + xlab(input$SNP_ID)+ ggtitle( paste0('tQTL: ', input$SNP_ID[t], ' & ', input$gene_name[t], ' (', input$gene_id[t], ')') )
graph <- graph + theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.title = element_text(size = 25 ), strip.text = element_text(size = 25, color = c("black"))) 

scientific_10 <- function(x, ...) {
  parse(text = gsub("e", "%*%10^", scales::label_scientific(...)(x)))
}

dat_text <- data.frame(
  label = c( paste0('P: ', scientific_10(ut[which(ut$phe_id == input$gene_id & ut$var_id == input$SNP_ID),'FDR'], digits = 3)), paste0('P: ', scientific_10(lps[which(lps$phe_id == input$gene_id & lps$var_id == input$SNP_ID),'FDR'], digits = 3)) , paste0('P: ', scientific_10(ifn[which(ifn$phe_id == input$gene_id & ifn$var_id == input$SNP_ID),'FDR'], digits = 3)) ),
  treatment = factor(c("Naive", "LPS", "IFNg"), levels = c('Naive','LPS','IFNg'))
)

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
fills <- c('#66CC99', '#FF9966', '#9999CC')
k <- 1

for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

library(grid)
grid.draw(g)

pdf(paste(input$SNP_ID[t], input$gene_id[t], input$gene_name[t],'.pdf',sep = '_'),width = 10, height = 8, useDingbats = F)
print(grid.draw(g))
dev.off()

#   }, error=function(e){})
# }
getwd()  

}
