library(dplyr)
library( ggpubr)
library(data.table)

REMOVE_OUTLIERS = TRUE
PC = TRUE

SNPs = fread('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/MAF_imputed_Allsamples_revised.txt', stringsAsFactors = F)
colnames(SNPs)[3] = 'SNP_ID'

incorporate.PC = data.frame(IFN=6,LPS24=5,UT=7)

# destination folder of plots
folder = paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/cis-eQTL-Monocyte-Revisions/Figures_coExQTL/BoxPlots/tQTL_BoxPlots/')  

if(!dir.exists(folder))
{
  dir.create(folder, recursive = T)
}

#---------- single input

#-- nominal
ifn = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/tQTL/QTLtools_outputs-100kb-window/IFN/nominal_pass/tQTL_nominal_pass_1.txt', stringsAsFactors = F)
lps = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/tQTL/QTLtools_outputs-100kb-window/LPS24/nominal_pass/tQTL_nominal_pass_1.txt', stringsAsFactors = F)
ut = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/tQTL/QTLtools_outputs-100kb-window/UT/nominal_pass/tQTL_nominal_pass_1.txt', stringsAsFactors = F)

header_nominal = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/tQTL/header_nominal_tQTL.txt', stringsAsFactors = F)
colnames(ifn) = colnames(header_nominal)
colnames(lps) = colnames(header_nominal)
colnames(ut) = colnames(header_nominal)

library(qvalue)
ifn$FDR <- qvalue(ifn$nom_pval)$qvalues
lps$FDR <- qvalue(lps$nom_pval)$qvalues
ut$FDR <- qvalue(ut$nom_pval)$qvalues

treatment='LPS24' 

input1 = data.frame(SNP_ID='rs2760527',gene_name='RGS1',gene_id='ENST00000367459', stringsAsFactors = F)
input2 = data.frame(SNP_ID='rs807612',gene_name='DDX1',gene_id='ENST00000233084', stringsAsFactors = F)
input3 = data.frame(SNP_ID='rs93075',gene_name='SEPTIN9',gene_id='ENST00000423034', stringsAsFactors = F)
input4 = data.frame(SNP_ID='rs414340',gene_name='SEPTIN9',gene_id='ENST00000589250', stringsAsFactors = F)
input5 = data.frame(SNP_ID='rs2386849',gene_name='CTSC',gene_id='ENST00000227266', stringsAsFactors = F)
input6 = data.frame(SNP_ID='rs2439510',gene_name='SDC2',gene_id='ENST00000523877', stringsAsFactors = F)
input7 = data.frame(SNP_ID='rs2967172',gene_name='KIFC3',gene_id='ENST00000564204', stringsAsFactors = F)
input8 = data.frame(SNP_ID='rs2305789',gene_name='EIF3G',gene_id='ENST00000253108', stringsAsFactors = F)
input9 = data.frame(SNP_ID='rs835044',gene_name='NDUFA12',gene_id='ENST00000327772', stringsAsFactors = F)
input10 = data.frame(SNP_ID='rs835044',gene_name='NDUFA12',gene_id='ENST00000551991', stringsAsFactors = F)
input11 = data.frame(SNP_ID='rs2428486',gene_name='MICA',gene_id='ENST00000449934', stringsAsFactors = F)
input12 = data.frame(SNP_ID='rs2428486',gene_name='MICA',gene_id='ENST00000616296', stringsAsFactors = F)
input13 = data.frame(SNP_ID='rs9277533',gene_name='HLA-DPB1',gene_id='ENST00000471184', stringsAsFactors = F)
input14 = data.frame(SNP_ID='rs9277533',gene_name='HLA-DPA1',gene_id='ENST00000416804', stringsAsFactors = F)
input15 = data.frame(SNP_ID='rs4072037',gene_name='MUC1',gene_id='ENST00000612778', stringsAsFactors = F)
input16 = data.frame(SNP_ID='rs4072037',gene_name='MUC1',gene_id='ENST00000620103', stringsAsFactors = F)
input17 = data.frame(SNP_ID='rs10735079',gene_name='OAS1',gene_id='ENST00000202917', stringsAsFactors = F)
input18 = data.frame(SNP_ID='rs10735079',gene_name='OAS1',gene_id='ENST00000452357', stringsAsFactors = F)
input19 = data.frame(SNP_ID='rs2243999',gene_name='TRAPPC10',gene_id='ENST00000468864', stringsAsFactors = F)
input20 = data.frame(SNP_ID='rs61869825',gene_name='USMG5',gene_id='ENST00000337003', stringsAsFactors = F)
input21 = data.frame(SNP_ID='rs230519',gene_name='NFKB1',gene_id='ENST00000505458', stringsAsFactors = F)

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

for(i in 1:21)
{
  print(i)
  input = get(paste0('input', i))
  
  slop = round(query$slope[which(query$phe_id == input$gene_id & query$var_id == input$SNP_ID)], digits = 3) 
  
  print(input)
  print(paste0('Slop: ', slop))
  
  #visualization()
}

visualization()

visualization = function()
{
#----------
TP_id = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),select='id',stringsAsFactors = F) # includes IFB1
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

TP_GENOTYPE = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2_USED_for_BOXPLOT.txt'),skip = paste0(SNP, '_'), nrows = 1, stringsAsFactors = F)
TP_GENOTYPE$V1 = gsub('_', '', TP_GENOTYPE$V1)

GENOTYPE_sample_names = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),nrows = 1, stringsAsFactors = F)
colnames(TP_GENOTYPE) = colnames(GENOTYPE_sample_names)

#------------------------------------------- UT
treatment='UT'
library(data.table)
TP_EXPESSION = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),skip = TRANSCRIPT, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),nrows = 1,stringsAsFactors = F) # includes IFB1
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
expression.data = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),stringsAsFactors = F) # includes IFB1
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
TP_EXPESSION = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),skip = TRANSCRIPT, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),nrows = 1,stringsAsFactors = F) # includes IFB1
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
  expression.data = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),stringsAsFactors = F) # includes IFB1
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
TP_EXPESSION = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),skip = TRANSCRIPT, nrows = 1,stringsAsFactors = F) # includes IFB1
EXPESSION_sample_names = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),nrows = 1,stringsAsFactors = F) # includes IFB1
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
  expression.data = fread(paste0('/Users/isarnassiri/Documents//RESULTS_USED/Expression/Input_files_transcript/expression_',treatment,'.txt'),stringsAsFactors = F) # includes IFB1
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
