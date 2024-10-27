# module purge
# module add R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
# module add R-bundle-CRAN/2023.12-foss-2023a

.libPaths("/well/parkkinen/users/gbf362/R/4.1/skylake/")

library(data.table)
treat=c('LPS24','UT')
j=2
treatment=treat[j]

dir.create('/exafs1/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/tmQTL_coExQTL/', recursive=T)
setwd('/exafs1/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/tmQTL_coExQTL/')

#---genesymboltoEnsembletranscriptIDs
# library(rtracklayer)
# GTF<-rtracklayer::import('/exafs1/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/GENOMEANNOT/GRCh38.gtf')
# GTF<-as.data.frame(GTF)

setwd('/exafs1/well/parkkinen/users/gbf362/Analysis_data/ciseQTL_Monocyte/')

library(data.table)
library(dplyr)
expression=fread(paste0('Input_files_transcript/expression_',treatment,'_SAVER.txt'),stringsAsFactors=F)
expression=as.data.frame(expression)
row.names(expression)=expression$GeneID
expression=expression[,-1]
colnames(expression)=gsub('X','',colnames(expression))

# expression=expression[which(row.names(expression)%in%GTF$transcript_id),]
# GTF_subset=GTF[which(GTF$transcript_id%in%row.names(expression)),]
# GTF_subset=GTF_subset[!duplicated(GTF_subset$transcript_id),]
# 
# dim(GTF_subset)
# dim(expression)
# 
# GTF_subset=GTF_subset[match(row.names(expression),GTF_subset$transcript_id),]
# identical(GTF_subset$transcript_id,row.names(expression))

#ifyouwanttoreadanSNPyoushouldmakethefiletabseparated
TP_GENOTYPE=fread(paste0('Input_files_Genotype/Monocyte_imputed_matrixQTL_Allsamples_justSNPs_format2.txt'),stringsAsFactors=F)
TP_GENOTYPE=as.data.frame(TP_GENOTYPE)
TP_GENOTYPE=TP_GENOTYPE[!duplicated(TP_GENOTYPE$id),]
row.names(TP_GENOTYPE)=TP_GENOTYPE$id
TP_GENOTYPE=TP_GENOTYPE[,-1]

expression=expression[,which(colnames(expression)%in%colnames(TP_GENOTYPE))]
TP_GENOTYPE=TP_GENOTYPE[,which(colnames(TP_GENOTYPE)%in%colnames(expression))]

dim(expression)
dim(TP_GENOTYPE)

METHYLATION = fread(paste0('Input_files_methylation/Methylation_',treatment,'_hg38_QTLtools.bed.gz'), stringsAsFactors = F, header = T)
METHYLATION = as.data.frame(METHYLATION)
METHYLATION = METHYLATION[,-c(1:3)]
colnames(METHYLATION)[1] = 'id'
colnames(METHYLATION)[-1] = sub('X','',colnames(METHYLATION)[-1])
row.names(METHYLATION) = METHYLATION[,1]
METHYLATION = METHYLATION[,which(colnames(METHYLATION) %in% colnames(expression))]
expression = expression[,which(colnames(expression) %in%  colnames(METHYLATION))]

identical(colnames(METHYLATION), colnames(expression))
expression_methylation = rbind(METHYLATION, expression)

expression_methylation = expression_methylation[,which(colnames(expression_methylation) %in% colnames(TP_GENOTYPE))]
TP_GENOTYPE = TP_GENOTYPE[,which(colnames(TP_GENOTYPE) %in% colnames(expression_methylation))]
identical(colnames(expression_methylation), colnames(TP_GENOTYPE))
dim(expression_methylation)

METHYLATION_sites_details = fread(paste0('Input_files_methylation/Methylation_',treatment,'_hg38_QTLtools.bed.gz'), stringsAsFactors = F, header = T , select=c(1:4,6))
METHYLATION_sites_details = as.data.frame(METHYLATION_sites_details)
colnames(METHYLATION_sites_details)[which(colnames(METHYLATION_sites_details) == "pid")] = 'Gene1'
colnames(METHYLATION_sites_details)[-which(colnames(METHYLATION_sites_details) == "Gene1")] = paste0(colnames(METHYLATION_sites_details)[-which(colnames(METHYLATION_sites_details) == "Gene1")], '_cpg')
#--------
print(treat[j])
tQTL_coExQTL=fread(paste0('coExQTL/tQTL_coExQTL_subset_of_TP/', treat[j], '_coExQTL_tQTL_R2.txt'), stringsAsFactors=F, header=T)
tQTL_coExQTL=as.data.frame(tQTL_coExQTL)
tQTL_coExQTL=tQTL_coExQTL[which(tQTL_coExQTL$pValDiff_adj<=1e-5),]
tQTL_coExQTL=tQTL_coExQTL[order(tQTL_coExQTL$pValDiff_adj,decreasing=F),]
tQTL_coExQTL=tQTL_coExQTL[!duplicated(tQTL_coExQTL$Gene2),]
dim(tQTL_coExQTL)
#--

library(DGCA)
library(ggpubr)

i=1
for(i in 1:dim(tQTL_coExQTL)[1]){
  tryCatch({
    
    print(i)
    
    coExQTL_temp = tQTL_coExQTL[i,]
    Expression_subset = expression[which(row.names(expression) %in% coExQTL_temp$Gene2 ),]
    TP_GENOTYPE_sub = TP_GENOTYPE[which(row.names(TP_GENOTYPE) == coExQTL_temp$var_id),]
    
    Methylation = expression_methylation[-grep('ENST', row.names(expression_methylation)),]
    dim(Methylation)
    
    Expression_subset = rbind(Expression_subset, Methylation)
    identical(colnames(Expression_subset), colnames(TP_GENOTYPE_sub))
    
    #--- remove outliers - Expression
    Q <- quantile(as.vector(t(Expression_subset[1,])), probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR( as.vector(t(Expression_subset[1,])) )
    up <-  Q[2]+1.5*iqr # Upper Range
    low<- Q[1]-1.5*iqr  # Lower Range
    
    if(length(which(as.vector(t(Expression_subset[1,])) < up & as.vector(t(Expression_subset[1,])) > low)) != 0)
    {
      Expression_subset = Expression_subset[,which(as.vector(t(Expression_subset[1,])) < up & as.vector(t(Expression_subset[1,])) > low)]
    }
    
    dim(Expression_subset)
    #---
    
    #----------Refhomozygote
    expr0=Expression_subset[,which(colnames(Expression_subset)%in%colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,] == 0)])]
    expr1=Expression_subset[,which(colnames(Expression_subset)%in%colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,] == 1)])]
    expr2=Expression_subset[,which(colnames(Expression_subset)%in%colnames(TP_GENOTYPE_sub)[which(TP_GENOTYPE_sub[1,] == 2)])]
    #------------------------------------
    
    set.seed(1234)
    
    design_matrix<-cbind(rep(c(1,0),c(dim(expr2)[2],dim(expr0)[2])),rep(c(0,1),c(dim(expr2)[2],dim(expr0)[2])))
    colnames(design_matrix)<-c("G2","G0")
    INPUT=cbind(expr2,expr0)
    ddcor_G2_G0=ddcorAll(inputMat=INPUT,design=design_matrix,compare=c("G2","G0"),adjust="fdr",nPerm=0,corrType="pearson",splitSet=coExQTL_temp$Gene2,sigThresh=0.05,sortBy="pValDiff_adj",verbose=TRUE)

    design_matrix<-cbind(rep(c(1,0),c(dim(expr2)[2],dim(expr1)[2])),rep(c(0,1),c(dim(expr2)[2],dim(expr1)[2])))
    colnames(design_matrix)<-c("G2","G1")
    INPUT=cbind(expr2,expr1)
    ddcor_G2_G1=ddcorAll(inputMat=INPUT, design=design_matrix, compare=c("G2","G1"),adjust="fdr", nPerm=0, corrType="pearson", splitSet=coExQTL_temp$Gene2, sigThresh=0.05, sortBy="pValDiff_adj", verbose=TRUE)

    ddcor_G2_G0 = cbind(ddcor_G2_G0, coExQTL_temp[,c('Gene2_name', 'var_id', 'REF', 'ALT')])
    ddcor_G2_G1 = cbind(ddcor_G2_G1, coExQTL_temp[,c('Gene2_name', 'var_id', 'REF', 'ALT')])
    
    colnames(ddcor_G2_G0) = c('Gene1','Gene2','GM_cor','GM_pVal','GR_cor','GR_pVal','zScoreDiff','pValDiff','pValDiff_adj','Classes','index','var_id','REF','ALT')
    colnames(ddcor_G2_G1) = c('Gene1','Gene2','GM_cor','GM_pVal','GR_cor','GR_pVal','zScoreDiff','pValDiff','pValDiff_adj','Classes','index','var_id','REF','ALT')
    
    ddcor_G2_G0$index = 'G2_G0'
    ddcor_G2_G1$index = 'G2_G1'
    
    ddcor_G2_G0$Gene2_name = coExQTL_temp$Gene2_name
    ddcor_G2_G1$Gene2_name = coExQTL_temp$Gene2_name
    
    ddcor_G2_G1 = ddcor_G2_G1[which(ddcor_G2_G1$Classes != 'NonSig'),]
    ddcor_G2_G0 = ddcor_G2_G0[which(ddcor_G2_G0$Classes != 'NonSig'),]
    
    dim(ddcor_G2_G1)
    dim(ddcor_G2_G0)
    
    #--- add cpg details
    
    if(dim(ddcor_G2_G0)[1]>0)
    {
      ddcor_G2_G0 = merge(ddcor_G2_G0, METHYLATION_sites_details, by = 'Gene1')
    }
    
    if(dim(ddcor_G2_G1)[1]>0)
    {
      ddcor_G2_G1 = merge(ddcor_G2_G1, METHYLATION_sites_details, by = 'Gene1')
    }
    
    #---
    
    if(i==1){RESULT = rbind(ddcor_G2_G0, ddcor_G2_G1)}else{RESULT = rbind(RESULT, rbind(ddcor_G2_G0, ddcor_G2_G1))}
    
    print(i)
    print(dim(RESULT))
  
  },error=function(e){})
}

fwrite(RESULT, paste0('coExQTL/tmQTL_coExQTL/',treatment,'_coExQTL_tmQTL.txt'),quote=F,row.names=F,sep='\t')

# Note: The correction for multiple tests is not effective because I only have one test.
