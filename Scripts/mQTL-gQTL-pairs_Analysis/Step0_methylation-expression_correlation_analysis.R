
library(MEAL)
library(brgedata)
library(MultiDataSet)
library(missMethyl)
library(minfi)
library(GenomicRanges)
library(ggplot2)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

treatment = 'UT'
##============================== expression matrix
library(data.table)
expmatrix <- fread(paste0('/Users/isarnassiri/Documents/RESULTS_USED/Expression/Input_files_gene/expression_',treatment,'.txt'), stringsAsFactors = F, header = T)
gene_name_id <- fread('/Users/isarnassiri/Documents/RESULTS_USED/QTLtools_eQTL_Monocyte_Moloc/gQTL/QTLtools_inputs/Biomart_DB_gene_name_id.txt', stringsAsFactors = F, header = T)
colnames(gene_name_id)[1] = c('symbol')
expmatrix = merge(gene_name_id, expmatrix, by = 'id')

library('org.Hs.eg.db')
de = as.character(mapIds(org.Hs.eg.db, expmatrix$symbol, 'ENTREZID', 'SYMBOL'))
length(de)
dim(expmatrix)

expmatrix = as.data.frame(expmatrix)
expmatrix = cbind(de , expmatrix)
dim(expmatrix)

row.names(expmatrix) = make.names(as.character(expmatrix$symbol),unique=T) 
dim(expmatrix)

expmatrix = expmatrix[,-2]
colnames(expmatrix)[1:2] = c('gene_id','symbol')
expmatrix = expmatrix[!is.na(expmatrix$gene_id),]

#-------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
geneshg19 = genes(txdb) #https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
colnames(geneshg19)
expmatrix_sub = merge(geneshg19, expmatrix, by = 'gene_id')
dim(expmatrix_sub) 

#------------- hg19
library(dplyr)
featureData = data.frame(symbol=expmatrix_sub$symbol, entrez=expmatrix_sub$gene_id , chromosome=as.character(expmatrix_sub$seqnames), start=expmatrix_sub$start, end=expmatrix_sub$end, strand=expmatrix_sub$strand)
dim(featureData)

phenoData = data.frame(SampleID= colnames(expmatrix_sub)[-c(1:7)], Disease = rep(1,length(colnames(expmatrix_sub)[-c(1:7)])), stringsAsFactors = F)
identical(colnames(expmatrix_sub)[-c(1:7)], phenoData$SampleID)

dim(phenoData)
dim(featureData)
dim(expmatrix_sub[,-c(1:7)])
assayData = as.matrix(expmatrix_sub[,-c(1:7)])
class(assayData)
class(phenoData)
phenoData = as.data.frame(phenoData)
featureData = as.data.frame(featureData)
row.names(phenoData) = colnames(expmatrix_sub[,-c(1:7)])
row.names(featureData) = make.names(as.character(expmatrix_sub$symbol),unique=T) 
row.names(assayData) = make.names(as.character(expmatrix_sub$symbol),unique=T) 
head(featureData)
colnames(phenoData)[1] = 'id'

library(data.table)
# write.table(featureData, '/Volumes/Fairfaxlab1/Laptop/Article_figures/Step63_transQTL_QTLtools/Exp_Met_analysis/featureData.txt', quote = F, row.names = T, sep = '\t' )
# write.table(phenoData, '/Volumes/Fairfaxlab1/Laptop/Article_figures/Step63_transQTL_QTLtools/Exp_Met_analysis/phenoData.txt', quote = F, row.names = T, sep = '\t' )
# write.table(assayData, '/Volumes/Fairfaxlab1/Laptop/Article_figures/Step63_transQTL_QTLtools/Exp_Met_analysis/assayData.txt', quote = F, row.names = T, sep = '\t' )

eset <- ExpressionSet(assayData = assayData,
                      phenoData = AnnotatedDataFrame(phenoData),
                      featureData = AnnotatedDataFrame(featureData))

##============================== methylation matrix
beta <- fread(paste0("/Users/isarnassiri/Documents/RESULTS_USED/methylationQTL/QTLtools_PC_Selection/Methylation_profile_",treatment,".bed.gz"), stringsAsFactors = F, header = T)
colnames(beta)
View(beta[1:10,1:10])
rMean = rowMeans(beta[,-c(1:6)])
summary(rMean)
length(which(rMean>=0.1))
selectedbeta = beta[which(rMean>=0.1),]
identical(colnames(expmatrix_sub[-c(1:7)]), colnames(selectedbeta)[-c(1:6)])
class(selectedbeta)
selectedbeta = as.data.frame(selectedbeta)
row.names(selectedbeta) = make.names(selectedbeta$pid,unique=T)
View(selectedbeta[1:10,])

betaextract_annotated = selectedbeta[,-c(1:6)]
class(betaextract_annotated[1,1])
# betaextract_annotated <- apply(betaextract_annotated, 2, function(x) as.numeric(as.character(x)))
betaextract_annotated = as.matrix(betaextract_annotated)
row.names(betaextract_annotated) = selectedbeta$pid
View(betaextract_annotated[1:10,])

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

betaextract_annotated = betaextract_annotated[,which(colnames(betaextract_annotated) %in% row.names(phenoData))]
phenoData = phenoData[which(row.names(phenoData) %in% colnames(betaextract_annotated)),]
identical(colnames(betaextract_annotated), row.names(phenoData))

betagrs <- makeGenomicRatioSetFromMatrix(betaextract_annotated, pData =  phenoData, array = "IlluminaHumanMethylation450k", mergeManifest = FALSE, what = c("Beta"), annotation = "ilmn12.hg19")

# write.table(betaextract_annotated, '/Volumes/Fairfaxlab1/Laptop/Article_figures/Step63_transQTL_QTLtools/Exp_Met_analysis/betaextract_annotated.txt', quote = F, row.names = T, sep = '\t' )
# write.table(cpg_featureData, '/Volumes/Fairfaxlab1/Laptop/Article_figures/Step63_transQTL_QTLtools/Exp_Met_analysis/cpg_featureData.txt', quote = F, row.names = T, sep = '\t' )
dim(beta)
dim(betaextract_annotated)

# View(beta[1:10,])
# View(assayData[1:10,])

##============================== target region

## ----New Multi Meth Exp-------------------------------------------------------
multi <- createMultiDataSet()
multi <- add_genexp(multi, eset)
multi <- add_methy(multi, betagrs)

## -----------------------------------------------------------------------------
# multi.filt <- multi[, , targetRange]

## ----Corr Meth Exp------------------------------------------------------------
methExprs <- correlationMethExprs(multi, verbose = T)
dim(methExprs)

#-- For figure 4 - correlation of methylation and expression
View(methExprs[1:5,1:5])
methExprs[which(methExprs$cpg == 'cg22687766' & methExprs$exprs == 'CD55'),]
# cpg exprs      Beta        se
# 810058 cg22687766  CD55 0.5483072 0.7677337
# P.Value adj.P.Val
# 810058 0.4761192 0.8553127
#--

table(methExprs$adj.P.Val<0.01)
methExprs_sub = methExprs[which(methExprs$adj.P.Val<0.01),]
head(methExprs_sub)
dim(methExprs_sub)
#  beta: coefficient of the methylation change
fwrite(methExprs, paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/eQTL_mQTL_pairs_analysis/',treatment,'correlationMethExprs_results.txt'), quote = F, row.names = F, sep = '\t')
getwd()

