# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("MRCIEU/MRInstruments")
# remotes::install_github("izhbannikov/haplor")
#------------------------ backend data
library(TwoSampleMR)
library(MRInstruments)
library(haploR)

# List available GWASs
ao <- available_outcomes()
ao = ao[!is.na(ao$ncase),]
dim(ao)

setwd('D:/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/')
write.table(ao, 'List_of_10141_GWASsummaryStat.txt', quote = F, row.names = F, sep = '\t')

# List of available DBs
DBs = as.data.frame(table(gsub('-.*','',ao$id)))

UKBIO = ao[grep('ukb', ao$id),]
UKBIOtrait = as.data.frame(table(UKBIO$trait))
UKBIO=UKBIO[grep('covid|COVID|Cancer|cancer|immun|Immun|inflammat|Inflammat', UKBIO$trait),]
UKBIO=UKBIO[-grep('non|Non|Illnesses of', UKBIO$trait),]
UKBIO <- UKBIO[grepl("European", UKBIO$population), ]
dim(UKBIO)

BBJ = ao[grep('bbj', ao$id),]
BBJtrait = as.data.frame(table(BBJ$trait))
BBJ <- BBJ[grepl("European", BBJ$population), ]   #mainly East Asian
dim(BBJ)

EBI = ao[grep('ebi', ao$id),]
EBItrait = as.data.frame(table(EBI$trait))
EBI <- EBI[grepl("European", EBI$population), ]
dim(EBI)

ieu = ao[grep('ieu', ao$id),]
ieutrait = as.data.frame(table(ieu$trait))
ieu <- ieu[grepl("European", ieu$population), ]
dim(ieu)

#The terms seems redundant and ambiguous
# FINN = ao[grep('finn', ao$id),]
# FINNtrait = as.data.frame(table(FINN$trait))
# FINN=FINN[grep('covid|COVID|Cancer|cancer|immun|Immun|inflammat|Inflammat', FINN$trait),]
# FINN <- FINN[grepl("European", FINN$population), ]
# dim(FINN)

table(ao$subcategory)  # Includes Psychiatric / neurological
toMatch <- c("Autoimmune / inflammatory",
             'Cofactors and vitamins',
             "Haemotological", "Cancer",
             "Immune system", "Immune cell subset frequency", "Metabolites ratio", "Nucleotide","Paediatric disease")
subcategories <- ao[grepl(paste(toMatch,collapse="|"), ao$subcategory),]
subcategories <- subcategories[grepl("European", subcategories$population), ]
dim(subcategories)

allDBs = do.call("rbind", list(subcategories, UKBIO, BBJ, EBI, ieu))
table(duplicated(allDBs$id))
allDBs = allDBs[!duplicated(allDBs$id),]

allDBs = allDBs[-grep('Illnesses of|self-reported|Self-reported|Non-cancer|non-cancer|noninflammatory|Noninflammatory', UKBIO$trait),]
dim(allDBs)

table(allDBs$consortium)

#--------------------------- query
treatment = 'LPS24'
if(treatment == 'IFN'){query_nCases = 139} else if (treatment == 'LPS24') {query_nCases = 176 } else {query_nCases = 176 }

library(data.table)
library(coloc) 
library(qvalue)
#============================== read gQTLs

SNPs = fread('D:/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/QTLtools_inputs/MAF_imputed_Allsamples_TYPED.txt', stringsAsFactors = F, header = T, fill=TRUE)
SNPs = as.data.frame(SNPs)
SNPs = SNPs[,c('ID', 'REF', 'ALT', '0', '1', '2', 'mac', 'maf', 'TYPED')]
colnames(SNPs)[which(colnames(SNPs) %in% 'ID')] = 'var_id'
SNPs = SNPs[order(SNPs$maf, decreasing = T),]
SNPs = SNPs[!duplicated(SNPs$var_id),]
table(duplicated(SNPs$var_id))
head(SNPs)

header = fread('D:/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/gQTL_nominal_Header.txt', stringsAsFactors = F, header = F)
headerC = fread('D:/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/gQTL_conditional_Header.txt', stringsAsFactors = F, header = F)
headerP = fread('D:/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/gQTL_permutation_Header.txt', stringsAsFactors = F, header = F)

setwd('D:/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/')

Transcript_conditional_LPS24 = fread('QTLtools_outputs/LPS24/conditional_pass/gQTL_conditional_pass_1.txt', stringsAsFactors = F)
colnames(Transcript_conditional_LPS24) = as.character(headerC)
Transcript_conditional_LPS24 = Transcript_conditional_LPS24[-which(Transcript_conditional_LPS24$var_id == '.'),]

# bwd_pval: The nominal backward p-value of the association between the most significant variant and the phenotype.
Transcript_conditional_LPS24$FDR <- qvalue(p = Transcript_conditional_LPS24$bwd_pval)$qvalues
Transcript_conditional_LPS24_sub = Transcript_conditional_LPS24[Transcript_conditional_LPS24$FDR < 0.01,]
length(unique(Transcript_conditional_LPS24_sub$phe_id))

table(Transcript_conditional_LPS24_sub$bwd_best_hit)

Transcript_conditional_UT = fread('QTLtools_outputs/UT/conditional_pass/gQTL_conditional_pass_1.txt', stringsAsFactors = F)
colnames(Transcript_conditional_UT) = as.character(headerC)
Transcript_conditional_UT = Transcript_conditional_UT[-which(Transcript_conditional_UT$var_id == '.'),]

# bwd_pval: The nominal backward p-value of the association between the most significant variant and the phenotype.
Transcript_conditional_UT$FDR <- qvalue(p = Transcript_conditional_UT$bwd_pval)$qvalues
Transcript_conditional_UT_sub = Transcript_conditional_UT[Transcript_conditional_UT$FDR < 0.01,]
length(unique(Transcript_conditional_UT_sub$phe_id))

(table(Transcript_conditional_UT_sub$bwd_best_hit))

Transcript_conditional_IFN = fread('QTLtools_outputs/IFN/conditional_pass/gQTL_conditional_pass_1.txt', stringsAsFactors = F)
colnames(Transcript_conditional_IFN) = as.character(headerC)
Transcript_conditional_IFN = Transcript_conditional_IFN[-which(Transcript_conditional_IFN$var_id == '.'),]

# bwd_pval: The nominal backward p-value of the association between the most significant variant and the phenotype.
Transcript_conditional_IFN$FDR <- qvalue(p = Transcript_conditional_IFN$bwd_pval)$qvalues
Transcript_conditional_IFN_sub = Transcript_conditional_IFN[Transcript_conditional_IFN$FDR < 0.01,]
length(unique(Transcript_conditional_IFN_sub$phe_id))

Transcript_conditional_LPS24_best_hit = Transcript_conditional_LPS24_sub[which(Transcript_conditional_LPS24_sub$bwd_best_hit == 1),]
Transcript_conditional_LPS24_best_hit$gQTLs = paste(Transcript_conditional_LPS24_best_hit$phe_id, Transcript_conditional_LPS24_best_hit$var_id, sep = '_')
Transcript_conditional_UT_best_hit = Transcript_conditional_UT_sub[which(Transcript_conditional_UT_sub$bwd_best_hit == 1),]
Transcript_conditional_UT_best_hit$gQTLs = paste(Transcript_conditional_UT_best_hit$phe_id, Transcript_conditional_UT_best_hit$var_id, sep = '_')
Transcript_conditional_IFN_best_hit = Transcript_conditional_IFN_sub[which(Transcript_conditional_IFN_sub$bwd_best_hit == 1),]
Transcript_conditional_IFN_best_hit$gQTLs = paste(Transcript_conditional_IFN_best_hit$phe_id, Transcript_conditional_IFN_best_hit$var_id, sep = '_')

Transcript_conditional_LPS24_without_best_hit = Transcript_conditional_LPS24[-which(Transcript_conditional_LPS24$phe_id %in% Transcript_conditional_LPS24_best_hit$phe_id),]
Transcript_conditional_LPS24_without_best_hit = Transcript_conditional_LPS24_without_best_hit[order(Transcript_conditional_LPS24_without_best_hit$FDR, decreasing = F),]
Transcript_conditional_LPS24_without_best_hit = Transcript_conditional_LPS24_without_best_hit[!duplicated(Transcript_conditional_LPS24_without_best_hit$phe_id),]
Transcript_conditional_LPS24_without_best_hit$gQTLs = paste(Transcript_conditional_LPS24_without_best_hit$phe_id, Transcript_conditional_LPS24_without_best_hit$var_id, sep = '_')

Transcript_conditional_UT_without_best_hit = Transcript_conditional_UT[-which(Transcript_conditional_UT$phe_id %in% Transcript_conditional_UT_best_hit$phe_id),]
Transcript_conditional_UT_without_best_hit = Transcript_conditional_UT_without_best_hit[order(Transcript_conditional_UT_without_best_hit$FDR, decreasing = F),]
Transcript_conditional_UT_without_best_hit = Transcript_conditional_UT_without_best_hit[!duplicated(Transcript_conditional_UT_without_best_hit$phe_id),]
Transcript_conditional_UT_without_best_hit$gQTLs = paste(Transcript_conditional_UT_without_best_hit$phe_id, Transcript_conditional_UT_without_best_hit$var_id, sep = '_')

Transcript_conditional_IFN_without_best_hit = Transcript_conditional_IFN[-which(Transcript_conditional_IFN$phe_id %in% Transcript_conditional_IFN_best_hit$phe_id),]
Transcript_conditional_IFN_without_best_hit = Transcript_conditional_IFN_without_best_hit[order(Transcript_conditional_IFN_without_best_hit$FDR, decreasing = F),]
Transcript_conditional_IFN_without_best_hit = Transcript_conditional_IFN_without_best_hit[!duplicated(Transcript_conditional_IFN_without_best_hit$phe_id),]
Transcript_conditional_IFN_without_best_hit$gQTLs = paste(Transcript_conditional_IFN_without_best_hit$phe_id, Transcript_conditional_IFN_without_best_hit$var_id, sep = '_')

length(unique(Transcript_conditional_IFN_without_best_hit$phe_id))
length(unique(Transcript_conditional_UT_without_best_hit$phe_id))
length(unique(Transcript_conditional_LPS24_without_best_hit$phe_id))

#--- CS gQTLs
gQTLs <- list(A = Transcript_conditional_IFN$phe_id,
              B = Transcript_conditional_UT$phe_id,
              C = Transcript_conditional_LPS24$phe_id)
setdiffgQTLs = lapply(1:length(gQTLs), function(n) setdiff(gQTLs[[n]], unlist(gQTLs[-n])))
names(setdiffgQTLs) = c('IFN', 'UT', 'LPS24')

# moloc_gQTL = fread('moloc_all/moloc_nominal_complete_output.txt', stringsAsFactors = F, header = T)
# moloc_gQTL = as.data.frame(moloc_gQTL)
# moloc_gQTL_IFN = moloc_gQTL[which(moloc_gQTL$IFN.unique > 0.8),]
# moloc_gQTL_UT = moloc_gQTL[which(moloc_gQTL$UT.unique > 0.8),]
# moloc_gQTL_LPS24 = moloc_gQTL[which(moloc_gQTL$LPS24.unique > 0.8),]

#---  peak CS gQTLs
# LPS24 = rbind(Transcript_conditional_LPS24_best_hit, Transcript_conditional_LPS24_without_best_hit) 
# LPS24 = LPS24[which(LPS24$phe_id %in% c(setdiffgQTLs[['LPS24']], moloc_gQTL_LPS24$Gene )),]
# 
# IFN = rbind(Transcript_conditional_IFN_best_hit, Transcript_conditional_IFN_without_best_hit) 
# IFN = IFN[which(IFN$phe_id %in% c(setdiffgQTLs[['IFN']], moloc_gQTL_IFN$Gene )),]
# 
# UT = rbind(Transcript_conditional_UT_best_hit, Transcript_conditional_UT_without_best_hit) 
# UT = UT[which(UT$phe_id %in% c(setdiffgQTLs[['UT']], moloc_gQTL_UT$Gene )),]

LPS24 = rbind(Transcript_conditional_LPS24_best_hit, Transcript_conditional_LPS24_without_best_hit) 
LPS24 = LPS24[which(LPS24$phe_id %in% c(setdiffgQTLs[['LPS24']] )),]

IFN = rbind(Transcript_conditional_IFN_best_hit, Transcript_conditional_IFN_without_best_hit) 
IFN = IFN[which(IFN$phe_id %in% c(setdiffgQTLs[['IFN']] )),]

UT = rbind(Transcript_conditional_UT_best_hit, Transcript_conditional_UT_without_best_hit) 
UT = UT[which(UT$phe_id %in% c(setdiffgQTLs[['UT']] )),]

#------ per state

Query = get(treatment)
Query = as.data.frame(Query)

library(dplyr)
Query = merge(Query, SNPs, by = 'var_id')
Query = Query[order(Query$FDR, decreasing = F),]
Query = Query[!duplicated(Query$var_id),]

#-- If you donot use the LD peak SNPs
#-- the old version of the clumping function was removing some chromosomes

# library(ieugwasr) #R/3.5.0-newgcc
# clump_input = cbind.data.frame(rsid = Query$SNP_ID, chr_name = Query$SNP_CHROM, chrom_start = Query$SNP_POS, pval = Query$FDR)
# i <- sapply(clump_input, is.factor)
# clump_input[i] <- lapply(clump_input[i], as.character)
# str(clump_input)
# clumped_input <- ld_clump(clump_input)

# Query = Query[which(Query$SNP_ID %in% clumped_input$rsid),]

# INPUT = QUERY
INPUT = cbind.data.frame(SNP = Query$var_id,
                         beta.exposure = Query$bwd_slope,
                         se.exposure = Query$bwd_slope_se,
                         effect_allele.exposure = Query$REF,
                         other_allele.exposure = Query$ALT,
                         eaf.exposure = Query$maf,
                         pval.exposure = Query$FDR,
                         gene.exposure = Query$phe_id)

INPUT = cbind.data.frame(INPUT, rep('gQTL', dim(INPUT)[1]), rep('gQTL', dim(INPUT)[1]))
colnames(INPUT) = c('SNP','beta.exposure','se.exposure','effect_allele.exposure','other_allele.exposure','eaf.exposure','pval.exposure','gene.exposure','exposure','id.exposure')

i <- sapply(INPUT, is.factor)
INPUT[i] <- lapply(INPUT[i], as.character)

str(INPUT)
dim(INPUT)

dir.create(paste0('ciseQTL_trait_enrichment/CSgQTL_trait_colocalization_',treatment,'/mr'), recursive = T)
setwd(paste0('ciseQTL_trait_enrichment/CSgQTL_trait_colocalization_',treatment,'/'))

getwd()
ok=ok2=0
RESULTS = data.frame()

library(biomaRt)
snpdetail=useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp",host="grch37.ensembl.org", path="/biomart/martservice",archive=FALSE)
i=1

for(i in 1:length(allDBs$id))
{
  input_enrichment = chd_out_dat = NULL
  
  tryCatch({
    
    chd_out_dat <- extract_outcome_data(
      snps = INPUT$SNP,
      outcomes = allDBs$id[i]
    )
    
  }, error=function(e){})
  
  if(!is.null(chd_out_dat))
  {
    input_enrichment <- harmonise_data(
      exposure_dat = INPUT,
      outcome_dat = chd_out_dat
    )
    # colnames(INPUT)
    # chd_out_dat[,grep('proxy', colnames(chd_out_dat))]#[c(1:3,7:10,12,16)]
    
    input_enrichment = input_enrichment[which(input_enrichment$pval.outcome<1e-5),] # - report 1e-6 because 1e-5 < 7.03493e-06 

    #-------------- coloc
    if(dim(input_enrichment)[1]>0)
    {
      
      proxies = input_enrichment#[which(input_enrichment$proxy.outcome),]
      proxies$SNP = as.character(proxies$SNP)
      
      proxies$target_snp.outcome[is.na(proxies$target_snp.outcome)] = proxies$SNP[is.na(proxies$target_snp.outcome)]
      proxies$proxy_snp.outcome[is.na(proxies$proxy_snp.outcome)] = proxies$SNP[is.na(proxies$proxy_snp.outcome)]
      
      #----- coloc per SNP
      
      library(coloc)
      j=1
      for(j in 1:(dim(proxies)[1]))
      {
        
        if(is.na(proxies$eaf.outcome[j]))
        {
          SNPs_rsID_MAF = getBM(attributes=c('chr_name', 'chrom_start',"minor_allele_freq","refsnp_id"),filters="snp_filter",value=as.character(proxies$SNP[j]), snpdetail)
          if(dim(SNPs_rsID_MAF)[1]>0){proxies$eaf.outcome[j] = SNPs_rsID_MAF$minor_allele_freq}
        }
        
        if(is.na(proxies$eaf.exposure[j]))
        {
          SNPs_rsID_MAF = getBM(attributes=c('chr_name', 'chrom_start',"minor_allele_freq","refsnp_id"),filters="snp_filter",value=as.character(proxies$SNP[j]), snpdetail)
          if(dim(SNPs_rsID_MAF)[1]>0){proxies$eaf.outcome[j] = SNPs_rsID_MAF$minor_allele_freq}
        }
        
        tryCatch({
          
          COLOC_Per_SNP=NULL
          COLOC_Per_SNP <- coloc.abf(list(snp = proxies$SNP[j], pvalues=as.numeric(proxies$pval.exposure[j]), N=query_nCases, MAF=as.numeric(proxies$eaf.exposure[j]), type="quant"),
                                     list(snp = proxies$SNP[j], pvalues=as.numeric(proxies$pval.outcome[j]), N=proxies$samplesize.outcome[j], MAF=as.numeric(proxies$eaf.outcome[j]), type="quant"))
          
        }, error=function(e){})
        
        if(!is.null(COLOC_Per_SNP))
        {
          
          if(COLOC_Per_SNP$summary['PP.H4.abf']<0.5){print('1k'); input_enrichment = input_enrichment[-which(input_enrichment$proxy_snp.outcome == proxies$proxy_snp.outcome[j]),]}else{print('0k')}
          
          res = data.frame(proxies[j,], PP.H4.abf = as.numeric(COLOC_Per_SNP$summary['PP.H4.abf']))
          res = res[,which(colnames(res)%in%c('SNP','beta.exposure','beta.outcome','eaf.exposure','eaf.outcome','id.outcome', 'chr','pos','se.outcome','pval.outcome','outcome','originalname.outcome','target_snp.outcome','proxy_snp.outcome','pval.exposure', 'PP.H4.abf'))]
          
          print(paste0('dim(res): ', dim(res)[2] ))
          
          if(dim(res)[2]==16)
          {
            if(i==1 & ok == 0){RESULTS = res; ok=1}
            if(i!=1 & ok == 0){RESULTS = res; ok=1}
            if(i!=1 & ok == 1){RESULTS=rbind(res,RESULTS)}
          }
          
          print(paste0('dim(RESULTS): ', dim(RESULTS) ))
          
        }
        #print(RESULTS)
      }
    }

  #-------------- Mendelian randomization test
  if(dim(input_enrichment)[1]>0)
  {
    # remove duplicate summaries by selecting most informative one
    # input_enrichment<-power.prune(input_enrichment, method.size=F)
    
    # enrichment
    res2 <- mr(input_enrichment)
    
    if(i==1 & ok2 == 0){RESULTS2 = res2; ok2=1}
    if(i!=1 & ok2 == 0){RESULTS2 = res2; ok2=1}
    if(i!=1 & ok2 == 1){RESULTS2=rbind(res2, RESULTS2)}
   
    write.table(res2, paste0('mr/', allDBs$id[i], '_', treatment, '_Enrichment_traits_coloc_MRT.txt'), quote = F, row.names = F, sep = '\t')
  }
 }
}

fwrite(RESULTS2, paste0('CSgQTL_',treatment,'_Enrichment_traits_ALL.txt'), quote = F, row.names = F, sep = '\t')
fwrite(RESULTS, paste0('CSgQTL_',treatment,'_traits_colocalization.txt'), quote = F, row.names = F, sep = '\t')

length(unique(RESULTS[which(RESULTS$PP.H4.abf > 0.8), 'SNP']))
length(unique(RESULTS[which(RESULTS$PP.H4.abf > 0.8 & RESULTS$pval.outcome < 1e-8), 'SNP']))
RESULTS$originalname.outcome[grep('COVID', RESULTS$originalname.outcome)]



