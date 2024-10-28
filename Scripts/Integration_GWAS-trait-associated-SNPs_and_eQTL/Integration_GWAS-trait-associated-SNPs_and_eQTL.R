

############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
library(coloc)
library(biomaRt)
library(dplyr)
library(data.table)
library(coloc) 
library(qvalue)
library(TwoSampleMR)
library(MRInstruments)
library(haploR)

#----------------- Description -----------------
# The 380 GWAS summary statistics from UK-Biobank112 and MR Base GWAS databases8 (European-ancestry individuals) containing significant associations were prepared as instruments of GWAS trait enrichment analysis. 
# A causal LD proxy SNP was defined as an SNP in LD with the query SNP and causal in both GWAS and g/tQTL/mQTL studies. 
# We used co-localisation tests of two genetic traits (coloc) with default parameter for prior probabilities to estimate the Bayes factor as the posterior probability that an eSNP is causal in both GWAS and g/tQTL studies.
# The Mendelian randomisation (MR) test was applied to enrich the harmonized list of query SNPs and GWAS summary statistics. 

#----------------- Output -----------------
# Results of integration of GWAS trait-associated SNP and g/tQTL as text files.

#----------------- Examples -----------------
treatment = 'LPS24'
if(treatment == 'IFN'){query_nCases = 139} else if (treatment == 'LPS24') {query_nCases = 176 } else {query_nCases = 176 }

# Execute the rest of the code line by line.

##########################################################################################################################################################

# The list of available GWAS summary statistics
ao <- available_outcomes()
ao = ao[!is.na(ao$ncase),]

# List of available databases providing GWAS summary statistics.
DBs = as.data.frame(table(gsub('-.*','',ao$id)))

UKBIO = ao[grep('ukb', ao$id),]
UKBIOtrait = as.data.frame(table(UKBIO$trait))
UKBIO=UKBIO[grep('covid|COVID|Cancer|cancer|immun|Immun|inflammat|Inflammat', UKBIO$trait),]
UKBIO=UKBIO[-grep('non|Non|Illnesses of', UKBIO$trait),]
UKBIO <- UKBIO[grepl("European", UKBIO$population), ]

BBJ = ao[grep('bbj', ao$id),]
BBJtrait = as.data.frame(table(BBJ$trait))
BBJ <- BBJ[grepl("European", BBJ$population), ]   #mainly East Asian

EBI = ao[grep('ebi', ao$id),]
EBItrait = as.data.frame(table(EBI$trait))
EBI <- EBI[grepl("European", EBI$population), ]

ieu = ao[grep('ieu', ao$id),]
ieutrait = as.data.frame(table(ieu$trait))
ieu <- ieu[grepl("European", ieu$population), ]

table(ao$subcategory)  # Includes Psychiatric / neurological
toMatch <- c("Autoimmune / inflammatory",
             'Cofactors and vitamins',
             "Haemotological", "Cancer",
             "Immune system", "Immune cell subset frequency", "Metabolites ratio", "Nucleotide","Paediatric disease")
subcategories <- ao[grepl(paste(toMatch,collapse="|"), ao$subcategory),]
subcategories <- subcategories[grepl("European", subcategories$population), ]

allDBs = do.call("rbind", list(subcategories, UKBIO, BBJ, EBI, ieu))
allDBs = allDBs[!duplicated(allDBs$id),]

allDBs = allDBs[-grep('Illnesses of|self-reported|Self-reported|Non-cancer|non-cancer|noninflammatory|Noninflammatory', UKBIO$trait),]

#--------------------------- query
# treatment = 'LPS24'
# if(treatment == 'IFN'){query_nCases = 139} else if (treatment == 'LPS24') {query_nCases = 176 } else {query_nCases = 176 }
#---------------------------


# Obtain information from gQTL profiles
SNPs = fread('~/gQTL/QTLtools_inputs/MAF_imputed_Allsamples_TYPED.txt', stringsAsFactors = F, header = T, fill=TRUE)
SNPs = as.data.frame(SNPs)
SNPs = SNPs[,c('ID', 'REF', 'ALT', '0', '1', '2', 'mac', 'maf', 'TYPED')]
colnames(SNPs)[which(colnames(SNPs) %in% 'ID')] = 'var_id'
SNPs = SNPs[order(SNPs$maf, decreasing = T),]
SNPs = SNPs[!duplicated(SNPs$var_id),]

header = fread('~/gQTL/gQTL_nominal_Header.txt', stringsAsFactors = F, header = F)
headerC = fread('~/gQTL/gQTL_conditional_Header.txt', stringsAsFactors = F, header = F)
headerP = fread('~/gQTL/gQTL_permutation_Header.txt', stringsAsFactors = F, header = F)

Transcript_conditional_LPS24 = fread('QTLtools_outputs/LPS24/conditional_pass/gQTL_conditional_pass_1.txt', stringsAsFactors = F)
colnames(Transcript_conditional_LPS24) = as.character(headerC)
Transcript_conditional_LPS24 = Transcript_conditional_LPS24[-which(Transcript_conditional_LPS24$var_id == '.'),]
Transcript_conditional_LPS24$FDR <- qvalue(p = Transcript_conditional_LPS24$bwd_pval)$qvalues
Transcript_conditional_LPS24_sub = Transcript_conditional_LPS24[Transcript_conditional_LPS24$FDR < 0.01,]

Transcript_conditional_UT = fread('QTLtools_outputs/UT/conditional_pass/gQTL_conditional_pass_1.txt', stringsAsFactors = F)
colnames(Transcript_conditional_UT) = as.character(headerC)
Transcript_conditional_UT = Transcript_conditional_UT[-which(Transcript_conditional_UT$var_id == '.'),]
Transcript_conditional_UT$FDR <- qvalue(p = Transcript_conditional_UT$bwd_pval)$qvalues
Transcript_conditional_UT_sub = Transcript_conditional_UT[Transcript_conditional_UT$FDR < 0.01,]

Transcript_conditional_IFN = fread('QTLtools_outputs/IFN/conditional_pass/gQTL_conditional_pass_1.txt', stringsAsFactors = F)
colnames(Transcript_conditional_IFN) = as.character(headerC)
Transcript_conditional_IFN = Transcript_conditional_IFN[-which(Transcript_conditional_IFN$var_id == '.'),]
Transcript_conditional_IFN$FDR <- qvalue(p = Transcript_conditional_IFN$bwd_pval)$qvalues
Transcript_conditional_IFN_sub = Transcript_conditional_IFN[Transcript_conditional_IFN$FDR < 0.01,]

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

#--- Locate gQTLs that are specific to each context across states
gQTLs <- list(A = Transcript_conditional_IFN$phe_id,
              B = Transcript_conditional_UT$phe_id,
              C = Transcript_conditional_LPS24$phe_id)
setdiffgQTLs = lapply(1:length(gQTLs), function(n) setdiff(gQTLs[[n]], unlist(gQTLs[-n])))
names(setdiffgQTLs) = c('IFN', 'UT', 'LPS24')

LPS24 = rbind(Transcript_conditional_LPS24_best_hit, Transcript_conditional_LPS24_without_best_hit) 
LPS24 = LPS24[which(LPS24$phe_id %in% c(setdiffgQTLs[['LPS24']] )),]

IFN = rbind(Transcript_conditional_IFN_best_hit, Transcript_conditional_IFN_without_best_hit) 
IFN = IFN[which(IFN$phe_id %in% c(setdiffgQTLs[['IFN']] )),]

UT = rbind(Transcript_conditional_UT_best_hit, Transcript_conditional_UT_without_best_hit) 
UT = UT[which(UT$phe_id %in% c(setdiffgQTLs[['UT']] )),]

#------ Prepare query
Query = get(treatment)
Query = as.data.frame(Query)

Query = merge(Query, SNPs, by = 'var_id')
Query = Query[order(Query$FDR, decreasing = F),]
Query = Query[!duplicated(Query$var_id),]
#------ 

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

dir.create(paste0('ciseQTL_trait_enrichment/CSgQTL_trait_colocalization_',treatment,'/mr'), recursive = T)
setwd(paste0('ciseQTL_trait_enrichment/CSgQTL_trait_colocalization_',treatment,'/'))

ok=ok2=0
RESULTS = data.frame()

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
 
    input_enrichment = input_enrichment[which(input_enrichment$pval.outcome<1e-5),] # - report 1e-6 because 1e-5 < 7.03493e-06 

    #-------------- coloc analysis
    if(dim(input_enrichment)[1]>0)
    {
      
      proxies = input_enrichment#[which(input_enrichment$proxy.outcome),]
      proxies$SNP = as.character(proxies$SNP)
      
      proxies$target_snp.outcome[is.na(proxies$target_snp.outcome)] = proxies$SNP[is.na(proxies$target_snp.outcome)]
      proxies$proxy_snp.outcome[is.na(proxies$proxy_snp.outcome)] = proxies$SNP[is.na(proxies$proxy_snp.outcome)]
      
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
