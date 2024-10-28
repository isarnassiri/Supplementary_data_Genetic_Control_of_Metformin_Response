
############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(network))
suppressPackageStartupMessages(library(sna))
suppressPackageStartupMessages(library(ggnetwork))

#----------------- Description -----------------
# A coExQTL's functional impact would be vary depending on the specific context and the biological significance of the gene within the pathway. 
# We used a pathway-based differential co-expression analysis to determine the permissible range of expression variation within a pathway that would classify a coExQTL as functionally disruptive. 
# Pathway-based differential co-expression analysis was performed on a subset of genes that have allele-specific co-expression relationships and enriched with curated gene sets from online pathway databases (FDR < 0.01). 
# The analysis uncovers the average shift in correlation between gene expression among two genotype classes, as well as its statistical significance.

#----------------- Output -----------------
# Results of pathway enrichment analysis and Pathway-based differential co-expression analysis as text files.

#----------------- Examples -----------------
treatment='LPS24'

# It is necessary to call the functions in lines 39-172. To do that, select the lines and hit the Run button/Enter.
# Execute the rest of the code line by line.

##########################################################################################################################################################

## -------- FUNCTIONS (start) --------

get_p_from_r2n <- function(r2, n)
{
  fval <- r2 * (n-2) / (1 - r2)
  pval <- pf(fval, 1, n-1, low=FALSE)
  return(pval)
}

get_r_from_pn <- function(p, n)
{
  optim.get_p_from_rn <- function(x, sample_size, pvalue)
  {
    abs(-log10(get_p_from_r2n(x, sample_size)) - -log10(pvalue))
  }
  
  if(length(p) > 1 & length(n) == 1)
  {
    message("Assuming n the same for all p values")
    n <- rep(n, length(p))
  }
  
  Fval <- qf(p, 1, n-1, low=FALSE)
  R2 <- Fval / (n - 2 + Fval)
  index <- !is.finite(Fval)
  if(any(index))
  {
    index <- which(index)
    for(i in 1:length(index))
    {
      R2[index[i]] <- optim(0.001, optim.get_p_from_rn, sample_size=n[index[i]], pvalue=p[index[i]])$par
    }
  }
  return(sqrt(R2))
}


steiger_sensitivity <- function(rgx_o, rgy_o, ...)
{
  if(rgy_o > rgx_o)
  {
    a <- rgy_o
    b <- rgx_o
  } else {
    a <- rgx_o
    b <- rgy_o
  }
  
  mycolors.trans = rgb(c(0,0), c(0,0), 
                       c(0,255),alpha = c(70,255), maxColorValue = 255) 
  
  d <- expand.grid(rxx_o=seq(rgx_o,1,length.out=50), ryy_o=seq(rgy_o,1,length.out=50), type=c("A","B"))
  d$rgy <- rgy_o / d$ryy_o
  d$rgx <- rgx_o / d$rxx_o
  d$z <- d$rgy - d$rgx
  d$z[d$type=="A"] <- 0
  temp <- wireframe(
    z ~ rxx_o * ryy_o, 
    groups=type, 
    data=d, 
    scales=list(arrows=FALSE), 
    col.groups = mycolors.trans, 
    drape=FALSE, 
    ylab=expression(rho[xx[o]]), 
    xlab=expression(rho[yy[o]]),
    zlab=expression(rho[gy]-rho[gx]),
    par.settings = list(axis.line=list(col="transparent")),
    ...
  )
  
  vz <- a * log(a) - b * log(b) + a*b*(log(b)-log(a))
  vz0 <- -2*b - b * log(a) - a*b*log(a) + 2*a*b
  
  vz1 <- abs(vz - vz0)
  
  sensitivity <- vz0 / (2 * vz0 + abs(vz))
  sensitivity_ratio <- vz1 / vz0
  
  return(list(
    vz = vz,
    vz0 = vz0,
    vz1 = vz1,
    sensitivity = sensitivity,
    sensitivity_ratio = sensitivity_ratio,
    pl = temp
  ))
}

mr_steiger <- function(p_exp, p_out, n_exp, n_out, r_xxo = 1, r_yyo=1, ...) 
{
  require(psych)
  index <- any(is.na(p_exp)) | any(is.na(p_out)) | any(is.na(n_exp)) | any(is.na(n_out))
  p_exp <- p_exp[!index]
  p_out <- p_out[!index]
  n_exp <- n_exp[!index]
  n_out <- n_out[!index]
  r_exp <- get_r_from_pn(p_exp, n_exp)
  r_out <- get_r_from_pn(p_out, n_out)
  
  r_exp_adj <- sqrt(r_exp^2 / r_xxo^2)
  r_out_adj <- sqrt(r_out^2 / r_yyo^2)
  
  sensitivity <- steiger_sensitivity(r_exp, r_out, ...)
  
  rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
  rtest_adj <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp_adj, r34 = r_out_adj)
  l <- list(
    r2_exp = r_exp^2, 
    r2_out = r_out^2, 
    r2_exp_adj = r_exp_adj^2, 
    r2_out_adj = r_out_adj^2, 
    correct_causal_direction = r_exp > r_out, 
    steiger_test = rtest$p,
    correct_causal_direction_adj = r_exp_adj > r_out_adj, 
    steiger_test_adj = rtest_adj$p,
    vz = sensitivity$vz,
    vz0 = sensitivity$vz0,
    vz1 = sensitivity$vz1,
    sensitivity = sensitivity$sensitivity,
    sensitivity_ratio = sensitivity$sensitivity_ratio,
    sensitivity_plot = sensitivity$pl
  )
  return(l)
}

mr_wald_ratio <- function(b_exp, b_out, se_exp, se_out, parameters) 
{
  if (length(b_exp) > 1) {
    return(list(b = NA, se = NA, pval = NA, nsnp = NA))
  }
  b <- b_out/b_exp
  se <- se_out/abs(b_exp)
  pval <- pnorm(abs(b)/se, lower.tail = F) * 2
  return(list(b = b, se = se, pval = pval, nsnp = 1))
}

## -------- FUNCTIONS (end) --------

# read in met-exp correlation
library(data.table)
methExprs = fread(paste0(treatment,'correlationMethExprs_results.txt'), stringsAsFactors = F, header = T)
methExprs_sub = methExprs[which(methExprs$adj.P.Val<0.001),]
#  beta: coefficient of the methylation change

# read in metQTL profile
metQTL_LPS24 = fread('~/nominal_imputedGenotype_all_LPS_annotated.txt', stringsAsFactors = F)
dim(metQTL_LPS24)
colnames(metQTL_LPS24)[1] = c('cpg_ID')

metQTL_LPS24_sub = metQTL_LPS24[metQTL_LPS24$FDR < 0.001,]
dim(metQTL_LPS24_sub)
head(metQTL_LPS24_sub)
metQTL_LPS24_sub$SNP_ID = gsub('.*;', '', metQTL_LPS24_sub$SNP_ID)
table(metQTL_LPS24_sub$chrID)

# read in gQTLs (LPS)
eQTL_LPS24 = fread('/Users/isarnassiri/Documents/RESULTS_USED/eQTL/eQTL_LPS24_annotated_clumped_leadSNPs.txt', stringsAsFactors = F)
eQTL_LPS24_sub = eQTL_LPS24[eQTL_LPS24$FDR < 0.001,]
dim(eQTL_LPS24_sub)

colnames(eQTL_LPS24_sub)[match(c("FDR", "slop"), colnames(eQTL_LPS24_sub))] = c("FDR_eQTL", "slope_eQTL")
colnames(metQTL_LPS24_sub)[match(c("FDR", "slope"), colnames(metQTL_LPS24_sub))] = c("FDR_mQTL", "slope_mQTL")

eQTL_mQTL_all = merge(eQTL_LPS24_sub, metQTL_LPS24_sub, by = 'SNP_ID')
eQTL_mQTL_all = as.data.frame(eQTL_mQTL_all)
eQTL_mQTL_all = eQTL_mQTL_all[,-which(colnames(eQTL_mQTL_all) == 'index')]

eQTL_mQTL_all$index = paste(eQTL_mQTL_all$cpg_ID, eQTL_mQTL_all$gene_name, sep = '_')
methExprs_sub$index = paste(methExprs_sub$cpg, methExprs_sub$exprs, sep = '_')

eQTL_mQTL_corr = merge(eQTL_mQTL_all, methExprs_sub, by = 'index')
eQTL_mQTL_corr$index2 = paste(eQTL_mQTL_corr$SNP_ID, eQTL_mQTL_corr$gene_id, sep = '_')

#== read in gQTL profile (Naive)
library(data.table)
eQTL_Naive = fread('/Users/isarnassiri/Documents/RESULTS_USED/eQTL/eQTL_UT_annotated_clumped_leadSNPs.txt', stringsAsFactors = F)
eQTL_Naive_sub = eQTL_Naive[eQTL_Naive$FDR < 0.001,]
dim(eQTL_Naive_sub)
eQTL_Naive_sub$index = paste(eQTL_Naive_sub$SNP_ID, eQTL_Naive_sub$gene_id, sep = '_')
dim(eQTL_mQTL_corr) 

# Find eQTL_mQTL_pairs with no correlation in Naive

met_exp_UT = fread(paste0('UTcorrelationMethExprs_results.txt'), stringsAsFactors = F ,header = T)
range(met_exp_UT$adj.P.Val)
met_exp_UT = met_exp_UT[met_exp_UT$adj.P.Val < 0.001,]

met_exp_UT$index = paste(met_exp_UT$cpg, met_exp_UT$exprs, sep = '_')

eQTL_mQTL_corr_CS = eQTL_mQTL_corr[-which(eQTL_mQTL_corr$index %in% met_exp_UT$index),]
dim(eQTL_mQTL_corr_CS)

library(data.table)
getwd()
fwrite(eQTL_mQTL_corr_CS, 'eQTL_mQTL_corr_CS_coloc_LPS_Supp5.txt', quote = F, row.names = F, sep = '\t')

eQTL_mQTL_corr_CS[(eQTL_mQTL_corr_CS$exprs == 'CD55'),]
eQTL_mQTL_corr_CS[grep('CD55', eQTL_mQTL_corr_CS$gene_name),]

# colocalisation analyses per SNP
library(coloc)

ok=0
j=1
for(j in 1:(dim(eQTL_mQTL_corr_CS)[1]))
{
COLOC_Per_SNP=NULL
COLOC_Per_SNP <- coloc.abf(list(snp = eQTL_mQTL_corr_CS$SNP_ID[j], pvalues=as.numeric(eQTL_mQTL_corr_CS$FDR_eQTL[j]), N=176, MAF=as.numeric(eQTL_mQTL_corr_CS$maf[j]), type="quant"), 
                           list(snp = eQTL_mQTL_corr_CS$SNP_ID[j], pvalues=as.numeric(eQTL_mQTL_corr_CS$FDR_mQTL[j]), N=176, MAF=as.numeric(eQTL_mQTL_corr_CS$maf[j]), type="quant"))

temp = data.frame(cbind(eQTL_mQTL_corr_CS[j,], t(as.data.frame(COLOC_Per_SNP$summary[-1])) ), SNP.PP.H4 = COLOC_Per_SNP$results$SNP.PP.H4)

  if(j==1 & ok == 0){RESULTS = temp; ok=1}
  if(j!=1 & ok == 0){RESULTS = temp; ok=1}
  if(j!=1 & ok == 1){RESULTS=rbind(temp,RESULTS)} 

print(j)
}

colnames(RESULTS) = gsub('\\.y', '_mQTL', colnames(RESULTS))
colnames(RESULTS) = gsub('\\.x', '_gQTL', colnames(RESULTS))
colnames(RESULTS)

fwrite(RESULTS, paste0('eQTL_mQTL_corr_CS_coloc_',treatment,'.txt'), quote = F, row.names = F, sep = '\t')

eQTL_mQTL_corr_CS_coloc = RESULTS
eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8),]

# Read the gQTL profile (nominal pass) to extract standard errors for gQTLs.
header = fread('~/gQTL/gQTL_nominal_Header.txt', stringsAsFactors = F, header = F)
Gene_nominal_LPS24 = fread('~/gQTL/QTLtools_outputs/LPS24/nominal_pass/gQTL_nominal_pass_1.txt', stringsAsFactors = F)
Gene_nominal_LPS24 = Gene_nominal_LPS24[which(Gene_nominal_LPS24$var_id %in% eQTL_mQTL_corr_CS_coloc$SNP_ID),]
Gene_nominal_LPS24$index = paste(Gene_nominal_LPS24$var_id, Gene_nominal_LPS24$phe_id, sep = '_')

eQTL_mQTL_corr_CS_coloc$index = paste(eQTL_mQTL_corr_CS_coloc$SNP_ID, eQTL_mQTL_corr_CS_coloc$gene_name, sep = '_')
eQTL_mQTL_corr_CS_coloc_interset = merge(eQTL_mQTL_corr_CS_coloc, Gene_nominal_LPS24[,c('slope_se', 'index')], by = 'index')

eQTL_mQTL_corr_CS_coloc$slope_se = summary(eQTL_mQTL_corr_CS_coloc_interset$slope_se)[4]

eQTL_mQTL_corr_CS_coloc_interset = as.data.frame(eQTL_mQTL_corr_CS_coloc_interset)
eQTL_mQTL_corr_CS_coloc = as.data.frame(eQTL_mQTL_corr_CS_coloc)

eQTL_mQTL_corr_CS_coloc = rbind(eQTL_mQTL_corr_CS_coloc_interset, eQTL_mQTL_corr_CS_coloc[-which(eQTL_mQTL_corr_CS_coloc$index %in% eQTL_mQTL_corr_CS_coloc_interset$index),])

#------------
qtl_orig = eQTL_mQTL_corr_CS_coloc[,c("SNP_ID", "gene_name", "cpg_ID", "FDR_eQTL", "FDR_mQTL", "REF", "ALT", "slope_eQTL", "slope_mQTL", 'slope_se', 'se')]
colnames(qtl_orig) = c("SNP_ID", "exp_ind" , "meth_ind", "PVALUE_expr", "PVALUE_meth", "AL1_expr", "AL2_expr", 'EFFECT_expr', 'EFFECT_meth', 'SE_expr', 'SE_meth')

qtl_orig = as.data.frame(qtl_orig)
me <- expand.grid(rxx_o=c(1), ryy_o=c(1))
qtl <- group_by(me, rxx_o, ryy_o) %>%
  do({
    x <- .
    y <- qtl_orig
    y$rxx_o <- x$rxx_o[1]
    y$ryy_o <- x$ryy_o[1]
    y
  })

qtl$r_exp <- NA
qtl$r_meth <- NA
qtl$dir <- NA
qtl$dir_p <- NA
qtl$mr_eff <- NA
qtl$mr_se <- NA
qtl$mr_p <- NA

for(i in 1:nrow(qtl))
{
  l <- mr_steiger(
    qtl$PVALUE_meth[i], 
    qtl$PVALUE_expr[i], 
    176, 176,
    qtl$rxx_o[i],
    qtl$ryy_o[i],
    screen = list(z = 70, x = -60, y = 3)
  )
  
  l$sensitivity_plot
  qtl$r_exp[i] <- l$r2_out_adj
  qtl$r_meth[i] <- l$r2_exp_adj
  qtl$steiger[i] <- l$steiger_test
  qtl$dir[i] <- l$correct_causal_direction_adj
  qtl$dir_p[i] <- l$steiger_test_adj
  qtl$sensitivity[i] <- l$sensitivity
  qtl$reliability[i] <- l$sensitivity_ratio
  
  #The Wald test is a rough approximation of the Likelihood Ratio Test.
  if(qtl$dir[i])
  {
    a <- mr_wald_ratio(qtl$EFFECT_meth[i], qtl$EFFECT_expr[i], qtl$SE_meth[i], qtl$SE_expr[i])
    qtl$mr_eff[i] <- a$b
    qtl$mr_se[i] <- a$se
    qtl$mr_p[i] <- a$pval
  } else {
    a <- mr_wald_ratio(qtl$EFFECT_expr[i], qtl$EFFECT_meth[i], qtl$SE_expr[i], qtl$SE_meth[i])
    qtl$mr_eff[i] <- a$b
    qtl$mr_se[i] <- a$se
    qtl$mr_p[i] <- a$pval
  }
}

# Does expression or methylation most likely have the causal effect?
qtl$signSlopMeth = sign(qtl$EFFECT_meth)
qtl$signSlopExp = sign(qtl$EFFECT_expr)

fwrite(qtl, paste0('eQTL_mQTL_corr_CS_mr_steiger_',treatment,'.txt'), quote = F, row.names = F, sep = '\t')
