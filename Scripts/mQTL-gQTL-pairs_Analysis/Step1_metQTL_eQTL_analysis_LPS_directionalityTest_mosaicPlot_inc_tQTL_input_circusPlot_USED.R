
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/eQTL_mQTL_pairs_analysis/')

#==== read in
treatment='LPS24'

#=== met-exp correlation
library(data.table)
methExprs = fread(paste0(treatment,'correlationMethExprs_results.txt'), stringsAsFactors = F, header = T)
methExprs_sub = methExprs[which(methExprs$adj.P.Val<0.001),]
#  beta: coefficient of the methylation change

#== metQTL all
metQTL_LPS24 = fread('/Users/isarnassiri/Documents/RESULTS_USED/methylationQTL/QTLtools_PC_Selection/nominal_imputedGenotype_all_LPS_annotated.txt', stringsAsFactors = F)
dim(metQTL_LPS24)
colnames(metQTL_LPS24)[1] = c('cpg_ID')

metQTL_LPS24_sub = metQTL_LPS24[metQTL_LPS24$FDR < 0.001,]
dim(metQTL_LPS24_sub)
head(metQTL_LPS24_sub)
metQTL_LPS24_sub$SNP_ID = gsub('.*;', '', metQTL_LPS24_sub$SNP_ID)
table(metQTL_LPS24_sub$chrID)

#== eQTL peak LPS
eQTL_LPS24 = fread('/Users/isarnassiri/Documents/RESULTS_USED/eQTL/eQTL_LPS24_annotated_clumped_leadSNPs.txt', stringsAsFactors = F)
eQTL_LPS24_sub = eQTL_LPS24[eQTL_LPS24$FDR < 0.001,]
dim(eQTL_LPS24_sub)

#== eQTL_mQTL_pairs_all
colnames(eQTL_LPS24_sub)[match(c("FDR", "slop"), colnames(eQTL_LPS24_sub))] = c("FDR_eQTL", "slope_eQTL")
colnames(metQTL_LPS24_sub)[match(c("FDR", "slope"), colnames(metQTL_LPS24_sub))] = c("FDR_mQTL", "slope_mQTL")

eQTL_mQTL_all = merge(eQTL_LPS24_sub, metQTL_LPS24_sub, by = 'SNP_ID')
dim(eQTL_mQTL_all)
head(eQTL_mQTL_all)

eQTL_mQTL_all = as.data.frame(eQTL_mQTL_all)
eQTL_mQTL_all = eQTL_mQTL_all[,-which(colnames(eQTL_mQTL_all) == 'index')]

#== eQTL_mQTL_pairs_correlated all
eQTL_mQTL_all$index = paste(eQTL_mQTL_all$cpg_ID, eQTL_mQTL_all$gene_name, sep = '_')
methExprs_sub$index = paste(methExprs_sub$cpg, methExprs_sub$exprs, sep = '_')

eQTL_mQTL_corr = merge(eQTL_mQTL_all, methExprs_sub, by = 'index')
dim(eQTL_mQTL_corr)

eQTL_mQTL_corr$index2 = paste(eQTL_mQTL_corr$SNP_ID, eQTL_mQTL_corr$gene_id, sep = '_')

#== eQTL peak Naive
library(data.table)
eQTL_Naive = fread('/Users/isarnassiri/Documents/RESULTS_USED/eQTL/eQTL_UT_annotated_clumped_leadSNPs.txt', stringsAsFactors = F)
eQTL_Naive_sub = eQTL_Naive[eQTL_Naive$FDR < 0.001,]
dim(eQTL_Naive_sub)
eQTL_Naive_sub$index = paste(eQTL_Naive_sub$SNP_ID, eQTL_Naive_sub$gene_id, sep = '_')
dim(eQTL_mQTL_corr) 

#=============================== eQTL_mQTL_pairs_correlated + no corr in naive

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

#=============================== coloc per SNP
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

library(data.table)
fwrite(RESULTS, paste0('eQTL_mQTL_corr_CS_coloc_',treatment,'.txt'), quote = F, row.names = F, sep = '\t')

#== In total, we find 1416  eQTL/mQTL pairs 
dim(RESULTS)

#== of which 0.899% (1273/1416) are observed to share one causal variant affecting both expression and methylation (posterior probabilities of H4 > 0.75)
table(RESULTS$PP.H4.abf>0.8)


View(RESULTS[which(RESULTS$SNP_ID == 'rs2914937'),])

#----------------- do I have any trans-effect or no
# head(featureData)
# dim(featureData)
# featureData_sub2 = featureData
# colnames(featureData_sub2) = paste0('gene_',colnames(featureData_sub2))
# colnames(featureData_sub2)[1] = 'exprs'
# 
# methExprs_annotated = merge(methExprs_sub, featureData_sub2, by = 'exprs')
# 
# dim(cpg_featureData)
# cpg_featureData_sub2 = cpg_featureData
# colnames(cpg_featureData_sub2) = paste0('cpg_',colnames(cpg_featureData_sub2))
# colnames(cpg_featureData_sub2)[1] = 'cpg'
# 
# methExprs_annotated = merge(methExprs_annotated, cpg_featureData_sub2, by = 'cpg')
# 
# methExprs_annotated$cpg_chromosome = gsub('_.*','',methExprs_annotated$cpg_chromosome)
# table(methExprs_annotated$gene_chromosome == methExprs_annotated$cpg_chromosome)
# methExprs_annotated_distal = methExprs_annotated[which(methExprs_annotated$gene_chromosome != methExprs_annotated$cpg_chromosom),]
# 

#--- CHECK PVALUE FOR FIGURE
# eQTL_mQTL_corr_CS_coloc = fread('eQTL_mQTL_corr_CS_coloc.txt', stringsAsFactors = F ,header = T)
# eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$cpg_ID == 'cg22687766' &  eQTL_mQTL_corr_CS_coloc$SNP_ID == 'rs2914937'),]
# 
# 
# MethyQTL_UT = fread('nominal_imputedGenotype_all_UT_annotated.txt', stringsAsFactors = F, header = T, fill = T)
# MethyQTL_UT$SNP_ID = gsub('.*;', '', MethyQTL_UT$SNP_ID)
# MethyQTL_UT[which(MethyQTL_UT$CpG == 'cg22687766' &  MethyQTL_UT$SNP_ID == 'rs2914937'),]


#=============================== steiger test ===============================
## ---- FUNCTIONS ----
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

##========================================================= shakhbazov =========================================================
library(data.table)
eQTL_mQTL_corr_CS_coloc = fread(paste0('eQTL_mQTL_corr_CS_coloc_',treatment,'.txt'), stringsAsFactors = F ,header = T)
# eQTL_mQTL_corr_CS_coloc = fread(paste0('/Users/isar.nassiri/Desktop/Analysis_monocyte/KEY-FILES/Article_figures/FIGURES_TABLES-ciseQTLs/Exp_Met_analysis/eQTL_mQTL_PermutatedPeakSNPs_CS_coloc',treatment,'.txt'), stringsAsFactors = F ,header = T)
colnames(eQTL_mQTL_corr_CS_coloc)
dim(eQTL_mQTL_corr_CS_coloc)
eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8),]
table(eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8)
colnames(eQTL_mQTL_corr_CS_coloc)
dim(eQTL_mQTL_corr_CS_coloc)

table(eQTL_mQTL_corr_CS_coloc$FDR_eQTL < 1e-40)
eQTL_mQTL_corr_CS_coloc[grep('CD55', eQTL_mQTL_corr_CS_coloc$gene_name),]

dim(table(eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$FDR_eQTL < 1e-40),'chrID']))

#------------ add se ------------
#== metQTL nominal

library(data.table)
header = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/gQTL_nominal_Header.txt', stringsAsFactors = F, header = F)

#-- LPS nominal
Gene_nominal_LPS24 = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/QTLtools_outputs/LPS24/nominal_pass/gQTL_nominal_pass_1.txt', stringsAsFactors = F)
colnames(Gene_nominal_LPS24) = as.character(header)
Gene_nominal_LPS24 = Gene_nominal_LPS24[which(Gene_nominal_LPS24$var_id %in% eQTL_mQTL_corr_CS_coloc$SNP_ID),]
Gene_nominal_LPS24$index = paste(Gene_nominal_LPS24$var_id, Gene_nominal_LPS24$phe_id, sep = '_')

eQTL_mQTL_corr_CS_coloc$index = paste(eQTL_mQTL_corr_CS_coloc$SNP_ID, eQTL_mQTL_corr_CS_coloc$gene_name, sep = '_')
eQTL_mQTL_corr_CS_coloc_interset = merge(eQTL_mQTL_corr_CS_coloc, Gene_nominal_LPS24[,c('slope_se', 'index')], by = 'index')
dim(eQTL_mQTL_corr_CS_coloc_interset)

#-- I do not have se for all gQTLs, therefore I use mean
eQTL_mQTL_corr_CS_coloc$slope_se = summary(eQTL_mQTL_corr_CS_coloc_interset$slope_se)[4]

identical(colnames(eQTL_mQTL_corr_CS_coloc_interset), colnames(eQTL_mQTL_corr_CS_coloc))
class(eQTL_mQTL_corr_CS_coloc_interset)
class(eQTL_mQTL_corr_CS_coloc)

eQTL_mQTL_corr_CS_coloc_interset = as.data.frame(eQTL_mQTL_corr_CS_coloc_interset)
eQTL_mQTL_corr_CS_coloc = as.data.frame(eQTL_mQTL_corr_CS_coloc)

eQTL_mQTL_corr_CS_coloc = rbind(eQTL_mQTL_corr_CS_coloc_interset, eQTL_mQTL_corr_CS_coloc[-which(eQTL_mQTL_corr_CS_coloc$index %in% eQTL_mQTL_corr_CS_coloc_interset$index),])
colnames(eQTL_mQTL_corr_CS_coloc)

#------------
qtl_orig = eQTL_mQTL_corr_CS_coloc[,c("SNP_ID", "gene_name", "cpg_ID", "FDR_eQTL", "FDR_mQTL", "REF", "ALT", "slope_eQTL", "slope_mQTL", 'slope_se', 'se')]
colnames(qtl_orig) = c("SNP_ID", "exp_ind" , "meth_ind", "PVALUE_expr", "PVALUE_meth", "AL1_expr", "AL2_expr", 'EFFECT_expr', 'EFFECT_meth', 'SE_expr', 'SE_meth')

# qtl_orig$SE_meth = rep(0.01, dim(qtl_orig)[1])
# qtl_orig$SE_expr = rep(0.01, dim(qtl_orig)[1])

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
i=1
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
dim(qtl)

# Focusing upon the post-LPS state we find that 853/1273 (67%) of these loci demonstrate dependent causal effects on methylation linked to expression or vice versa (Figure 3b). 
table(qtl$dir_p<0.05) # steiger is equal to dir_p
# c('FALSE', 'TRUE'))] = c("Expression causes Methylation", "Methylation causes Expression")
table(qtl$mr_p<0.05) 
# qtl = qtl[which(qtl$mr_p<0.05),]
dim(qtl)

View(qtl[which(qtl$SNP_ID == 'rs2914937'),])

getwd()
library(data.table)
fwrite(qtl, paste0('eQTL_mQTL_corr_CS_mr_steiger_',treatment,'.txt'), quote = F, row.names = F, sep = '\t')

##=========================================================
INDEP = qtl[which(qtl$dir_p>0.05),]
table(qtl$dir)
length(which(qtl$dir_p>0.05)) # We find that commonly a SNP affects gene expression and DNA methylation dependently 
dim(qtl)
table((INDEP$mr_eff>=0))
table((INDEP$mr_eff>=0))
dim(INDEP)
table(INDEP$signSlopMeth == INDEP$signSlopExp)
INDEP = data.frame(table(INDEP$signSlopMeth == INDEP$signSlopExp))
INDEP$Var1 = as.character(INDEP$Var1)
INDEP$Var1[which(INDEP$Var1==F)] = 'Negative_Effect'
INDEP$Var1[which(INDEP$Var1==T)] = 'Possitive_Effect'
row.names = INDEP$Var1
INDEP = data.frame(INDEP = INDEP[,-1])
row.names(INDEP) = row.names

nonINDEP = subset(qtl, dir_p <= 0.05) %>%
  group_by(dir) %>% mutate(pos = ifelse(signSlopMeth == signSlopExp, 'Possitive_Effect', 'Negative_Effect'))
dim(nonINDEP)
nonINDEP_table <- table(nonINDEP$pos, nonINDEP$dir) 
colnames(nonINDEP_table)[which(colnames(nonINDEP_table) == c('FALSE', 'TRUE'))] = c("Expression causes Methylation", "Methylation causes Expression")
nonINDEP = data.frame(rbind(nonINDEP_table))
class(nonINDEP)

sum(nonINDEP)
sum(INDEP)

input2 = cbind(INDEP, nonINDEP)
input2

dev.off()
library("vcd")
mosaicplot(input2, shade=T, color=T, las=1, type="pearson", cex.axis=1, main="")

## Loading required package: grid
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/Figures_draft')
pdf(paste0('mosaicplot_revised',treatment,'.pdf'), width = 10, height = 8, useDingbats = F)
mosaicplot(input2, shade=F, color=F, las=1, type="pearson", cex.axis=1, main="")
dev.off()
getwd()
#

round(input2/sum(input2), digits = 3)*100


#=============================== circus plot visualization ===============================
library(data.table)
treatment='LPS24'
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/eQTL_mQTL_pairs_analysis/')

eQTL_mQTL_corr_CS_coloc = fread(paste0('eQTL_mQTL_corr_CS_coloc_',treatment,'.txt'), stringsAsFactors = F ,header = T)
length(unique(eQTL_mQTL_corr_CS_coloc$SNP_ID))
table(eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8)
dim(eQTL_mQTL_corr_CS_coloc)
table(paste0('chr', eQTL_mQTL_corr_CS_coloc$seqnames[order(eQTL_mQTL_corr_CS_coloc$seqnames)]))

eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8,]
eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[order(eQTL_mQTL_corr_CS_coloc$bpval, decreasing = F),]
eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[!duplicated(eQTL_mQTL_corr_CS_coloc$gene_id),]
dim(eQTL_mQTL_corr_CS_coloc) # I remove duplicated genes

eQTL_mQTL_corr_CS_coloc = as.data.frame(eQTL_mQTL_corr_CS_coloc)
table(eQTL_mQTL_corr_CS_coloc[,colnames(eQTL_mQTL_corr_CS_coloc) == 'seqnames'])
eQTL_mQTL_corr_CS_coloc[,colnames(eQTL_mQTL_corr_CS_coloc) == 'seqnames'] = factor(eQTL_mQTL_corr_CS_coloc[,colnames(eQTL_mQTL_corr_CS_coloc) == 'seqnames'], levels = c(1:22))
colnames(eQTL_mQTL_corr_CS_coloc)
dim(eQTL_mQTL_corr_CS_coloc)
table(eQTL_mQTL_corr_CS_coloc$seqnames)
eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[with(eQTL_mQTL_corr_CS_coloc, order(seqnames, start_gQTL)),]
colnames(eQTL_mQTL_corr_CS_coloc)

summary(eQTL_mQTL_corr_CS_coloc$FDR_eQTL)

dim(eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$FDR_eQTL < 1e-40), ])

#=============================== input files for visualization

#-------------------------------------- labels
colnames(eQTL_mQTL_corr_CS_coloc)
gene_label = eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$FDR_eQTL < 1e-40),c('seqnames','start_gQTL','end_gQTL',"gene_name")]
CPG_label = eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$FDR_eQTL < 1e-40),c('chrID','start_mQTL','end_mQTL',"cpg")]

colnames(gene_label) = c('chr',	'start',	'end', 'label')
colnames(CPG_label) = c('chr',	'start',	'end',	'label')

dim(gene_label)
dim(CPG_label)

head(gene_label)

# eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[!is.na(eQTL_mQTL_corr_CS_coloc$gene_id.y),]
# eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$FDR<1e-12),]
# transcript_tQTL = eQTL_mQTL_corr_CS_coloc[,c('seqnames.y','start','end', "gene_id.y")]
# colnames(transcript_tQTL) = c('chr1',	'start1',	'end1',	'gene_name')

labels = do.call(rbind, list(gene_label, CPG_label)) #transcript_tQTL,

labels = labels[!duplicated(labels$label),]
labels$chr = paste0('chr',labels$chr)
dim(labels)
length(labels$label)
dim(table(labels$chr))
dim(table(labels$chr))

# add color
library(RColorBrewer)
# listOfcolors = row.names(brewer.pal.info[which(brewer.pal.info$colorblind),])
# brewer.pal(3,listOfcolors[1])[1]

listOfcolors = brewer.pal(n = 8, name = "Dark2")
i=it=1
l=length(listOfcolors)

Genes = unique(gene_label$label)

for (i in 1:length(Genes)) {
  
  GL = Genes[i] # 'MRPS7'
  CpGL = unique(eQTL_mQTL_corr_CS_coloc$cpg_ID[which(eQTL_mQTL_corr_CS_coloc$gene_name == gene_label$label[i])])
  # print(length(GL)) # all equal to 1
  
  labels[which(labels$label %in% c(GL, CpGL)),'color'] = c(listOfcolors[it], listOfcolors[it])
  
  # I use each color multiple times
  if(it==length(listOfcolors))
  {it=1}else{it=it+1}
  
  print(paste(i, it, sep = '_'))
  
}

table(labels$end>labels$start)

labels$end = labels$start+1

colnames(labels)
dim(labels)

write.csv(labels,paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/FIGURES/Figure-3/supportive/label_data.csv'), row.names = F)
# (labels[,-which(colnames(labels) == c("color"))]

#--------------------------------------  tracks
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/Figures_draft')

#--- add tQTLs
tQTL_LPS24 = fread('/Users/isarnassiri/Documents/RESULTS_USED/tQTL/tQTL_LPS24_annotated_clumped_leadSNPs.txt', stringsAsFactors = F)
tQTL_LPS24_sub = tQTL_LPS24[tQTL_LPS24$FDR < 0.001,]
dim(tQTL_LPS24)

tQTL_LPS24$index_tQTL = paste(tQTL_LPS24$SNP_ID, tQTL_LPS24$gene_name, sep = '_')
eQTL_mQTL_corr_CS_coloc$index_tQTL = paste(eQTL_mQTL_corr_CS_coloc$SNP_ID, eQTL_mQTL_corr_CS_coloc$gene_name, sep = '_')

dim(eQTL_mQTL_corr_CS_coloc)
dim(tQTL_LPS24)

library(dplyr)
eQTL_mQTL_corr_CS_coloc_tQTL = left_join(eQTL_mQTL_corr_CS_coloc, tQTL_LPS24[,c("gene_id", "FDR", "seqnames", "start", "end", "index_tQTL")], by = 'index_tQTL')
dim(eQTL_mQTL_corr_CS_coloc_tQTL)

tQTL_LPS24$index = paste(tQTL_LPS24$SNP_ID, tQTL_LPS24$gene_name, sep = '_')
eQTL_mQTL_corr_CS_coloc$index = paste(eQTL_mQTL_corr_CS_coloc$SNP_ID, eQTL_mQTL_corr_CS_coloc$gene_name, sep = '_')

length(intersect(tQTL_LPS24$SNP_ID, eQTL_mQTL_corr_CS_coloc$SNP_ID))

library(dplyr)
eQTL_mQTL_corr_CS_coloc_tQTL = left_join(eQTL_mQTL_corr_CS_coloc, tQTL_LPS24[,c("gene_id", "FDR", "seqnames", "start", "end", "index")], by = 'index')
dim(eQTL_mQTL_corr_CS_coloc_tQTL)
dim(eQTL_mQTL_corr_CS_coloc)
colnames(eQTL_mQTL_corr_CS_coloc_tQTL)

test = eQTL_mQTL_corr_CS_coloc_tQTL[duplicated(eQTL_mQTL_corr_CS_coloc_tQTL$gene_name),]

#===  add correlation coefficient
library(data.table)
Expression <- fread(paste0('/Users/isarnassiri/Documents/RESULTS_USED/Expression/Input_files_gene/expression_',treatment,'.txt'), stringsAsFactors = F, header = T)

# gene_name_id <- fread('/Users/isarnassiri/Documents/RESULTS_USED/DataBases/Biomart_DB_gene_name_id.txt', stringsAsFactors = F, header = T)
# colnames(gene_name_id)[c(1,2)] = c('symbol', 'id')
# Expression = merge(gene_name_id, expmatrix, by = 'id')

Methylation <- fread(paste0("/Users/isarnassiri/Documents/RESULTS_USED/methylationQTL/QTLtools_PC_Selection/Methylation_profile_LPS.bed.gz"), stringsAsFactors = F, header = T)

colnames(Methylation)
colnames(Expression)

Methylation = as.data.frame(Methylation)
row.names(Methylation) = make.names(Methylation$pid, unique = T)

Expression = as.data.frame(Expression)
row.names(Expression) = make.names(Expression$id, unique = T)

Expression = Expression[,which(colnames(Expression) %in% colnames(Methylation))]
Methylation = Methylation[,which(colnames(Methylation) %in% colnames(Expression))]

identical(colnames(Expression), colnames(Methylation))

for(i in 1:dim(eQTL_mQTL_corr_CS_coloc)[1])
{
 
  Methylation_sub = Methylation[which(row.names(Methylation) == eQTL_mQTL_corr_CS_coloc$cpg_ID[i]),]
  Expression_sub = Expression[which(row.names(Expression) == eQTL_mQTL_corr_CS_coloc$gene_id[i]),]
  eQTL_mQTL_corr_CS_coloc$CC[i] = cor(as.numeric(Expression_sub), as.numeric(Methylation_sub))
  
  if(is.na(eQTL_mQTL_corr_CS_coloc$CC[i])){print(i)}
  
}

range(eQTL_mQTL_corr_CS_coloc_tQTL$FDR[!is.na(eQTL_mQTL_corr_CS_coloc_tQTL$FDR)])
range(eQTL_mQTL_corr_CS_coloc_tQTL$FDR_eQTL[!is.na(eQTL_mQTL_corr_CS_coloc_tQTL$FDR_eQTL)])
range(eQTL_mQTL_corr_CS_coloc_tQTL$FDR_mQTL[!is.na(eQTL_mQTL_corr_CS_coloc_tQTL$FDR_mQTL)])
range(eQTL_mQTL_corr_CS_coloc$CC)

#--- point_tQTL  
table(is.na(eQTL_mQTL_corr_CS_coloc_tQTL$FDR))

data_point=data.frame(chr=paste0('chr',eQTL_mQTL_corr_CS_coloc_tQTL$seqnames.y),
                      start=abs(eQTL_mQTL_corr_CS_coloc_tQTL$start),
                      end= abs(eQTL_mQTL_corr_CS_coloc_tQTL$start+1),
                      value1=-log10(eQTL_mQTL_corr_CS_coloc_tQTL$FDR))

data_point = data_point[!is.na(data_point$value1),]
range(data_point$value1)
dim(data_point)
summary(data_point$end-data_point$start)

write.csv(data_point,paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/FIGURES/Figure-3/supportive/point_tQTL.csv'), row.names = F)

#--- point_mQTL  
data_point=data.frame(chr=paste0('chr',eQTL_mQTL_corr_CS_coloc$seqnames),
                      start=abs(eQTL_mQTL_corr_CS_coloc$start_mQTL),
                      end= abs(eQTL_mQTL_corr_CS_coloc$start_mQTL+1),
                      value1=-log10(eQTL_mQTL_corr_CS_coloc$FDR_mQTL))
data_point = data_point[!is.na(data_point$value1),]
range(data_point$value1)
dim(data_point)
summary(data_point$end-data_point$start)

write.csv(data_point,paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/FIGURES/Figure-3/supportive/point_metQTL.csv'), row.names = F)

#--- point_gQTL  
data_point=data.frame(chr=paste0('chr',eQTL_mQTL_corr_CS_coloc$seqnames),
                      start=abs(eQTL_mQTL_corr_CS_coloc$start_gQTL),
                      end= abs(eQTL_mQTL_corr_CS_coloc$start_gQTL+1),
                      value1=-log10(eQTL_mQTL_corr_CS_coloc$FDR_eQTL))
data_point = data_point[!is.na(data_point$value1),]
range(data_point$value1)
dim(data_point)
write.csv(data_point,paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/FIGURES/Figure-3/supportive/point_eQTL.csv'), row.names = F)

data_point=data.frame(chr=paste0('chr',eQTL_mQTL_corr_CS_coloc$seqnames),
                      start=abs(eQTL_mQTL_corr_CS_coloc$start_gQTL),
                      end= abs(eQTL_mQTL_corr_CS_coloc$start_gQTL+1),
                      value1=eQTL_mQTL_corr_CS_coloc$CC)
range(eQTL_mQTL_corr_CS_coloc$CC)
write.csv(data_point,paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/FIGURES/Figure-3/supportive/point_correlation.csv'), row.names = F)

#---- CytoBand
library(RCircos)
data("UCSC.HG38.Human.CytoBandIdeogram")
chro.bands = UCSC.HG38.Human.CytoBandIdeogram
colnames(chro.bands) = c('chr',	'start',	'end',	'value1',	'value2')
chro.bands$chr = as.character(chro.bands$chr)
# chro.bands = chro.bands[-which(chro.bands$chr %in% c('chrX')),]
table(chro.bands$chr)

chro.bands$start[which(chro.bands$start == 0)] = 1

write.csv(chro.bands,'/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/FIGURES/Figure-3/supportive/Hg38_chromosome_cytoband.csv',row.names = F)

#-- parameters
#99CC33
#FFCC33
#CCCC99

