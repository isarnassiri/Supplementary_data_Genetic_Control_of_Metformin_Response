
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/KEY-FILES/Articles_Monocyte_eQTL/Article_ciseQTLs/eQTL_mQTL_pairs_analysis/')

#==== read in
treatment='UT'

#=== met-exp correlation
library(data.table)
methExprs = fread(paste0(treatment,'correlationMethExprs_results.txt'), stringsAsFactors = F, header = T)
methExprs_sub = methExprs[which(methExprs$adj.P.Val<0.05),]
#  beta: coefficient of the methylation change

#== metQTL all
metQTL = fread(paste0('/Users/isarnassiri/Documents/RESULTS_USED/methylationQTL/QTLtools_PC_Selection/nominal_imputedGenotype_all_',treatment,'_annotated.txt'), stringsAsFactors = F)
dim(metQTL)
colnames(metQTL)[1] = c('cpg_ID')

metQTL_sub = metQTL[metQTL$FDR < 0.001,]
dim(metQTL_sub)
head(metQTL_sub)
metQTL_sub$SNP_ID = gsub('.*;', '', metQTL_sub$SNP_ID)
table(metQTL_sub$chrID)

#== eQTL peak LPS
eQTL = fread(paste0('/Users/isarnassiri/Documents/RESULTS_USED/eQTL/eQTL_',treatment,'_annotated_clumped_leadSNPs.txt'), stringsAsFactors = F)
eQTL_sub = eQTL[eQTL$FDR < 0.001,]
dim(eQTL_sub)

#== eQTL_mQTL_pairs_all
colnames(eQTL_sub)[match(c("FDR", "slop"), colnames(eQTL_sub))] = c("FDR_eQTL", "slope_eQTL")
colnames(metQTL_sub)[match(c("FDR", "slope"), colnames(metQTL_sub))] = c("FDR_mQTL", "slope_mQTL")

eQTL_mQTL_all = merge(eQTL_sub, metQTL_sub, by = 'SNP_ID')
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

eQTL_mQTL_corr_CS = eQTL_mQTL_corr # correlation in UT is weak I omitted the eQTL_mQTL_pairs_correlated + no corr in LPS part

#=============================== eQTL_mQTL_pairs_correlated + no corr in LPS
# 
# met_exp_LPS = fread(paste0('LPS24correlationMethExprs_results.txt'), stringsAsFactors = F ,header = T)
# range(met_exp_LPS$adj.P.Val)
# met_exp_LPS = met_exp_LPS[met_exp_LPS$adj.P.Val < 0.001,]
# 
# met_exp_LPS$index = paste(met_exp_LPS$cpg, met_exp_LPS$exprs, sep = '_')
# 
# eQTL_mQTL_corr_CS = eQTL_mQTL_corr[-which(eQTL_mQTL_corr$index %in% met_exp_LPS$index),]
# dim(eQTL_mQTL_corr_CS)
# 
# library(data.table)
# fwrite(eQTL_mQTL_corr_CS, 'eQTL_mQTL_corr_CS_coloc_UT.txt', quote = F, row.names = F, sep = '\t')
# 
# eQTL_mQTL_corr_CS[(eQTL_mQTL_corr_CS$exprs == 'CD55'),]
# eQTL_mQTL_corr_CS[grep('CD55', eQTL_mQTL_corr_CS$gene_name),]

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


(RESULTS[which(RESULTS$gene_name == 'CD55'),])


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
table(eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8)
eQTL_mQTL_corr_CS_coloc = eQTL_mQTL_corr_CS_coloc[which(eQTL_mQTL_corr_CS_coloc$PP.H4.abf>0.8),]
colnames(eQTL_mQTL_corr_CS_coloc)

#------------ add se ------------
#== metQTL nominal

library(data.table)
header = fread('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/gQTL_nominal_Header.txt', stringsAsFactors = F, header = F)

#-- nominal
Gene_nominal = fread(paste0('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/gQTL/QTLtools_outputs/',treatment,'/nominal_pass/gQTL_nominal_pass_1.txt'), stringsAsFactors = F)
colnames(Gene_nominal) = as.character(header)
Gene_nominal = Gene_nominal[which(Gene_nominal$var_id %in% eQTL_mQTL_corr_CS_coloc$SNP_ID),]
Gene_nominal$index = paste(Gene_nominal$var_id, Gene_nominal$phe_id, sep = '_')

eQTL_mQTL_corr_CS_coloc$index = paste(eQTL_mQTL_corr_CS_coloc$SNP_ID, eQTL_mQTL_corr_CS_coloc$gene_name, sep = '_')
eQTL_mQTL_corr_CS_coloc_interset = merge(eQTL_mQTL_corr_CS_coloc, Gene_nominal[,c('slope_se', 'index')], by = 'index')
dim(eQTL_mQTL_corr_CS_coloc_interset)

#-- I do not have se for all gQTLs, therefore I use mean
eQTL_mQTL_corr_CS_coloc$slope_se = summary(eQTL_mQTL_corr_CS_coloc_interset$slope_se)[4]

identical(colnames(eQTL_mQTL_corr_CS_coloc_interset), colnames(eQTL_mQTL_corr_CS_coloc))
class(eQTL_mQTL_corr_CS_coloc_interset)
class(eQTL_mQTL_corr_CS_coloc)

eQTL_mQTL_corr_CS_coloc_interset = as.data.frame(eQTL_mQTL_corr_CS_coloc_interset)
eQTL_mQTL_corr_CS_coloc = as.data.frame(eQTL_mQTL_corr_CS_coloc)

eQTL_mQTL_corr_CS_coloc = rbind(eQTL_mQTL_corr_CS_coloc_interset, eQTL_mQTL_corr_CS_coloc[-which(eQTL_mQTL_corr_CS_coloc$index %in% eQTL_mQTL_corr_CS_coloc_interset$index),])
dim(eQTL_mQTL_corr_CS_coloc)

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
table(qtl$dir_p<0.05) # steiger is equal to dir_p
# c('FALSE', 'TRUE'))] = c("Expression causes Methylation", "Methylation causes Expression")
table(qtl$mr_p<0.05) 
# qtl = qtl[which(qtl$mr_p<0.05),]
dim(qtl)

(qtl[which(qtl$SNP_ID == 'rs2914937'),])

library(data.table)
fwrite(qtl, paste0('eQTL_mQTL_corr_CS_mr_steiger_',treatment,'.txt'), quote = F, row.names = F, sep = '\t')
getwd()

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
mosaicplot(input2, shade=F, color=F, las=1, type="pearson", cex.axis=1, main="")

#-- proportion for table
sum(input2)
round(input2/sum(input2), digits = 3)*100

## Loading required package: grid
setwd('/Users/isarnassiri/Documents/Analysis_FairfaxLab/New_Analysis_eQTL_Monocyte/Figures_draft')
pdf(paste0('mosaicplot_revised',treatment,'.pdf'), width = 10, height = 8, useDingbats = F)
mosaicplot(input2, shade=F, color=F, las=1, type="pearson", cex.axis=1, main="")
dev.off()
getwd()
#
 
