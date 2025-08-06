


#---------------------- Compare the count of mitocondrias with age and sex ----------------------

# Reviewer's comment
# 7. The manuscript would also benefit from an expanded discussion on the potential confounding effects of mtDNA count variability due to age, disease states, and oxidative stress. Specifically, the implications of these factors on the study’s findings, such as the association between rs2910792 and ERAP2, should be addressed. Furthermore, the authors should critically evaluate the representativeness of their sample set concerning the patient cohorts in the cited GWA studies (such as IBD) where rs2910792 co-loc was detected.

#---------------------- Figure 1. Association_between_mitochondrial_count_and_age

setwd('./') # The location of the Mitochondria_count folder

library(data.table)
mit_count = fread('IFN.gamma.Mitochondria_count.txt', stringsAsFactors = F, header = T, data.table = F)

mit_count = mit_count[order(mit_count$plt.norm),]
mit_count = mit_count[!duplicated(mit_count$ID),]
mit_count = mit_count[order(mit_count$ID),]

Metadata_array = fread('Age_sex_BMI_432_information.txt', stringsAsFactors = F, header = T, data.table = F)

colnames(Metadata_array)[1] = "ID"

mit_count_aanotated = merge(Metadata_array, mit_count, by = 'ID')
dim(mit_count_aanotated)

Age <- mit_count_aanotated$Age
Mitochondria_count <- mit_count_aanotated$mean.ratio

pdf(file = "Association_between_mitochondrial_count_and_age.pdf",    
    width = 5,                  
    height = 5)                 

plot(Age, Mitochondria_count, pch = 19, col = "lightblue")
abline(lm(Mitochondria_count ~ Age), col = "red", lwd = 3)
PC = cor.test(Age, Mitochondria_count)
text(paste("Correlation:", round(as.numeric(PC$estimate), digits = 2)), x = 35, y = 95)

dev.off()

#---------------------- Figure 2. Association_between_mitochondrial_count_and_sex

mit_count_aanotated$Type[which(mit_count_aanotated$Sex == 1)] = 'Male' # one X chromosome
mit_count_aanotated$Type[which(mit_count_aanotated$Sex == 2)] = 'Female' # two X chromosome

library(ggpubr)
p <- ggboxplot(mit_count_aanotated, x = "Type", y = 'mean.ratio', short.panel.labs = FALSE)

pdf('Association_between_mitochondrial_count_and_sex.pdf', width = 8, height = 6, useDingbats = F)
p + stat_compare_means(label =  "p.signif", label.x = 1.5) + theme(text = element_text(size = 20))
dev.off()
