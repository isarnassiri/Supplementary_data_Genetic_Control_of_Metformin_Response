############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
library(methods)
library(ballgown)
library(rtracklayer)

#----------------- Description -----------------
# Ballgown is a powerful tool for analyzing RNA-seq data, including transcript-level expression quantification. 
# Step-1: Import and organize RNA-seq data from StringTie.
# Step-2: Calculate transcript-level expression values, such as FPKM (Fragments Per Kilobase of transcript per Million mapped reads).

#----------------- Output -----------------
# The results are presented in text files including transcript, gene, exon, and intron expression level.

#----------------- Examples -----------------

# You need to specify a path to the output of the stringtie tool including the following files per sample:
# e2t.ctab	e_data.ctab	i2t.ctab	i_data.ctab	t_data.ctab

# Additionally, you need to provide a path to stringtie_merged.gtf, which represents the output of stringtie --merge for all samples.

##########################################################################################################################################################

# You need to specify a path to the output of the stringtie tool including the following files per sample:
# e2t.ctab	e_data.ctab	i2t.ctab	i_data.ctab	t_data.ctab
bg = ballgown(dataDir='~/', samplePattern='-', meas='all')

# You need to provide a path to stringtie_merged.gtf, which represents the output of stringtie --merge for all samples.
gtf <- rtracklayer::import('~/stringtie_merged.gtf')
gtf_df=as.data.frame(gtf)

#--- Transcript read count
rm(transcript_fpkm)
transcript_fpkm = texpr(bg, 'all')
dim(transcript_fpkm)
transcript_fpkm <- transcript_fpkm[,-which(grepl('cov.', colnames(transcript_fpkm), fixed=FALSE))]
dim(transcript_fpkm)

transcript_fpkm <- transcript_fpkm[-which(transcript_fpkm$gene_name == "."),]
dim(transcript_fpkm)

rsum = rowSums(transcript_fpkm[,11:dim(transcript_fpkm)[2]])
table((rsum > 1))
transcript_fpkm = transcript_fpkm[which(rsum >= 1),]
transcript_fpkm[is.na(transcript_fpkm)] <- 0
dim(transcript_fpkm)

write.table(transcript_fpkm[,c(2,4:6)], "~/GC_content_input_ISOFORMS.bed", row.names= FALSE, quote=FALSE, sep = "\t")
write.table(transcript_fpkm[,c(6,2,4,5)], "~/ISOFORMS_location.txt", row.names= FALSE, quote=FALSE, sep = "\t")

rownames(transcript_fpkm) <- transcript_fpkm$t_name
write.table(transcript_fpkm[,c(1:10)], "~/ISOFORMS_annotation.txt", row.names= FALSE, quote=FALSE, sep = "\t")

transcript_fpkm <- transcript_fpkm[,-c(1:10)]
colnames(transcript_fpkm) <- gsub('FPKM.','',colnames(transcript_fpkm))

#---- Save outputs Transcript read count
write.table(transcript_fpkm, "~/ISOFORMS_merged.txt", row.names= TRUE, quote=FALSE, sep = "\t")

#--- Exon read count 
exon = eexpr(bg, 'mrcount')
whole_exon_table = eexpr(bg, 'all')

exon_transcript_table = indexes(bg)$e2t
exon_transcript_table[,3] <- as.character(transcript_gene_table[as.numeric(exon_transcript_table[,2]),2])
for(i in 1:dim(exon)[1])
{
  rownames(exon)[i] <- paste(unlist(exon_transcript_table[i,]), collapse ="_");
}

#--- Gene read count
gene_expression = gexpr(bg) 

#--- Intron read count
Intron = iexpr(bg, 'mrcount')
whole_Intron_table = iexpr(bg, 'all')

rsum = rowSums(Intron)
Intron = Intron[which(rsum >= 1),]
Intron[is.na(Intron)] <- 0
dim(Intron)

whole_Intron_table = whole_Intron_table[which(rsum >= 1),]
rownames(Intron) <- paste(whole_Intron_table$i_id,'_',whole_Intron_table$chr,'_',whole_Intron_table$strand,'_',whole_Intron_table$start,'_',whole_Intron_table$end, sep='')

#---- Save output count for Intron reads.
write.table(Intron, '~/INTRON_merged.txt', row.names= TRUE, quote=FALSE, sep = "\t")
write.table(cbind(whole_Intron_table[,c(2,4,5)],rownames(Intron)), '~/GC_content_input_INTRON.bed', row.names= FALSE, quote=FALSE, sep = "\t")

location <- cbind(rownames(Intron),whole_Intron_table[,c(2,4,5)])
colnames(location) <- c('id',  'chr',	'start',	'end')
write.table(location, '~/intron_location.txt', row.names= FALSE, quote=FALSE, sep = "\t")

