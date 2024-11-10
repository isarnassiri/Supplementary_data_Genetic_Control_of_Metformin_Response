############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
module add stringtie
module add hisat

#----------------- Description -----------------
# StringTie can merge transcripts from multiple samples, creating a unified set of transcripts for comparative analysis. This allows for consistent quantification and identification of differentially expressed transcripts across samples.

#----------------- Output -----------------
# The script generates a count table, providing the number of reads mapped to each transcript.

#----------------- Examples -----------------
# Execute the command below with your list of samples and gene annotation file (gtf) as inputs.

##########################################################################################################################################################

stringtie --merge -p 8 -G ~/GRCh38.gtf -o ~/stringtie_merged.gtf ~/List_of_GTFs_per_Sample.txt # OUTPUT_nodup_properPairs_NH.gtf
