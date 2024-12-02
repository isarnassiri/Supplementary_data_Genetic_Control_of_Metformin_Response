############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
module add python 
module add HTSeq/0.6.1p1

#----------------- Description -----------------
# By quantifying the number of reads mapping to each exon, DEXSeq_count enables the detection of subtle changes in splicing patterns that may not be apparent at the gene level.

#----------------- Output -----------------
# The results are presented in text files including exon expression level.

#----------------- Examples -----------------
# You need to specify a path to the output of the bam and reference genome files.

##########################################################################################################################################################

python ~/dexseq_count.py -f bam -s no -r pos ~/EXON.gff INPUT.bam OUTPUT_exon_read_count.txt
