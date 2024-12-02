############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
module add python-cbrg

#----------------- Description -----------------
# HTSeq.scripts.count tool counts reads mapped to specific genomic features, providing crucial information for downstream analyses like differential expression analysis.

#----------------- Output -----------------
# The script generates a count table, providing the number of reads mapped to each gene.

#----------------- Examples -----------------
# Use a suitable aligner (e.g., HISAT2) to map sequencing reads to a reference genome.
# Create a GTF or GFF file detailing gene structures and exons.
# Execute the command below with your aligned reads (bam files) and gene annotation file (gtf) as inputs.

##########################################################################################################################################################

python -m HTSeq.scripts.count -r pos --format=bam --minaqual=0 --stranded=no --mode=union INPUT.bam ~/GRCh38.gtf > OUTPUT_gene_read_count.txt 
