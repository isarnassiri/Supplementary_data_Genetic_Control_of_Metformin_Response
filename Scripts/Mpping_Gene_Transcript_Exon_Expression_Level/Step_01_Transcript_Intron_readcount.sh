############################## Code Description ##########################################################################################################

#----------------- Packages that need to be installed -----------------
module add hisat2 
module add samtools
module add picard-tools
module add stringtie
module add bamtools

#----------------- Description -----------------
# StringTie is a tool for reconstructing and quantifying transcript abundance from RNA-seq data. 
# Step-01: StringTie leverages a reference genome and annotation to assemble transcripts from RNA-seq reads, identifying isoforms and their expression levels.
# Step-02: StringTie calculates Fragments Per Kilobase of transcript per Million mapped reads (FPKM) and Transcripts Per Million (TPM) to estimate transcript abundance, normalizing for transcript length and sequencing depth.
# Step-03: StringTie generates input files for Ballgown, enabling downstream expression analysis.

#----------------- Output -----------------
# The script generates a count table, providing the number of reads mapped to each transcript.

#----------------- Examples -----------------
# Use a suitable aligner (e.g., HISAT2) to map sequencing reads to a reference genome.
# Create a GTF or GFF file detailing transcript structures.
# Execute the command below with your aligned reads (bam files) and gene annotation file (gtf) as inputs.

##########################################################################################################################################################


# Use a HISAT2 aligner to map sequencing reads to a reference genome.
hisat2 -x ~/GRCh38.dna -p 8 -1 INPUT_1.fastq.gz -2 INPUT_2.fastq.gz --rg-id=NAME --rg SM: OUTPUT | samtools view -bS - > OUTPUT.bam

# sort bam files
samtools sort --threads 8 OUTPUT.bam > OUTPUT_sorted.bam 

# identify duplicated reads as sequencing artifacts
MarkDuplicates I=OUTPUT_sorted.bam O=OUTPUT_nodup.bam M=OUTPUT_marked_dup_metrics.txt REMOVE_DUPLICATES=true

# filter out duplicated reads
bamtools filter -tag NH:1 -in OUTPUT_nodup.bam | bamtools filter -isProperPair true -out OUTPUT_nodup_properPairs_NH.bam

# generate index file for bam file
BuildBamIndex I=OUTPUT_nodup_properPairs_NH.bam; 

# generate summary QC for the bam file
samtools flagstat OUTPUT_nodup_properPairs_NH.bam > OUTPUT_nodup_properPairs_NH_flagstat.txt;

# The script generates a count table, providing the number of reads mapped to each transcript.
stringtie -p 8 -G ~/GRCh38.gtf -o OUTPUT_nodup_properPairs_NH.gtf OUTPUT_nodup_properPairs_NH.bam;
