#!/bin/bash



##################################################################
# Call SNPs from the mapped RNA-seq reads
# following this biostar post: https://www.biostars.org/p/327558/
#################################################################

bcftools mpileup --threads 30 --max-depth 10000 --redo-BAQ --min-BQ 30 --per-sample-mF -f S_lycopersicum_chromosomes.2.50.fa \
--skip-indels --output-type v ../results_lyco_genome/star/F2-28_concat_Aligned.sortedByCoord.out.bam | \
bcftools call --multiallelic-caller --variants-only -Ov - > F2-28_concat.vcf 

# -d INT    minimum read depth [2]
vcfutils.pl varFilter -d 10 F2-28_concat.vcf  > F2-28_concat.filtered.vcf

vcf2bed < F2-28_concat.filtered.vcf > F2-28_concat.filtered.bed

##############################################
# Make windows from the S. lycopersicum genome 
##############################################
samtools faidx S_lycopersicum_chromosomes.2.50.fa
cp S_lycopersicum_chromosomes.2.50.fa.fai chromsizes_SL2.50.txt
awk '{print $1 "\t" $2}' chromsizes_SL2.50.txt > chromsizes.txt

WINDOW_SIZE=1000000 # 1Mb
bedtools makewindows -g chromsizes.txt -w $WINDOW_SIZE > genome.bed


#####################################
# Count number of SNPs per 1Mb window
######################################

bedmap --echo --count --delim "\t" genome.bed F2-28_concat.filtered.bed> F2-28_snp_counts.tsv 


