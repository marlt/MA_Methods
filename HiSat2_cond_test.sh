#!/bin/bash
## REMAPPING ON TEST CONDITION ASSEMBLIES
# Spaghetti Code for remapping reads to the condition test assemblies, converting into bam and writing statistics as well as base depth reports

# Paths
Rbase='/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades/cond_test'
uREADS='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/'

# Create folders
mkdir -p ${Rbase}/auto_kmer/corrected/remapping
mkdir -p ${Rbase}/auto_kmer/uncorrected/remapping
mkdir -p ${Rbase}/set_kmer/corrected/remapping
mkdir -p ${Rbase}/set_kmer/uncorrected/remapping

# INDEXING
echo 'Create Index auto corrected'
nice hisat2-build -p 40 ${Rbase}/auto_kmer/corrected/contigs.fasta ${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test &> ${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_indexing.log
echo 'Create Index auto uncorrected'
nice hisat2-build -p 40 ${Rbase}/auto_kmer/uncorrected/contigs.fasta ${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test &> ${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_indexing.log
echo 'Create Index set corrected'
nice hisat2-build -p 40 ${Rbase}/set_kmer/corrected/contigs.fasta ${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test &> ${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_indexing.log
echo 'Create Index set uncorrected'
nice hisat2-build -p 40 ${Rbase}/set_kmer/uncorrected/contigs.fasta ${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test &> ${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_indexing.log

# REMAPPING
echo 'Remap auto corrected...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test -U "${Rbase}/auto_kmer/corrected/corrected/all_viral_pooled_trimmed_filterGC.00.0_0.cor.fastq" -S "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test.sam" --un "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_unmapped.fastq" &> "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_remapping.log"
echo 'Remap auto uncorrected...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test -U "${uREADS}/all_viral_pooled_trimmed_filterGC.fastq" -S "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test.sam" --un "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_unmapped.fastq" &> "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_remapping.log"
echo 'Remap set corrected...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test -U "${Rbase}/set_kmer/corrected/corrected/all_viral_pooled_trimmed_filterGC.00.0_0.cor.fastq" -S "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test.sam" --un "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_unmapped.fastq" &> "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_remapping.log"
echo 'Remap set uncorrected...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test -U "${uREADS}/all_viral_pooled_trimmed_filterGC.fastq" -S "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test.sam" --un "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_unmapped.fastq" &> "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_remapping.log"

# Convert sam to bam
echo "Make bam..."
samtools view -b "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test.sam" > "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test.bam"
samtools view -b "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test.sam" > "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test.bam"
samtools view -b "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test.sam" > "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test.bam"
samtools view -b "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test.sam" > "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test.bam"

# Sort the bam files
echo "Sort bam..."
samtools sort "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test.bam" > "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_sorted.bam"
samtools sort "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test.bam" > "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_sorted.bam"
samtools sort "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test.bam" > "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_sorted.bam"
samtools sort "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test.bam" > "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_sorted.bam"

# Index the sorted bams
echo "Index sorted bam..."
samtools index "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_sorted.bam"
samtools index "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_sorted.bam"
samtools index "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_sorted.bam"
samtools index "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_sorted.bam"

# Get the statistics (contig lengths and average depth)
echo "Output statistics... "
samtools idxstats "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_sorted.bam" > "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_sorted.stat"
samtools idxstats "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_sorted.bam" > "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_sorted.stat"
samtools idxstats "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_sorted.bam" > "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_sorted.stat"
samtools idxstats "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_sorted.bam" > "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_sorted.stat"

# Compute per base depth
echo "Compute cov depth..."
samtools depth -a "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test_sorted.bam"  > "${Rbase}/auto_kmer/corrected/remapping/auto_cor_cond_test.depth.tsv"
samtools depth -a "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test_sorted.bam" > "${Rbase}/auto_kmer/uncorrected/remapping/auto_uncor_cond_test.depth.tsv"
samtools depth -a "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test_sorted.bam" > "${Rbase}/set_kmer/corrected/remapping/set_cor_cond_test.depth.tsv"
samtools depth -a "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test_sorted.bam" > "${Rbase}/set_kmer/uncorrected/remapping/set_uncor_cond_test.depth.tsv"


# hisat2 Parameters:
#   -p  threads to use
#   -k  Primary alignments to report; the best alignment will be chosen from them
#   -x  specify the path and base name to an reference seq index

# samtools Parameters:
#   -b  input is BAM
#   -a  output all positions in contig, even if depth is 0