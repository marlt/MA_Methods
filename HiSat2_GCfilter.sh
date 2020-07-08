#!/bin/bash
# Spaghetti code for the HiSat2 remapping against contigs of GC filtered read sets

# Paths
Rbase='/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades/'
uREADS='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/'

# hisat2
# samtools
mkdir -p ${Rbase}/GCfilter_test/bigger25/remapping 
mkdir -p ${Rbase}/GCfilter_test/bigger45/remapping 
mkdir -p ${Rbase}/GCfilter_test/25_45/remapping 
mkdir -p ${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping 
mkdir -p ${Rbase}/GCfilter_test/all/remapping 

# INDEX
echo 'Create Index for bigger25'
nice hisat2-build -p 40 ${Rbase}/GCfilter_test/bigger25/contigs.fasta ${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25 &> ${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_indexing.log
echo 'Create Index for bigger45'
nice hisat2-build -p 40 ${Rbase}/GCfilter_test/bigger45/contigs.fasta ${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45 &> ${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_indexing.log
echo 'Create Index 25_45'
nice hisat2-build -p 40 ${Rbase}/GCfilter_test/25_45/contigs.fasta ${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45 &> ${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45_indexing.log
echo 'Create Index merged_25_45_+_bigger45'
nice hisat2-build -p 40 ${Rbase}/GCfilter_test/merged_25_45_+_bigger45/contigs.fasta ${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged &> ${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_indexing.log
echo 'Create Index all'
nice hisat2-build -p 40 ${Rbase}/GCfilter_test/all/contigs.fasta ${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all &> ${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_indexing.log


# REMAPPING
echo 'Remapping for bigger25'
nice hisat2 -p 40 -k 15 -x "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25" -2 "${Rbase}/GCfilter_test/bigger25/corrected/all_viral_pair1_trimmed_filterGC.00.0_0.cor.fastq"  -1 "${Rbase}/GCfilter_test/bigger25/corrected/all_viral_pair2_trimmed_filterGC.00.0_0.cor.fastq" -U "${Rbase}/GCfilter_test/bigger25/corrected/all_viral_allSE_trimmed_filterGC.00.0_1.cor.fastq"  -S "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25.sam" --un "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_unmapped.fastq" &> "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_remapping.log"
echo 'Remapping for bigger45'
nice hisat2 -p 40 -k 15 -x "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45" -2 "${Rbase}/GCfilter_test/bigger45/corrected/all_viral_pair1_trimmed_filterGC_45.00.0_0.cor.fastq"  -1 "${Rbase}/GCfilter_test/bigger45/corrected/all_viral_pair2_trimmed_filterGC_45.00.0_0.cor.fastq" -U "${Rbase}/GCfilter_test/bigger45/corrected/all_viral_allSE_trimmed_filterGC_45.00.0_1.cor.fastq"  -S "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45.sam" --un "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_unmapped.fastq" &> "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_remapping.log"
echo 'Remapping for 25_45'
nice hisat2 -p 40 -k 15 -x "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45" -2 "${Rbase}/GCfilter_test/25_45/corrected/all_viral_pair1_trimmed_filterGC_25_45.00.0_0.cor.fastq"  -1 "${Rbase}/GCfilter_test/25_45/corrected/all_viral_pair2_trimmed_filterGC_25_45.00.0_0.cor.fastq" -U "${Rbase}/GCfilter_test/25_45/corrected/all_viral_allSE_trimmed_filterGC_25_45.00.0_1.cor.fastq"  -S "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45.sam" --un "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45_unmapped.fastq" &> "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45_remapping.log"
echo 'Remapping for merged_25_45_+_bigger45'
nice hisat2 -p 40 -k 15 -x "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged" -2 "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/corrected/all_viral_pair1_trimmed_filterGC_merged.00.0_0.cor.fastq" -1 "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/corrected//all_viral_pair2_trimmed_filterGC_merged.00.0_0.cor.fastq" -U "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/corrected/all_viral_allSE_trimmed_filterGC_merged.00.0_1.cor.fastq"  -S "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged.sam" --un "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_unmapped.fastq" &> "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_remapping.log"
echo 'Remapping for all contigs'
nice hisat2 -p 40 -k 15 -x "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all" -2 "${Rbase}/GCfilter_test/all/corrected/all_viral_pair1_trimmed.00.0_0.cor.fastq"  -1 "${Rbase}/GCfilter_test/all/corrected/all_viral_pair2_trimmed.00.0_0.cor.fastq" -U "${Rbase}/GCfilter_test/all/corrected/all_viral_allSE_trimmed.00.0_1.cor.fastq" -S "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all.sam" --un "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_unmapped.fastq" &> "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_remapping.log"

# Convert sam to bam
echo "Make bam..."
samtools view -b "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25.sam" > "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25.bam"
samtools view -b "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45.sam" > "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45.bam"
samtools view -b "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_25_45.sam" > "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45.bam"
samtools view -b "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged.sam" > "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged.bam"
samtools view -b "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all.sam" > "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all.bam"

# Sort the bam files
echo "Sort bam..."
samtools sort "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25.bam" > "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_sorted.bam"
samtools sort "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45.bam" > "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_sorted.bam"
samtools sort "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45.bam" > "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45_sorted.bam"
samtools sort "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged.bam" > "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_sorted.bam"
samtools sort "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all.bam" > "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_sorted.bam"

# Index the sorted bams
echo "Index sorted bam..."
samtools index "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_sorted.bam"
samtools index "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_sorted.bam"
samtools index "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45_sorted.bam"
samtools index "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_sorted.bam"
samtools index "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_sorted.bam"

# Get the statistics (contig lengths and average depth)
echo "Output statistics..."
samtools idxstats "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_sorted.bam" > "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_sorted.stat"
samtools idxstats "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_sorted.bam" > "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_sorted.stat"
samtools idxstats "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45_sorted.bam" > "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45_sorted.stat"
samtools idxstats "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_sorted.bam" > "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_sorted.stat"
samtools idxstats "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_sorted.bam" > "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_sorted.stat"

# Compute per base depth
echo "Output depth..."
samtools depth -a "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_sorted.bam" > "${Rbase}/GCfilter_test/bigger25/remapping/GCfilter_test_bigger25_sorted.depth.tsv"
samtools depth -a "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_sorted.bam" > "${Rbase}/GCfilter_test/bigger45/remapping/GCfilter_test_bigger45_sorted.depth.tsv"
samtools depth -a "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45_sorted.bam" > "${Rbase}/GCfilter_test/25_45/remapping/GCfilter_test_bigger25_45_sorted.depth.tsv"
samtools depth -a "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_sorted.bam"  > "${Rbase}/GCfilter_test/merged_25_45_+_bigger45/remapping/GCfilter_test_merged_sorted.depth.tsv" 
samtools depth -a "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_sorted.bam"  > "${Rbase}/GCfilter_test/all/remapping/GCfilter_test_all_sorted.depth.tsv" 

# hisat2 Parameters:
#   -p  threads to use
#   -k  Primary alignments to report; the best alignment will be chosen from them
#   -x  specify the path and base name to an reference seq index

# samtools Parameters:
#   -b  input is BAM
#   -a  output all positions in contig, even if depth is