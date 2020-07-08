#!/bin/bash
# REMAPPING ON TEST CONDITION ASSEMBLIES
# Spaghetti Code for remapping reads to assemblies generated from different input read sets and/or generated with metaspades or spades

# Paths
Rbase='/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades/test_ass_type'
uREADS='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/'

mkdir -p ${Rbase}/spades/allSE/remapping
mkdir -p ${Rbase}/spades/PE/remapping
mkdir -p ${Rbase}/spades/PE_SE/remapping

mkdir -p ${Rbase}/metaspades/PE/remapping
mkdir -p ${Rbase}/metaspades/PE_SE/remapping

# ###### INDEXING
# spades
echo 'Create spades Index allSE'
nice hisat2-build -p 40 ${Rbase}/spades/allSE/contigs.fasta ${Rbase}/spades/allSE/remapping/test_ass_type_allSE &> ${Rbase}/spades/allSE/remapping/test_ass_type_allSE_indexing.log
echo 'Create spades Index PE'
nice hisat2-build -p 40 ${Rbase}/spades/PE/contigs.fasta ${Rbase}/spades/PE/remapping/test_ass_type_PE &> ${Rbase}/spades/PE/remapping/test_ass_type_PE_indexing.log
echo 'Create spades Index PE_SE'
nice hisat2-build -p 40 ${Rbase}/spades/PE_SE/contigs.fasta ${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE &> ${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_indexing.log

# metaspades
echo 'Create metaspades Index PE'
nice hisat2-build -p 40 ${Rbase}/metaspades/PE/contigs.fasta ${Rbase}/metaspades/PE/remapping/test_ass_type_PE &> ${Rbase}/metaspades/PE/remapping/test_ass_type_PE_indexing.log
echo 'Create metaspades Index PE_SE'
nice hisat2-build -p 40 ${Rbase}/metaspades/PE_SE/contigs.fasta ${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE &> ${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_indexing.log

# ####### REMAPPING
echo 'Remap allSE to spades contigs...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/spades/allSE/remapping/test_ass_type_allSE -U "${Rbase}/spades/allSE/corrected/all_viral_allSE_trimmed_filterGC.00.0_0.cor.fastq" -S "${Rbase}/spades/allSE/remapping/test_ass_type_allSE.sam" --un "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_unmapped.fastq" &> "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_remapping.log"
echo 'Remap PE to spades contigs...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/spades/PE/remapping/test_ass_type_PE -2 "${Rbase}/spades/PE/corrected/all_viral_pair1_trimmed_filterGC.00.0_0.cor.fastq" -1 "${Rbase}/spades/PE/corrected/all_viral_pair2_trimmed_filterGC.00.0_0.cor.fastq" -U "${Rbase}/spades/PE/corrected/all_viral_pair_unpaired.00.0_0.cor.fastq" -S "${Rbase}/spades/PE/remapping/test_ass_type_PE.sam" --un "${Rbase}/spades/PE/remapping/test_ass_type_PE_unmapped.fastq" &> "${Rbase}/spades/PE/remapping/test_ass_type_PE_remapping.log"
echo 'Remap PE_SE t o spades contigs...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE -2 "${Rbase}/spades/PE_SE/corrected/all_viral_pair1_trimmed_filterGC.00.0_0.cor.fastq" -1 "${Rbase}/spades/PE_SE/corrected/all_viral_pair2_trimmed_filterGC.00.0_0.cor.fastq" -U "${Rbase}/spades/PE_SE/corrected/all_viral_allSE_trimmed_filterGC.00.0_1.cor.fastq" -S "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE.sam" --un "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_unmapped.fastq" &> "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_remapping.log"

echo 'Remap PE to metaspades contigs...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/metaspades/PE/remapping/test_ass_type_PE -2 "${Rbase}/metaspades/PE/corrected/all_viral_pair1_trimmed_filterGC.00.0_0.cor.fastq" -1 "${Rbase}/metaspades/PE/corrected/all_viral_pair2_trimmed_filterGC.00.0_0.cor.fastq" -U "${Rbase}/metaspades/PE/corrected/all_viral_allSE.00.0_1.cor.fastq" -S "${Rbase}/metaspades/PE/remapping/test_ass_type_PE.sam" --un "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_unmapped.fastq" &> "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_remapping.log"
echo 'Remap PE + SE to metaspades contigs...'
nice hisat2 -p 40 -k 15 -x ${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE -2 "${Rbase}/metaspades/PE_SE/corrected/all_viral_pair1_trimmed_filterGC.00.0_0.cor.fastq" -1 "${Rbase}/metaspades/PE_SE/corrected/all_viral_pair2_trimmed_filterGC.00.0_0.cor.fastq" -U "${Rbase}/metaspades/PE_SE/corrected/all_viral_allSE_trimmed_filterGC.00.0_1.cor.fastq" -S "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE.sam" --un "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_unmapped.fastq" &> "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_remapping.log"


# Convert sam to bam
echo "Make bam..."
samtools view -b "${Rbase}/spades/allSE/remapping/test_ass_type_allSE.sam" > "${Rbase}/spades/allSE/remapping/test_ass_type_allSE.bam"
samtools view -b "${Rbase}/spades/PE/remapping/test_ass_type_PE.sam" > "${Rbase}/spades/PE/remapping/test_ass_type_PE.bam"
samtools view -b "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE.sam" > "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE.bam"
samtools view -b "${Rbase}/metaspades/PE/remapping/test_ass_type_PE.sam" > "${Rbase}/metaspades/PE/remapping/test_ass_type_PE.bam"
samtools view -b "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE.sam" > "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE.bam"

# Sort the bam files
echo "Sort bam..."
samtools sort "${Rbase}/spades/allSE/remapping/test_ass_type_allSE.bam" > "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_sorted.bam"
samtools sort "${Rbase}/spades/PE/remapping/test_ass_type_PE.bam" > "${Rbase}/spades/PE/remapping/test_ass_type_PE_sorted.bam"
samtools sort "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE.bam" > "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam"
samtools sort "${Rbase}/metaspades/PE/remapping/test_ass_type_PE.bam" > "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_sorted.bam"
samtools sort "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE.bam" > "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam"

# Index the sorted bams
echo "Index sorted bam..."
samtools index "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_sorted.bam"
samtools index "${Rbase}/spades/PE/remapping/test_ass_type_PE_sorted.bam"
samtools index "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam"
samtools index "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_sorted.bam"
samtools index "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam"

# Get the statistics (contig lengths and average depth)
echo "Output statistics... "
samtools idxstats "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_sorted.bam" > "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_sorted.stat"
samtools idxstats "${Rbase}/spades/PE/remapping/test_ass_type_PE_sorted.bam" > "${Rbase}/spades/PE/remapping/test_ass_type_PE_sorted.stat"
samtools idxstats "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam" > "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_sorted.stat"
samtools idxstats "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_sorted.bam" > "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_sorted.stat"
samtools idxstats "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam" > "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_sorted.stat"

# Compute per base depth
echo "Compute cov depth..."
samtools depth -a "${Rbase}/spades/allSE/remapping/test_ass_type_allSE_sorted.bam"  > "${Rbase}/spades/allSE/remapping/test_ass_type_allSE.depth.tsv"
samtools depth -a "${Rbase}/spades/PE/remapping/test_ass_type_PE_sorted.bam" > "${Rbase}/spades/PE/remapping/test_ass_type_PE.depth.tsv"
samtools depth -a "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam" > "${Rbase}/spades/PE_SE/remapping/test_ass_type_PE_SE.depth.tsv"
samtools depth -a "${Rbase}/metaspades/PE/remapping/test_ass_type_PE_sorted.bam" > "${Rbase}/metaspades/PE/remapping/test_ass_type_PE.depth.tsv"
samtools depth -a "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE_sorted.bam" > "${Rbase}/metaspades/PE_SE/remapping/test_ass_type_PE_SE.depth.tsv"


# hisat2 Parameters:
#   -p  threads to use
#   -k  Primary alignments to report; the best alignment will be chosen from them
#   -x  specify the path and base name to an reference seq index

# samtools Parameters:
#   -b  input is BAM
#   -a  output all positions in contig, even if depth is 0