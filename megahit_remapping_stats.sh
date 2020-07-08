#!/bin/bash

# Remapping statistics for bwa applied to megahit contigs using viral and NA read sets

# Input paths to viral and NA directories
Virus="/data/fass2/projects/ma_neander_assembly/MA_data/assembly/megahit/viral/bigger25"
NA="/data/fass2/projects/ma_neander_assembly/MA_data/assembly/megahit/NA/bigger25"
# Type of the mapping tool
Remapper="bwa"

# Index sorted bam files
echo "Index sorted bam..."
echo "NA PE..."
samtools index "${NA}/remapping/${Remapper}/NA_bigger25_PE_Psorted.bam"
echo "NA SE..."
samtools index "${NA}/remapping/${Remapper}/NA_bigger25_SE_Psorted.bam"
echo "Viral PE..."
samtools index "${Virus}/remapping/${Remapper}/viral_bigger25_PE_Psorted.bam"
echo "Viral SE..."
samtools index "${Virus}/remapping/${Remapper}/viral_bigger25_SE_Psorted.bam"

# Get average depth and contig length information from idxstats
echo "Output statistics... "
echo "NA PE..."
samtools idxstats "${NA}/remapping/${Remapper}/NA_bigger25_SE_Psorted.bam" > "${NA}/remapping/${Remapper}/NA_bigger25_SE_Psorted.stat"
echo "NA SE..."
samtools idxstats "${NA}/remapping/${Remapper}/NA_bigger25_PE_Psorted.bam" > "${NA}/remapping/${Remapper}/NA_bigger25_PE_Psorted.stat"
echo "Viral PE..."
samtools idxstats "${Virus}/remapping/${Remapper}/viral_bigger25_SE_Psorted.bam" > "${Virus}/remapping/${Remapper}/viral_bigger25_SE_Psorted.stat"
echo "Viral SE..."
samtools idxstats "${Virus}/remapping/${Remapper}/viral_bigger25_PE_Psorted.bam" > "${Virus}/remapping/${Remapper}/viral_bigger25_PE_Psorted.stat"

# Comput per base depth
echo "Compute cov depth..."
echo "NA PE..."
samtools depth -a "${NA}/remapping/${Remapper}/NA_bigger25_SE_Psorted.bam" > "${NA}/remapping/${Remapper}/NA_bigger25_SE_Psorted.depth.tsv"
echo "NA SE..."
samtools depth -a "${NA}/remapping/${Remapper}/NA_bigger25_PE_Psorted.bam" > "${NA}/remapping/${Remapper}/NA_bigger25_PE_Psorted.depth.tsv"
echo "Viral PE..."
samtools depth -a "${Virus}/remapping/${Remapper}/viral_bigger25_SE_Psorted.bam" > "${Virus}/remapping/${Remapper}/viral_bigger25_SE_Psorted.depth.tsv"
echo "Viral SE..."
samtools depth -a "${Virus}/remapping/${Remapper}/viral_bigger25_PE_Psorted.bam" > "${Virus}/remapping/${Remapper}/viral_bigger25_PE_Psorted.depth.tsv"


# samtools Parameters:
#   -b  input is BAM
#   -o  output file
#   -a  output all positions in contig, even if depth is 0