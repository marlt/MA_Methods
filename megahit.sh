#!/bin/bash

# metagenome assembly with megahit on NA and viral read sets

# input directory to fastq files
IN="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/"
# output dir for contigs
OUT="/data/fass2/projects/ma_neander_assembly/MA_data/assembly/megahit/"

# generate output dir if not already present
mkdir -p ${OUT}

# generate the contigs from viral reads; output all contigs
nice  megahit -t 30 --min-contig-len 0 -1 "${IN}/all_NA_pair1_trimmed_filterGC.fastq" -2 "${IN}/all_NA_pair2_trimmed_filterGC.fastq" --12 "${IN}/all_NA_allSE_trimmed_filterGC.fastq" --k-list 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127 -o ${OUT}/bigger25/NA
# generate the contigs from NA reads; output all contigs
nice  megahit -t 30 --min-contig-len 0 -1 "${IN}/all_viral_pair1_trimmed_filterGC.fastq" -2 "${IN}/all_viral_pair2_trimmed_filterGC.fastq" --12 "${IN}/all_viral_allSE_trimmed_filterGC.fastq" --k-list 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127 -o ${OUT}/bigger25/viral

# MEGAHIT Parameteres:
#   -t  number of thread to use
#   --min-contig-len    specify length INT for that the contig is qualified
#   -1  pair1 input reads
#   -2  pair2 input reads
#   -12 merged input reads
#   --k-list    k-mers to be used for assembly graph construction
#   -o  output directories