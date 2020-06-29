#!/bin/bash

# Load fastq entries to python dictionary and look up the headers from txt one by one in the dictionary.
# Search in dictionary should be in O(1) therefore it should be fast
# This script shall iterate over all fastq files and load the according viral, NA, and bacterial headers:

# specify scripts directory
SCRIPTS="/data/mahlzeitlocal/projects/ma_neander_assembly/anc_virus_MA/scripts"
# specify the directory containing the mapping fo read header ID to phylogenetic path
TXT="/data/fass2/projects/ma_neander_assembly/MA_data/kmer_class/clark/19mer/processed/"
# fastq files to extract from
FASTQ='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming'

regex1='.*allSE.*'
regex2='.*pair1.*'
regex3='.*pair2.*'

# iterate over all quality enriched fastq files (allSE, pair1, pair2 per input bam) and get reads that were assigned to a viral, bacterial or no taxon on the first assignment 
for file in ${FASTQ}/*_trimmed.fastq
do
    echo ${file}
    if [[ ${file} =~ ${regex1} ]]           # allSE
    then
        BN=$(basename ${file} _trimmed.fastq)
        # get the allSE bacterial and NA reads in one file per CLARK output
        ${SCRIPTS}/extract_reads_dicts.py ${file} "${TXT}/${BN}_viral_hit1.txt" "${TXT}/${BN}_NA_hit1.txt" "${TXT}/${BN}_bacterial_hit1.txt" "${FASTQ}/extracted/"
    
    elif [[ ${file} =~ ${regex2} ]]         # pair1
    then
        BN=$(basename ${file} _pair1_trimmed.fastq)
        # get the pair1 bacterial and NA reads 
        ${SCRIPTS}/extract_reads_dicts.py ${file} "${TXT}/${BN}_PE_viral_hit1.txt" "${TXT}/${BN}_PE_NA_hit1.txt" "${TXT}/${BN}_PE_bacterial_hit1.txt" "${FASTQ}/extracted/"
    
    elif [[ ${file} =~ ${regex3} ]]         # pair2
    then
        BN=$(basename ${file} _pair2_trimmed.fastq)
        # get the pair2 bacterial and NA reads 
        ${SCRIPTS}/extract_reads_dicts.py ${file} "${TXT}/${BN}_PE_viral_hit1.txt" "${TXT}/${BN}_PE_NA_hit1.txt" "${TXT}/${BN}_PE_bacterial_hit1.txt" "${FASTQ}/extracted/"
    fi
done 
# Merge all allSE, pair1 and pair2 separately
# Bacteria: allSE, pair1, pair2:
echo "Merge bacterial reads..."
for file in ${FASTQ}/extracted/*bacteria*.fastq
do
    if [[ ${file} =~ ${regex1} ]]
    then
    
        cat ${file} >> all_bacterial_allSE.fastq
    
    elif [[ ${file} =~ ${regex2} ]]
    then
    
        cat ${file} >> all_bacteria_pair1.fastq
    
    elif [[ ${file} =~ ${regex3} ]]
    then
        
        cat ${file} >> all_bacteria_pair2.fastq
    fi
done

# viral: allSE, pair2, pair2
echo "Merge viral reads..."
for file in ${FASTQ}/extracted/*viral*.fastq
do
    if [[ ${file} =~ ${regex1} ]]
    then
    
        cat ${file} >> all_viral_allSE.fastq
    
    elif [[ ${file} =~ ${regex2} ]]
    then
    
        cat ${file} >> all_viral_pair1.fastq
    
    elif [[ ${file} =~ ${regex3} ]]
    then
        
        cat ${file} >> all_viral_pair2.fastq
    fi
done

# NA: allSE, pair1, pair2:
echo 'Merge NA reads...'
for file in ${FASTQ}/extracted/*NA*.fastq
do
    if [[ ${file} =~ ${regex1} ]]
    then
    
        cat ${file} >> all_NA_allSE.fastq
    
    elif [[ ${file} =~ ${regex2} ]]
    then
    
        cat ${file} >> all_NA_pair1.fastq
    
    elif [[ ${file} =~ ${regex3} ]]
    then
        
        cat ${file} >> all_NA_pair2.fastq
    fi
done
