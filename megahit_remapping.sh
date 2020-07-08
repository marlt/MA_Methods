#!/bin/bash

# Script for remapping viral and NA reads with higher GC content than 25% on megahit contigs, using bwa mem

INDIR="/data/fass2/projects/ma_neander_assembly/MA_data/assembly/megahit"
INDIR_B="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped"

# VIRAL

IN1="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_viral_pair1_trimmed_filterGC.fastq"
IN2="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_viral_pair2_trimmed_filterGC.fastq"
IN_M="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_viral_allSE_trimmed_filterGC.fastq"
REF="${INDIR}/viral/bigger25/final_contigs_sorted.fa"

# BWA

mkdir -p "${INDIR}/viral/bigger25/remapping/bwa/"

echo 'Remapping with BWA-mem...'
echo 'Index...'
nice bwa index -p "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25" "${REF}"
echo 'Remap...'
nice bwa mem -t 30 -k 11 "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25" "${IN1}" "${IN2}" > "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25_PE.sam" 2> "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25_PE.log"
nice bwa mem -t 30 -k 11 "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25" "${IN_M}" > "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25_SE.sam" 2> "${INDIR}/viral/bigger25/remapping/bwa/viral_bigger25_SE.log"


# NA

IN1="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_NA_pair1_trimmed_filterGC.fastq"
IN2="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_NA_pair2_trimmed_filterGC.fastq"
IN_M="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_NA_allSE_trimmed_filterGC.fastq"
REF="${INDIR}/NA/bigger25/final_contigs_sorted.fa"

## BWA

mkdir -p "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25"

echo 'Remapping with BWA-mem...'
echo 'Index...'
nice bwa index -p "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25" "${REF}"
echo 'Remap...'
nice bwa mem -t 30 -k 11 "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25" "${IN1}" "${IN2}" > "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25_PE.sam" 2> "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25_PE.log"
nice bwa mem -t 30 -k 11 "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25" "${IN_M}" > "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25_SE.sam" 2> "${INDIR}/NA/bigger25/remapping/bwa/NA_bigger25_SE.log"


# Remapping statistics
# samtools flagstats to print the remapping statistics (mapped/unmapped read counts)
for DIR in "viral" "NA"
do
    echo ${DIR}
    DN="${INDIR}/${DIR}"
    if [ "${DIR}" == 'NA' ]
    then
        BN="NA_bigger25"
    elif [ "${DIR}" == 'viral' ]
    then
        BN="viral_bigger25"
    fi
    echo ${BN}
    echo ${DN}
    # convert into bam
    samtools view -b "${DN}/bigger25/remapping/bwa/${BN}_PE.sam" -o "${DN}/bigger25/remapping/bwa/${BN}_PE.bam"
    samtools view -b "${DN}/bigger25/remapping/bwa/${BN}_SE.sam" -o "${DN}/bigger25/remapping/bwa/${BN}_SE.bam"

    # sort according to position in genome -> left most coordinates first
    samtools sort "${DN}/bigger25/remapping/bwa/${BN}_PE.bam" -o "${DN}/bigger25/remapping/bwa/${BN}_PE_Psorted.bam"
    samtools sort "${DN}/bigger25/remapping/bwa/${BN}_SE.bam" -o "${DN}/bigger25/remapping/bwa/${BN}_SE_Psorted.bam"

    # flagstat reports comprising statistics on mapped/unmapped reads on target
    samtools flagstat "${DN}/bigger25/remapping/bwa/${BN}_PE_Psorted.bam" > "${DN}/bigger25/remapping/bwa/${BN}_PE_Psorted.flagstat"
    samtools flagstat "${DN}/bigger25/remapping/bwa/${BN}_SE_Psorted.bam" > "${DN}/bigger25/remapping/bwa/${BN}_SE_Psorted.flagstat"

done

# HiSat2 Parameters:
#   -p  number of threads to use
#   -k  Primary alignments to report; the best alignment will be chosen from them
#   -x  specify the path and base name to an reference seq index

# Segemehl Parameters.
#   -t  number of threads to use
#   -q  specify pair1 input file. if -p is not specified, this is treated as single-end input file 
#   -p  specify pair2 input file
#   -x  path to write index (segemehl indexing)
#   -i  path to find index (segemehl remapping)
#   -d  path to reference 

# Bwa mem/index Parameters:
#   -p  path to write index
#   -t  number of threads to use
#   -k  k-mer size in INT for seeding

# samtools Parameters:
#   -b  input is BAM
#   -o  output file
