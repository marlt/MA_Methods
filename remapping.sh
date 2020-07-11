#!/bin/bash

# 07.04.20
# Script for comparable remapping using hisat2, segemehl and bwa mem of corrected reads on assembled viral, NA and bacterial contigs

# General input paths
INDIR="/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades"
INDIR_B="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped"


####### VIRAL #######################

# corrected input reads
IN1="${INDIR}/viral/bigger25/corrected/all_viral_pair1_trimmed_filterGC.00.0_0.cor.fastq"
IN2="${INDIR}/viral/bigger25/corrected/all_viral_pair2_trimmed_filterGC.00.0_0.cor.fastq"
IN_M="${INDIR}/viral/bigger25/corrected/all_viral_allSE_trimmed_filterGC.00.0_1.cor.fastq"
REF="${INDIR}/viral/bigger25/contigs_len300_cov10.fasta"

# create own directory for each mapper
mkdir -p ${INDIR}/viral/bigger25/remapping/hisat2
mkdir -p ${INDIR}/viral/bigger25/remapping/segemehl
mkdir -p ${INDIR}/viral/bigger25/remapping/bwa/len300_cov10

##### the HiSat2 mappings for viral reads > 25% GC had been created during assembly testings
# HiSat2
# echo "Processing VIRAL"
# echo "Remapping with HISAT2..."
# echo "Indexing..."
# nice hisat2-build -p 40 ${REF} ${INDIR}/viral/bigger25/remapping/hisat2/viral_bigger25 &> ${INDIR}/viral/bigger25/remapping/hisat2/viral_bigger25_indexing.log
# echo "Remapping..."
# hisat2 --version
# which hisat2
# nice hisat2 -p 40 -k 15 -x "${INDIR}/viral/bigger25/remapping/hisat2/viral_bigger25" -1 "${IN1}"  -2 "${IN2}" -U "${IN_M}"  -S "${INDIR}/viral/bigger25/remapping/hisat2/viral_bigger25.sam" --un "${INDIR}/viral/bigger25/remapping/hisat2/viral_bigger25_unmapped.fastq" &> "${INDIR}/viral/bigger25/remapping/hisat2/viral_bigger25_remapping.log"

# Segemehl
echo "Remapping with SEGEMEHL.."
# which segemehl.x 
echo " Indexing..."
nice segemehl.x -t 20 -x "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25.idx" -d ${REF} &> "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25_indexing.log"
echo "Remapping..."
# paired-end data
nice segemehl.x -t 30 -i "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25.idx" -d ${REF} -q "${IN1}" -p "${IN2}" > "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25_PE.sam" 2> "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25_PE.log"
# merged/single-end reads
nice segemehl.x -t 30 -i "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25.idx" -d ${REF} -q "${IN_M}"  > "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25_SE.sam" 2> "${INDIR}/viral/bigger25/remapping/segemehl/viral_bigger25_SE.log"

# BWA
echo 'Remapping with BWA-mem...'
echo 'Index...'
nice bwa index -p "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25" "${REF}"
echo 'Remap...'
nice bwa mem -t 30 -k 11 "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25" "${IN1}" "${IN2}" > "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25_PE.sam" 2> "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25_PE.log"
nice bwa mem -t 30 -k 11 "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25" "${IN_M}" > "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25_SE.sam" 2> "${INDIR}/viral/bigger25/remapping/bwa/len300_cov10/viral_bigger25_SE.log"


####### NA ########################
# corrected input reads
IN1="${INDIR}/NA/bigger25/corrected/all_NA_pair1_trimmed_filterGC.00.0_0.cor.fastq"
IN2="${INDIR}/NA/bigger25/corrected/all_NA_pair2_trimmed_filterGC.00.0_0.cor.fastq"
IN_M="${INDIR}/NA/bigger25/corrected/all_NA_allSE_trimmed_filterGC.00.0_1.cor.fastq"
REF="${INDIR}/NA/bigger25/contigs_len300_cov10.fasta"

# create own directory for each mapper
mkdir -p ${INDIR}/NA/bigger25/remapping/hisat2
mkdir -p ${INDIR}/NA/bigger25/remapping/segemehl
mkdir -p ${INDIR}/NA/bigger25/remapping/bwa/


## HiSat2
echo "Processing NA"
echo "Remapping with HiSat2..."
echo "Indexing..."
nice hisat2-build -p 40 "${REF}" "${INDIR}/NA/bigger25/remapping/hisat2/NA_bigger25" &> "${INDIR}/NA/bigger25/remapping/hisat2/NA_bigger25_indexing.log"
echo "Remapping..."
nice hisat2 -p 40 -k 15 -x "${INDIR}/NA/bigger25/remapping/hisat2/NA_bigger25" -1 "${IN1}"  -2 "${IN2}" -U "${IN_M}"  -S "${INDIR}/NA/bigger25/remapping/hisat2/NA_bigger25.sam" --un "${INDIR}/NA/bigger25/remapping/hisat2/NA_bigger25_unmapped.fastq" &> "${INDIR}/NA/bigger25/remapping/hisat2/NA_bigger25_remapping.log"

# Segemehl
echo "Remapping with SEGEMEHL.."
# which segemehl.x 
echo " Indexing..."
nice segemehl.x -t 20 -x "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25.idx" -d ${REF} &> "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25_indexing.log"
echo "Remapping..."
nice segemehl.x -t 30 -i "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25.idx" -d ${REF} -q "${IN1}" -p "${IN2}" > "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25_PE.sam" 2> "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25_PE.log"
nice segemehl.x -t 30 -i "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25.idx" -d ${REF} -q "${IN_M}"  > "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25_SE.sam" 2> "${INDIR}/NA/bigger25/remapping/segemehl/NA_bigger25_SE.log"

# # BWA
echo 'Remapping with BWA-mem...'
echo 'Index...'
nice bwa index -p "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25" "${REF}"
echo 'Remap...'
nice bwa mem -t 30 -k 11 "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25" "${IN1}" "${IN2}" > "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25_PE.sam" 2> "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25_PE.log"
nice bwa mem -t 30 -k 11 "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25" "${IN_M}" > "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25_SE.sam" 2> "${INDIR}/NA/bigger25/remapping/bwa/len300_cov10/NA_bigger25_SE.log"


####### BACTERIA ####################
# input are the uncorrected reads - bayeshammer is too memory exhaustive to be run on this reads

IN1="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_bacteria_pair1_trimmed_filterGC.fastq"
IN2="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_bacteria_pair2_trimmed_filterGC.fastq"
IN_M="${INDIR_B}/fastq_MA/trimming/extracted/trimming/GC_filter/bigger_25/all_bacteria_allSE_trimmed_filterGC.fastq"
REF="${INDIR}/bacteria/bigger25/contigs_len300_cov10.fasta"

mkdir -p ${INDIR}/bacteria/bigger25/remapping/hisat2
mkdir -p ${INDIR}/bacteria/bigger25/remapping/segemehl
mkdir -p ${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/

# HiSat2
echo "Remapping bacteria reads..."
echo "Indexing..."
nice hisat2-build -p 40 "${REF}" "${INDIR}/bacteria/bigger25/remapping/hisat2/bacteria_bigger25" &> "${INDIR}/bacteria/bigger25/remapping/hisat2/bacteria_bigger25_indexing.log"
echo "Remapping..."
nice hisat2 -p 40 -k 15 -x "${INDIR}/bacteria/bigger25/remapping/hisat2/bacteria_bigger25" -1 "${IN1}"  -2 "${IN2}" -U "${IN_M}"  -S "${INDIR}/bacteria/bigger25/remapping/hisat2/bacteria_bigger25.sam" --un "${INDIR}/bacteria/bigger25/remapping/hisat2/bacteria_bigger25_unmapped.fastq" &> "${INDIR}/bacteria/bigger25/remapping/hisat2/bacteria_bigger25_remapping.log"

# Segemehl
echo "Remapping with SEGEMEHL.."
which segemehl.x 
echo " Indexing..."
nice segemehl.x -t 20 -x "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25.idx" -d ${REF} &> "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25_indexing.log"
echo "Remapping..."
nice segemehl.x -t 30 -i "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25.idx" -d ${REF} -q "${IN1}" -p "${IN2}" > "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25_PE.sam" 2> "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25_PE.log"
nice segemehl.x -t 30 -i "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25.idx" -d ${REF} -q "${IN_M}" > "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25_SE.sam" 2> "${INDIR}/bacteria/bigger25/remapping/segemehl/bacteria_bigger25_SE.log"

#  BWA
file ${REF}
file ${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25*
file ${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/
file ${IN_M}
echo 'Remapping with BWA-mem...'
echo 'Index...'
nice bwa index -p "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25" "${REF}"
echo 'Remap...'
nice bwa mem -t 30 -k 11 "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25" "${IN1}" "${IN2}" > "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25_PE.sam" 2> "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25_PE.log"
nice bwa mem -t 30 -k 11 "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25" "${IN_M}" > "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25_SE.sam" 2> "${INDIR}/bacteria/bigger25/remapping/bwa/len300_cov10/bacteria_bigger25_SE.log"



############################ FLAGSTATS #################################

# Remapping statistics computed for each assembly and remapper 
for DIR in "bacteria" "viral" "NA"
do
    DN=${INDIR}/${DIR}
    echo "Pfad:     ${INDIR}/${DIR}"
    if [ "${DIR}" == 'bacteria' ]
    then
        BN="bacteria_bigger25"
    elif [ "${DIR}" == 'viral' ]
    then
        BN="viral_bigger25"
    elif [ "${DIR}" == 'NA' ]
    then
        BN="NA_bigger25"
    fi
    # convert into bam
    samtools view -b "${DIR}/bigger25/remapping/hisat2/${BN}.sam" -o "${DIR}/bigger25/remapping/hisat2/${BN}.bam"
    samtools view -b "${DIR}/bigger25/remapping//segemehl/${BN}_PE.sam" -o "${DIR}/bigger25/remapping/segemehl/${BN}_PE.bam"
    samtools view -b "${DIR}/bigger25/remapping/segemehl/${BN}_SE.sam" -o "${DIR}/bigger25/remapping/segemehl/${BN}_PE.bam"
    samtools view -b "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_PE.sam" -o "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_PE.bam"
    samtools view -b "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_SE.sam" -o "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_SE.bam"

    # sort according to position in genome -> leftmost coordinates first
    samtools sort "${DIR}/bigger25/remapping/hisat2/${BN}.bam" -o "${DIR}/bigger25/remapping/hisat2/${BN}_Psorted.bam"
    samtools sort "${DIR}/bigger25/remapping/segemehl/${BN}_PE.bam" -o "${DIR}/bigger25/remapping/segemehl/${BN}_PE_Psorted.bam"
    samtools sort "${DIR}/bigger25/remapping/segemehl/${BN}_SE.bam" -o "${DIR}/bigger25/remapping/segemehl/${BN}_SE_Psorted.bam"
    samtools sort "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_PE.bam" -o "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_PE_Psorted.bam"
    samtools sort "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_SE.bam" -o "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_SE_Psorted.bam"

    # flagstat reports comprising statistics on mapped/unmapped reads on target
    samtools flagstat "${DIR}/hisat2/${BN}_Psorted.bam" > "${DIR}/hisat2/${BN}_Psorted.flagstat"
    samtools flagstat "${DIR}/segemehl/${BN}_PE_Psorted.bam" > "${DIR}/segemehl/${BN}_PE_Psorted.flagstat"
    samtools flagstat "${DIR}/segemehl/${BN}_SE_Psorted.bam" > "${DIR}/segemehl/${BN}_SE_Psorted.flagstat"
    samtools flagstat "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_PE_Psorted.bam" > "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_PE_Psorted.flagstat"
    samtools flagstat "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_SE_Psorted.bam" > "${DN}/bigger25/remapping/bwa/len300_cov10/${BN}_SE_Psorted.flagstat"
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
#   -a  output all positions in contig, even if depth is 0