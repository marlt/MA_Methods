#!/bin/bash


# testing input library configuration and assembly strategies of SPAdes and metaSPAdes
# input libraries: PE, PE + merged, and pooled single-end files
# parameters had been adapted from cond_test.sh

# input/output variables:
INPATH='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/'
OUT='/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades'

# Create a directory for each assembly
mkdir -p ${OUT}/test_ass_type/spades/allSE 
mkdir -p ${OUT}/test_ass_type/spades/PE
mkdir -p ${OUT}/test_ass_type/spades/PE_SE
mkdir -p ${OUT}/test_ass_type/pooled
mkdir -p ${OUT}/test_ass_type/metaspades/PE_SE

# Create all SPAdes assemblies
echo 'Assemble allSE...'
nice spades.py -t 40 -s "${INPATH}/bigger_25/all_viral_allSE_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/test_ass_type/spades/allSE" &> "${OUT}/test_ass_type/spades/allSE/spades_allSE_trimmed.log"
echo 'Assemble PE...'
nice spades.py -t 40  -2 "${INPATH}/bigger_25/all_viral_pair1_trimmed_filterGC.fastq" -1 "${INPATH}/bigger_25/all_viral_pair2_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127 -o "${OUT}/test_ass_type/spades/PE" &> "${OUT}/test_ass_type/spades/PE/spades_PE_trimmed.log"
echo 'Assemble PE and SE...'
nice spades.py -t 40  -2 "${INPATH}/bigger_25/all_viral_pair1_trimmed_filterGC.fastq" -1 "${INPATH}/bigger_25/all_viral_pair2_trimmed_filterGC.fastq" --merged "${INPATH}/bigger_25/all_viral_allSE_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/test_ass_type/spades/PE_SE" &> "${OUT}/test_ass_type/spades/PE_SE/spades_PE_SE_trimmed.log"
pooled assembly done during condition tests

# Create all metaSPAdes assemblies
echo 'Assemble meta PE...'
nice metaspades.py -t 40  -2 "${INPATH}/bigger_25/all_viral_pair1_trimmed_filterGC.fastq" -1 "${INPATH}/bigger_25/all_viral_pair2_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127 -o "${OUT}/test_ass_type/metaspades/PE" &> "${OUT}/test_ass_type/metaspades/PE/spades_PE_trimmed.log"
echo 'Assemble meta PE and SE...'
nice metaspades.py -t 40  -2 "${INPATH}/bigger_25/all_viral_pair1_trimmed_filterGC.fastq" -1 "${INPATH}/bigger_25/all_viral_pair2_trimmed_filterGC.fastq" --merged "${INPATH}/bigger_25/all_viral_allSE_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/test_ass_type/metaspades/PE_SE" &> "${OUT}/test_ass_type/metaspades/PE_SE/metaspades_PE_SE_trimmed.log"
# metaSPAdes only runs with minimum two paired end read files

# SPAdes/metaSPAdes parameters:
#   -t  number of threads to use
#   -s  input single-end reads
#   -1  input pair1 reads
#   -2  input pair2 reads
#   --merged    input merged reads 
#   -k  specify a list of k-mers to use for assembly
#   -o  output file for assembled contigs

