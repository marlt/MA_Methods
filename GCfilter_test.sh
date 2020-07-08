#!/bin/bash
# Test Assembly with GC filtered reads of different degrees

# input/output variables:
INPATH='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/'
OUT='/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades'

########### Test GC-content ######################
# create an own directory for each assembly
mkdir -p ${OUT}/GCfilter_test/bigger25        
mkdir -p ${OUT}/GCfilter_test/bigger45
mkdir -p ${OUT}/GCfilter_test/25_45          
mkdir -p ${OUT}/GCfilter_test/all
mkdir -p ${OUT}/GCfilter_test/merged_25_45_+_bigger45

# run the assemblies
echo 'Assemble filtered 45...'
nice metaspades.py -t 40  -2 "${INPATH}/bigger_45/all_viral_pair1_trimmed_filterGC_45.fastq" -1 "${INPATH}/bigger_45/all_viral_pair2_trimmed_filterGC_45.fastq" --merged "${INPATH}/bigger_45/all_viral_allSE_trimmed_filterGC_45.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/GCfilter_test/bigger45/" &> "${OUT}/GCfilter_test/bigger45/metaspades_GC_45_trimmed.log"
echo 'Assemble filtered 25 and 45...'
# Use spades since metaspades fails on insert size estimation
nice spades.py -t 40  -1 "${INPATH}/25_45/all_viral_pair1_trimmed_filterGC_25_45.fastq" -2 "${INPATH}/25_45/all_viral_pair2_trimmed_filterGC_25_45.fastq" --merged "${INPATH}/25_45/all_viral_allSE_trimmed_filterGC_25_45.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/GCfilter_test/25_45/" &> "${OUT}/GCfilter_test/25_45/metaspades_GC_25_45_trimmed.log"
echo 'Assemble unfiltered...'
nice metaspades.py -t 40  -2 "${INPATH}/../all_viral_pair1_trimmed.fastq" -1 "${INPATH}/../all_viral_pair2_trimmed.fastq" --merged "${INPATH}/../all_viral_allSE_trimmed.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/GCfilter_test/all/" &> "${OUT}/GCfilter_test/all/metaspades_GC_all_trimmed.log"
echo 'Assemble filtered 25...'
nice metaspades.py  -t 40  -2 "${INPATH}/bigger_25/all_viral_pair1_trimmed_filterGC.fastq" -1 "${INPATH}/bigger_25/all_viral_pair2_trimmed_filterGC.fastq" --merged "${INPATH}/bigger_25/all_viral_allSE_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${OUT}/GCfilter_test/bigger25/" &> "${OUT}/GCfilter_test/bigger25/metaspades_GC_25_trimmed.log"

# conditions derived from results of cond_test.sh and test_lib.sh

# SPAdes/metaSPAdes parameters:
#   -t  number of threads to use
#   -s  input single-end reads
#   -1  input pair1 reads
#   -2  input pair2 reads
#   --merged    input merged reads 
#   -k  specify a list of k-mers to use for assembly
#   -o  output file for assembled contigs