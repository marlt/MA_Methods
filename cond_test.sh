#!/bin/bash

# assembly condition testings
# generate several assemblies with SPAdes; compared are:
#   - k-mer sizes: auto-detection, set to 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127
#   - read correction: enabled, disabled
# tested on the pooled viral read bin, GC > 25%

# input/output variables:
INPATH='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming/extracted/trimming/GC_filter/'
OUT='/data/fass2/projects/ma_neander_assembly/MA_data/assembly/spades'

mkdir -p ${OUT}/viral/bigger25

########### Test conditions ##############
# create an own directory for each assembly
mkdir -p ${OUT}/cond_test/auto_kmer/corrected
mkdir -p ${OUT}/cond_test/auto_kmer/uncorrected
mkdir -p ${OUT}/cond_test/set_kmer/corrected
mkdir -p ${OUT}/cond_test/set_kmer/uncorrected

# define variables for each 
ac="${OUT}/cond_test/auto_kmer/corrected"
au="${OUT}/cond_test/auto_kmer/uncorrected"
sc="${OUT}/cond_test/set_kmer/corrected"
su="${OUT}/cond_test/set_kmer/uncorrected"

# run the commands:
echo 'Assemble pooled viral uncorrected with set kmers'
nice spades.py -t 40 --only-assembler -s "${INPATH}/bigger_25/all_viral_pooled_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${su}" > "${su}/spades_pooled_trimmed.log"
echo 'Assemble pooled viral uncorrected with auto kmers'
nice spades.py -t 40 --only-assembler -s "${INPATH}/bigger_25/all_viral_pooled_trimmed_filterGC.fastq" -o "${au}" > "${au}/spades_pooled_trimmed.log"
echo 'Assemble pooled viral corrected with set kmers'
nice spades.py -t 40 -s "${INPATH}/bigger_25/all_viral_pooled_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127  -o "${sc}" > "${sc}/spades_pooled_trimmed.log"
echo 'Assemble pooled viral corrected with auto kmers'
nice spades.py -t 40 -s "${INPATH}/bigger_25/all_viral_pooled_trimmed_filterGC.fastq" -o "${ac}" > "${ac}/spades_pooled_trimmed.log"

# SPAdes parameter:
#   -t  NUM Threads to use
#   -s  single-end input reads
#   --only-assembler    run only the assembler module, read correction is disabled
#   -k  specify a comma-separated list of k-mers to use for assembly graph construction; if not defined, default is 'auto'
#   -o  output folder to write to
