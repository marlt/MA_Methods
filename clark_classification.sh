#!/bin/sh

# classification of multi-fastq/fasta sequence file on a nucleotide database using exact k-mer matching 
# this script wraps scripts from the CLARK classification software (v1.2.6.1)
# paths and parameters were modified to classify on three different k-mer sizes: 19, 25, 31

# create the main output directory
mkdir -p /data/fass2/projects/ma_neander_assembly/MA_data/kmer_class/clark/25mer/05conf

DIR='/data/fass2/projects/ma_neander_assembly/MA_data/kmer_class/clark/25mer'                                   # output directory used in the script
DATA='/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming'               # specify directory of input data

# download RefSeq database resources and configure for k-mer database creation
nice /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/set_targets.sh /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/db/ bacteria viruses fungi human # configure database to classify on
echo "set_targets finished"

cd /data/prostlocal/programs/clark_1.2.6.1/ # move to directory of the CLARK software to ensure proper functioning

# db creation and classification 

while read file;
do
	BN=$(basename ${file} _allSE_trimmed.fastq)         # get the basename of the files, to specify pair1, pair2 and allSE files per original input bam
	echo ${BN}
    # database creation and k-mer classification commands
    # firstly, the command tries to find a k-mer database in the directory specified with set_targets.sh; if not found, it builds it prior to classification
	nice /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/classify_metagenome.sh -O "${DATA}/${BN}_allSE_trimmed.fastq" -k 25 -m 0 -R "${DIR}/${BN}_allSE_trimmed.results" &>> ${DIR}/allSE_25mer_class.log  # single-ende classification
	wait "$!"           # wait until the previous process finished
    nice /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/classify_metagenome.sh -P "${DATA}/${BN}_pair1_trimmed.fastq" "${DATA}/${BN}_pair2_trimmed.fastq" -k 25 -m 0 -R "${DIR}/${BN}_PE_trimmed.results" &>> ${DIR}/allPE_25mer_class.log	    # paired-end classification
    wait "$!"           # wait until the previous process finished
done < ${DIR}/../samples.txt    # previously written: text-file containing fastq file names to classify on

    
cd /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/

# create krona output and generate html format visualization     
for FILE in ${DIR}/*.csv            # for all classification results (csv tables are output of classify_metagenome.sh)
do
    BN=$(basename ${FILE} .results.csv) 
    # estimate abundance of taxa; write krona output
    /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/estimate_abundance.sh -c 0.5 -F "${FILE}" -D /data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/db/ --krona 
    # echo finish
    mv results.krn "${DIR}/05conf/${BN}25mer05conf.krn"         # estimate_abundance always output a file named "results.krn" - rename output after generation to avoid information loss
    ktImportTaxonomy -o "${DIR}/05conf/${BN}25mer05conf.html" -m 3 "${DIR}/05conf/${BN}25mer05conf.krn"  # create html report from krona output
done 

# mkdir parameters:
#   -p  generate STR as whole path 

# set-targets.sh parameters:
#   "bacteria/human/viruses/fungi"  specify which RefSeq data to download prior to database creation

# classify_metagenome.sh parameters:
#   -O  input multi-fastq file to classify - single-end mode
#   -P  input multi-fastq files to classify - paired-end mode
#   -k  specify INT as k-mer size to classify on (default: 31)
#   -R  specify file for results (csv)

# estimate_abundance.sh parameters:
#   c   specify confidence threshold in range 0.5-1.0 (default 0.5)
#   F   specify input file (csv results file from classify_metagenome.sh)
#   D   target database to use, configured with set_targets
#   krona   write output file in krona compatible format

# ktImportTaxonomy parameters:
#   -o   output in html format
#   -m  column of input file to use as magnitude (= read counts per taxon)