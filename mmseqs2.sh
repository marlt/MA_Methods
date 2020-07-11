#!/bin/bash

# Script comprises all commands to perform a homology search with Mmseqs2

# define variables
QUERY="$1"                      # insert different types of input queries: VIRAL, NA, BACTERIA
 # concatinated fasta file of viral, bacterial, fungi and human sequences - create the mmseqs target db on this
DBFASTA="/data/fass1/database/blast_bak_virus_fungi_human/blastdb/bacteria_virus_fungi_human/all_bac_vir_fungi_homoS.fna" 
# target db directory             
TARGETDB="/data/fass1/database/blast_bak_virus_fungi_human/blastdb/bacteria_virus_fungi_human/mmseqs_db/all_bac_vir_fungi"


BN="$2"                     # which input data set is given: bacteria, viral, NA
DN="$(dirname ${QUERY})"    # dirname of the query seqs
TN="$(dirname ${TARGETDB})" # dirname of the target seqs db


QUERYDB="${DN}/${BN}_querydb"                           # basename of the query db
RESULTSDB="${DN}/mmseqs2_out/${BN}/${BN}_resultsdb"     # basename of the results db
OUT="${DN}/mmseqs2_out/final/${BN}_mmseqs2.tsv"         # output path for the results table
TMP="${TN}/TMP"                                         # directory for all temporary files

###### ALIGNMENT OPTIONS##############
ALNLEN=50           # minimal alignment length of 50, unlikely to have a hit this long by chance
EVAL=0.0001         # minimal evalue 0.0001 - rather restrictive
MSEQS=500           # prefilter step - max hits per query seq to retain
THREADS=30          # number of threads to use

# generate directories, if not present
mkdir -p "${TMP}"
mkdir -p "${DN}/mmseqs2_out/${BN}/"
mkdir -p "${DN}/mmseqs2_out/final/"


nice mmseqs createdb ${QUERY} ${QUERYDB}            # preprocess k-mer lists for the query sequences, create database
nice mmseqs createdb ${DBFASTA} ${TARGETDB}         # preprocess k-mer lists for the target sequences, create database
nice mmseqs createindex ${TARGETDB} ${TMP} --search-type 3 --check-compatible 1 --threads 30            # index the target db
# run homology search; looking for queries in targets and writing them to the results db
nice mmseqs search ${QUERYDB} ${TARGETDB} ${RESULTSDB} ${TMP} --threads ${THREADS} -e ${EVAL} --max-seqs ${MSEQS} --min-aln-len ${ALNLEN} --search-type 3
# convert alignment information into tab-separated blast-like output format contain common aln metrices
nice mmseqs convertalis ${QUERYDB} ${TARGETDB} ${RESULTSDB} ${OUT} --threads ${THREADS} --format-mode 0 --format-output 'query,target,taxid,taxname,bits,alnlen,pident,evalue,qlen,qstart,qend,qcov,tlen,tstart,tend,tcov' 


# Mmseqs Parameters:
#   --search-type   search type can be auto (0), amino acid (1), translated (2), nucleotide (3), translated nucleotide alignment (4)
#   --check-compatible  always recreate (0), check for present index (1), or fail if index is incompatible (2)
#   --threads   number of threads to use
#   -e          e-value threshold
#   --max-seqs     maximum target hits to pass prefiltering step
#   --min-aln-len   minimim alignment length
#   --format-mode   choose output format: BLAST tab (0), SAM (1), Blast tab + query/db length (2)
#   --format-output choose output columns by key word from 
