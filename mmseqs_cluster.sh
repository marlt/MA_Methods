#!/bin/bash

# Script for clustering nucleotide sequences and extracting their representative sequences

INFILE="$1"     # Input Multi-Fasta File
QUERYDB="$2"    # Dir to write query db to
PRE="$3"        # Set the prefix for the querydb
OUT="$4"        # Dir to write all output files to
TMP="$5"        # dir for temporary files of clustering step
CMODE="$6"      # Clustering mode, see parameter explanation
COV="$7"        # Minimal Sequency Coverage of Alignment for Clustering Sequences
EVAL="$8"       # Set E-Value threshold for alignment

mkdir -p ${QUERYDB}         # Create query output directory
mkdir -p ${OUT}             # Create the output path if not already present
mkdir -p ${TMP}             # Create directory for temporary files
mkdir -p ${OUT}/seqfiledb   # Create the output path if not already present
mkdir -p ${OUT}/repdb       # Create the output path if not already present

# Cluster Sequences and convert mmseqs format into tsv
nice mmseqs createdb ${INFILE} ${QUERYDB}/${PRE}
nice mmseqs cluster ${QUERYDB}/${PRE} ${OUT}/clust ${TMP} -e ${EVAL} -c ${COV} --cov-mode ${CMODE} --threads 30
nice mmseqs createtsv ${QUERYDB}/${PRE} ${QUERYDB}/${PRE} ${OUT}/clust ${OUT}/clust.tsv

# Output sequences in clusters to multi fasta
nice mmseqs createseqfiledb ${QUERYDB}/${PRE} ${OUT}/clust ${OUT}/seqfiledb/seqfiledb
nice mmseqs result2flat ${QUERYDB}/${PRE} ${QUERYDB}/${PRE} ${OUT}/seqfiledb/seqfiledb ${OUT}/seqfiledb/seqfiledb.fasta

# Extract representatives - The sequence with highest connectivity of the cluster
nice mmseqs createsubdb ${OUT}/clust ${QUERYDB}/${PRE} ${OUT}/repdb/repdb
nice mmseqs convert2fasta ${OUT}/repdb/repdb ${OUT}/repdb/repdb.fasta

# MMseqs Parameters:
#   -e  e-value threshold for all-against-all sequence search
#   -c  coverage/similarity threshold to qualify for an edge
#   --cov-mode  where the coverage threshold needs to be passed; here: the smaller of both sequences need to meet the coverage criterion
#   --threads   number of threads to use