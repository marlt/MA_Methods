#!/bin/bash
# BLAST contigs against a DB

QUERY="$1"              # multi fastq file containing query sequences
DB="$2"                 # path and basename of the target database
OUT="$3"                # csv file to write results to


TMP=$(dirname ${OUT})   # Directory, where to save the blast output table
# Write csv headline
echo "qseqid	qlen	qstart	qend	qseqid	sacc	staxid	ssciname	slen	length	sstart	send	score	pident	evalue	qcovs" > ${OUT}

# run actual blastn command - comparing nucleotides to nucleotides
nice blastn -num_threads 40 -task megablast -query ${QUERY} -db ${DB} -out ${TMP}/tmp.csv -evalue 0.001 -max_target_seqs 1 -outfmt "6 qseqid qlen qstart qend qseq sacc staxid ssciname slen length sstart send score pident evalue qcovs" 

cat ${TMP}/tmp.csv >> ${OUT}
rm ${TMP}/tmp.csv

# blastn    parameters:
#   -outfmt         6 means tabular: 	
#	    qseqid	    Query Seq-id, 
#	    qlen		Query sequence length,
#	    ssciname    Subject Scientific Name,
#	    staxid	    Subject Taxonomy IDS
#	    sacc	    Subject Accession
#	    sseqid	    Subject Seq-id
#	    slen	    Subject length	
#	    score	    Raw score
#		pident	    Percentage of identical matches
#		evalue	    Expect value
#		qstart	    Start of alignment in query
#		qend	    End of alignment in query
#		qseqid	    Aligned part of query sequence
#		length	    Alignment length
#		sstart	    Start of alignment in subject
#		send	    End of alignment in subject
#		qcovs	    query coverage on subject
#   -evalue         threshold for the value, how likely a sequence of the query seqs length is found on the particular db-size by chance
#   -db             specify the path to the target db basename
#   -max_target_seqs    number of target seqs with similarity to report 
#   -task           blast task to execute, algorithms with different sensitivities and alignment settings
#   -num_threads    number of available threads to use
#   -query          fastq file containing the query sequences

