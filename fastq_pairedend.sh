#!/bin/bash

#10.01.20
#sort bam files, extract fastq and generate QC-reports with fastqc

DIR="$1"			# input
OUT="$2"			# output dir
OUT_QC="$3"			# outdir QC reports
cd ${DIR}

#sort and extract
for file in ${DIR}/*_unmapped.bam
do
	BN=$(basename ${file} unmapped.bam)
	samtools sort -n ${file} -O BAM -T tmp -o ${BN}_um_sorted.bam &> ${BN}_um_sorted.log					# sorting bam file is needed for proper multi fastq extraction 
	samtools flagstat ${file} > ${BN}_um_sorted.flagstat													# get number of paired and unpaired for validation of extracted read numbers
	samtools fastq -s ${OUT}/${BN}_single.fastq -0 ${OUT}/${BN}_unflag.fastq -1 ${OUT}/${BN}_pair1.fastq -2 ${OUT}/${BN}_pair2.fastq ${BN}_um_sorted.bam 		# extract pairs, singletons and unflagged
done


# fastqc reports
for file in ${OUT}/*.fastq
do
	BN=$(basename ${file} .fastq)
	fastqc -t 4 ${file} -o ${OUT_QC} > ${OUT_QC}/${BN}.log	 			# quality reports for all fastq files generated

done


# samtools sort parameters:
#	-n	sort by read name
# 	-O	output format (SAM, BAM, CRAM)
#	-T	Write temporary files with PREFIX.nnnn.bam
#	-o	Write final output to FILE instead of standard output

# samtools flagstat parameters:
# no parameters provided

# samtools fastq parameters:
#	-s	Write singletons to FILE
#	-O	write paired reads flagged both or neither READ1 and READ2 to FILE
#	-1	write paired reads flagged READ1 to FILE
#	-2	write paired reads flagged READ2 to FILE


# fastqc parameters:
#	-t	NUM of threads to use
#	-o	Create OUTDIR and write all output files here

