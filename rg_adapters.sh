#!/bin/bash
# bash script to count exact matches of strings (in this case: adapter sequences) in a multi-fasta file
# input files are the adapter fasta and text file as well as the directory that contains all multi-fastq files to search in

DIR="$1"			# directory containing all fastq to be searched and adapter sequence files in fasta (with header) and txt format (only sequences, one per line)
BN_ADAP="$2"		# specify basename of the adapter containing files: .fasta and .txt
OUT="$3"			# output path and file


headline="Filename	Read_Counts"			# Headline to include in output

while read line								
do
	if [[ ${line} =~ ^\>.+ ]]				# if header line in fasta encountered
	then
		headline+="	${line#>}"				# generate a string with all adapter IDs, separated by tab
	fi
done < ${DIR}/${BN_ADAP}.fasta				# read every line in the adapter fasta file

echo ${headline} > ${OUT}					# write the header line to the defined output, containing all adapter IDs


for fastq in ${DIR}/*.fastq;					
do
	content="$(basename ${fastq} .fastq)"		# start content line with the fastq filename
	reads=$(sed -n '1~4 p' ${fastq} | wc -l)	# count read number in each fastq file
	content+="	${reads}"						# add read counts to string

	while read line
	do
		hits=$(rg -c ${line} ${fastq})   
		if [[ ${hits} =~ [0-9]+ ]]				
		then
			content+="	${hits}"			# extend  with counts of exact matches for each adapter
		else
			content+="	0"					# rg returns NoneType, if no matches occur; write 0 instead
		fi
	done < ${DIR}/${BN_ADAP}.txt			# 
	echo ${content} >> ${OUT}				# write line of read counts and adapter matches 
done


# ripgrep parameters:
#	-c	count hits instead of writing to STDOUT