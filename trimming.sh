#!/bin/bash

# trims and filters all ssDNA sequencing library reads and creates fastqc reports


DIR="$1"									# specify the directory all input fastq files are located
OUT_F="$2"									# specify output directory for trimmed fastq files
QC_DIR="$3"									# specify output directory for fastqc files

cd ${DIR} # move into the directory that contains all fastqs
for file in ../S*__um_sorted.bam			# iterate over parent directory of the fastq files - files with the specified names provide the correct basename for the fastq files
do
	BN=$(basename ${file} __um_sorted.bam)
	
	echo "${BN} processing"
	fname1="${BN}__pair1.fastq"						# specify pair1 reads
	fname2="${BN}__pair2.fastq"						# specify pair2 reads
	fnameS="${BN}__allSE.fastq"						# specify unflagged + singleton reads (merged into one, priorly)
	fname1_out="${BN}__pair1_trimmed.fastq"			# specify trimmed pair1 reads
	fname2_out="${BN}__pair2_trimmed.fastq"			# specify trimmed pair2 reads
	fname1_out_UP="${BN}__pair1UP_trimmed.fastq"	# specify trimmed pair1 reads without correpsonding pair2 due to trimming
	fname2_out_UP="${BN}__pair2UP_trimmed.fastq"	# specify trimmed pair2 reads without correpsonding pair1 due to trimming
	fnameS_out="${BN}__allSE_trimmed.fastq"			# specify trimmed unflagged + singleton reads
	f_outP="${BN}_PE_failed.fastq"					# write pair reads that do not pass trimming
	f_outS="${BN}_SE_failed.fastq"					# write singleton/unflag reads that do not pass trimming

	### paired-end trimming
	echo "paired-end trimming..."					
	# run fastp for paired-end reads
	fastp -M 24 -l 20 -3 -5 -w 8 -W 4 -n 5 -q 24 -u 25 -y -x -c -i ./${fname1} -I ./${fname2} -o ./trimming/${fname1_out} -O ./trimming/${fname2_out}  --unpaired1 ./trimming/${fname1_out_UP} --unpaired2 ./trimming/${fname2_out_UP} --adapter_fasta ./adapter_seqs.fasta --failed_out ./qual_fail/${f_outP} -h ${QC_DIR}/trimmed/fastp/${BN}__PE_trimmed.html &> ./trimming/${BN}_PE_trimmed.log
	
	# generate fastqc reports for all output fastq files of paired-end trimming
	fastqc -t 8 ./trimming/${fname1_out} -o ${QC_DIR}/trimmed &> "${QC_DIR}/trimmed/${BN}__pair1_trimmed_fastqc.log"
	fastqc -t 8 ./trimming/${fname2_out} -o ${QC_DIR}/trimmed &> "${QC_DIR}/trimmed/${BN}__pair2_trimmed_fastqc.log"
	fastqc -t 8 ./trimming/${fname1_out_UP}  -o ${QC_DIR}/trimmed &> "${QC_DIR}/trimmed/${BN}__pair1UP_trimmed_fastqc.log"
	fastqc -t 8 ./trimming/${fname2_out_UP}  -o ${QC_DIR}/trimmed &> "${QC_DIR}/trimmed/${BN}__pair2UP_trimmed_fastqc.log"
	fastqc -t 8 ./qual_fail/${f_outP} -o ${QC_DIR}/qual_fail &> "${QC_DIR}/qual_fail/${BN}_PE_failed.log"
	
	### single-end trimming
	echo "single-end trimming..."					
	# run fastp for unflagged/singleton reads
	fastp -M 24 -l 20 -3 -5 -w 8 -W 4 -n 5 -q 24 -u 25 -y -x -c -i ./${fnameS} -o ./trimming/${fnameS_out}  --adapter_fasta ./adapter_seqs.fasta --failed_out ./qual_fail/${f_outS} -h ${QC_DIR}/fastp/${BN}__SE_trimmed.html &> ./trimming/${BN}_SE_trimmed.log
	
	# generate fastqc reports for all output fastq files of single-end trimming
	fastqc -t 8 ./trimming/${fnameS_out} -o ${QC_DIR}/trimmed &> "${QC_DIR}/trimmed/${BN}__allSE_trimmed_fastqc.log"
	fastqc -t 8 ./qual_fail/${f_outS} -o ${QC_DIR}/qual_fail &> "${QC_DIR}/qual_fail/${BN}_SE_failed.log"	
done


# fastp Parameters:
#	-M	mean quality requirement shared by cut-front (-5), cut-tail (-3) or cut_sliding 
#	-l	minimum read length to retain
#	-3	quality trimming from 3' end
#	-5	quality trimming from 5' end
#	-w	number of threads to use
#	-W	window size to use for cut-front, cut-tail or cut-sliding
#	-n	N base threshold, if INT is exceeded in READ, READ is discarded
#	-q	quality value that a base is qualified (equal or higher)
#	-u	percents of bases allowed to be unqualified
#	-y	low complexity filtering; complexity is defined as percentage of bases that are different from its next base (base[i] != base[i+1])
#	-x	trim polyX tails of minimum length 10
#	-c	base correction for overlapping regions of paired-end reads
#	-i	input file, in paired-end modus this is pair1
#	-I	only paired-end: input file pair2
#	-o	output file, in paired-end modus this is pair1
#	-O	only paired-end: output file pair2
#	--unpaired1	only paired-end: output pair1 reads that lack pair2 due to trimming
#	--unpaired2	only paired-end: output pair2 reads that lack pair1 due to trimming
#	--adapter_fasta	provide a FASTA containing adapter sequences for manual adapter clipping
#	--failed-out	write all reads that failed trimming
#	-h	specify name for html report file

# fastqc parameters:
#	-t	NUM of threads to use
#	-o	Create OUTDIR and write all output files here