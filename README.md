# MA_Methods
This file shall supplement the Master's Thesis ``Identifying Bacteria and Viruses
from an Ancient Metagenomic Sample`` with detailed descriptions of all relevant commands and custom scripts. The structure is oriented on the Methods section which I recommend to conduct in case of uncertainties. 
<!-- Additional chapters are introduced if this is necessary to keep the pipeline character of the script explanations. -->
The parameter functionality of utilized tools is described within the custom scripts, if present. Custom scripts implementing command line interfaces (CLI) are explained in place. If not described further, default settings apply.

The programs had been used in a conda environment for better version control and dependency management. 

# Contents
[Data Preprocessing](#data-preprocessing)
[Extraction of Input Reads](#extraction-of-input-reads)
[Adapter Content Estimation](#adapter-content-estimation)
[Quality Enrichment](#quality-enrichment)
[K-mer based Read Classification](#k-mer-based-read-classification)
[Read Binning](#read-binning)
[Assembly Evaluation: K-mer Sizes and Read Correction](#assembly-evaluation-k-mer-sizes-and-read-correction)
[Assembly Evaluation: Input Library and Assembler Subtype](#assembly-evaluation-input-library-and-assembler-subtype)
[Assembly Evaluation: Contiguity and GC Filtering](#assembly-evaluation-contiguity-and-gc-filtering)
[Final Assembly Set Up](#final-assembly-set-up)
[Contig Homology Search on Comprehensive Databases](#contig-homology-search-on-comprehensive-databases)
[Contig Clustering](#contig-clustering)
[Clustering on Target Hit IDs](#clustering-on-target-hit-ids)
[Selection of Ancient Candidate Reference Sequences](#selection-of-ancient-candidate-reference-sequences)
[Creating a Microbacterium Consensus Sequence](#creating-a-microbacterium-consensus-sequence)
[Damage Pattern Analysis](#damage-pattern-analysis)
[Visualization with Matplotlib](#visualization-with-matplotlib)

## Data Preprocessing
The input data had been downloaded using a recursive ``wget`` command:
```
cd <Data_Directory>
wget -r -np -nd http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/unmapped_qualfail/

# wget Parameters:
#   -r  recursive download
#   -np exclude parent directories from download
#   -nd     o not create directories on local drive
```

This command downloads all files present in the specified directory, in this case 33 bam files. Since the read names were unconveniently long, files had been given shorter, unique names (commands not shown). A filename mapping is provided [here](https://osf.io/qbreu/files/).

### Extraction of Input Reads

The script `fastq_pairedend.sh` includes the following steps: sorting bam files by read names (``samtools sort``), calculating read mapping statistics (``samtools flagstat``) and extracting the read sequences (``samtools fastq``). The latter comprises reads flagged as pair1, pair2, singleton, or both or none of pair1 and pair2 in each bam file. This leads to 4 fastq files per input bam. In addition, it creates a fastqc report in html format assessing several quality measures for all fastq files.

```
BAM="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/"
FASTQ="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA"
QC="/data/mahlzeitlocal/projects/ma_neander_assembly/anc_virus_MA/qual_reports/raw"

<path-to-script>/fastq_pairedend.sh  "${BAM}" "${FASTQ}" "${QC}"
```

Fastq files labelled 'unflag' correspond to merged read pairs generated for the [Neanderthal sequencing project](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4031459/) which is explained in the introductory section of the thesis. For convenience, singletons and unflagged reads had been merged and resulting fastq files receive the label "allSE" (commands not shown).

 
### Adapter Content Estimation

Adapter sequences had been summarized in a multi-fasta file (with header), and additionally written to a text file (without header) for pattern matching automation. The files are provided [elsewhere] (LINK).
The script ``rg_adapters.sh`` prints read counts and exact adapter matches per fastq file to a tab-separated table. Columns correspond to the adapter (and read counts) and lines to the fastq file. The script is executed as follows:

```
INDIR="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA"
AD_BN="adapter_seqs"
OUT="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/adapter_counts.tsv"
<path-to-script>/rg_adapters.sh "${INDIR}" "${AD_BN}" "${OUT}"
```


### Quality Enrichment

Quality enrichment was performed using different filter and trimming options provided by the preprocessing tool fastp. Those include:

1. minimum phred quality of 24 to qualify a base.
2. mean phred quality of 24 in the sliding window to qualify the bases.
3. the sliding window was set to size 4.
4. sliding window is moved over the read from both ends, 3' and 5'.
5. maximum N content of 5 to retain a read.
6. maximum 25% unqualified bases per read (lowered compared to default: 40%)
7. enabling minimum complexity per read of 30 % (default)
8. enabling trimming of polyX (X = {A,C,T,G}) tails that are 10 nt long or more 
9. enabling base correction of overlapping regions in paired-end modus (minimal overlap of 30 nt with a maximum of 5 mismatches: default)  
10. adapter trimming, provided a multi-fasta with adapter sequences
11. minimum read length of 20 nt (after trimming)

(The reasoning of the chosen settings can be read in the thesis' methods (LINK).)

The script ``trimming.sh`` conducts paired-end trimming (for intact pairs) and single-end trimming (singletons + merged) separately. Output fastq files are generated in accordance to input file names + the "trimmed" name extension. Resulting singletons are kept and quality-failed reads are written to an individual directory. This leads to a total of 5 qualified (allSE, pair1, pair2, unpair1, unpair2) and two unqualified multi-fasta files per input bam. After trimming, each file quality is reported using fastqc. 
The following command is demanded:

```
INDIR="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA"
OUT="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/trimming"
QCDIR="/data/mahlzeitlocal/projects/ma_neander_assembly/anc_virus_MA/qual_reports/trimmed"

<path_to_script>/trimming.sh "${INDIR}" "${OUT}" "${QCDIR}"
```
Note, that all unpaired trimmed pair1 or pair2 reads had been move to the quality failed directory since they did not meet quality criteria in manual check up.

## K-mer based Read Classification

The script ``clark_classification.sh`` wraps all steps for the k-mer classification with CLARK. Thereby it performs the following steps:
1. Database download (set_targets.sh)
2. Database creation, calculation of k-mer spaces and classification (classify_metagenome.sh)
3. Abundance estimation and conversion into KRONA output format (estimate_abundance.sh)
4. Abundance visualization in multi-layered pie charts (ktImportTaxonomy)

Steps 2. - 4. had been applied to every quality enriched fastq file. Thereby, Paired-end and single-end files had been conducted separately, leading to two separated classification files.
For all three k-mers tested, the confidence threshold was set to the lowest possible of 0.5. The value is simply calculated by $\frac{hitcount1 + hitcount2}{hitcount1}$ (SOURCE). A low confidence score therefore might lead to lower sensitivity and thus to classification of database related unknown organisms (possibly ancient).
The script can be applied as follows:
```
<path_to_script>/clark_classification.sh
```
``ktImportTaxonomy`` generates the pie charts in interactive file of html format that can be visited using standard web browsers (firefox, chrome, vivaldi).  

## Read Binning

The read extraction and binning consists of three steps:
1. Assign the phylogenetic path (scientific names) to the first assignment (NCBI TaxID) of each read (``TaxID_sci_names.py``)
2. Extract the viral, bacterial or non-assigned reads from quality enriched multip-fastq-files (``extract_reads_dicts.py``). This leads to three files per trimmed fastq file (e.g. SN7_1_1__allSE_trimmed.fastq --> SN7_1_1__allSE_viral.fastq, SN7_1_1__allSE_NA.fastq,SN7_1_1__allSE_bacteria.fastq)
3. Pool the extracted reads of same taxonomic branch (viral, NA, bacteria) and read type (allSE, pair1, pair2) (``fromdict_extract_reads.sh``). This reduces the overall file number to 9 in total for the subsequent analyses.

The ``TaxID_sci_names.py`` translates NCBI taxonomy IDs (TaxIDs) into the pyhlogenetic paths of the corresponding scientific names. It utilizes the ``names.dmp`` and ``nodes.dmp`` files from the NCBI taxonomy structure [taxdumb.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz). The functions ``get_dicts`` and ``get_tax_path`` were adapted from a custom python script provided by Florian Mock (Thanks!). Per classification output file, three text files had been generated collecting reads with a viral, NA, or bacterial first assignment.

```
# iterate over all output csv files and get phylogenetic path
for file in <path_to_classify_metagenome.sh_outputdir>/*.csv; do echo ${file}; ./scripts/TaxID_sci_names.py ${file}; done
```

Generated text files share the structure:
```
<read_header_ID>,<TaxID>,<taxon_path>

# for example:
SN928_0068_BB022WACXX:1:1103:19311:49395,1982588,Edwardsiella virus KF1 > Kafunavirus > Podoviridae > Caudovirales > Viruses
```

Subsequently, ``extract_reads_dicts.py`` and ``fromdict_extract_reads.sh`` use the header IDs from those text files to extract the corresponding entries from quality enriched multi-fastq files.
Note, that the shell script executes the python script for the per file extraction.

```
# all paths are specified in the file
<path-to-script>/fromdict_extract_reads.sh
```

### GC filtering

Complexity filtering of fastp refers to the base change frequency in a read at position i and i+1. Therefore, it might fail to filter overrepresented low GC reads due to repetitive or non-repetitive sequences of G and C. Therefore, GC filtering was adressed manually using the script ``trim_GC.py``. The approach simply counts G and C occurences in a read and divides it by the read length. Reads are qualified due to a customizable GC threshold. Script usage is:

```
 <path-to-script>/trim_GC.py -i <input-fastq> -o <passed-fastq> -f <failed-fastq> -g <threshold>
```

The script was applied on pair1, pair2 and allSE fastq files of viral, NA and bacteria read bin. The threshold was set to 25 and 45 filtering the entire reads and the >25% GC reads, respectively. For the latter, the quality failed fastqs comprised the 25% - 45% read set.
In case of paired-end sequences, singletons might arise due to diverging GC content in mate pairs. For that, output files with mate pairs had been compared on their included header IDs. If no match could be found, reads were written to a temporary singleton file and merged with the corresponding allSE file. This is comprised in the script ``filter_pairs.py``:

```
<path-to-script>/filter_pairs.py -i1 <in_pair1> -i2 <in_pair2> -o1 <out_pair1> -o2 <out_pair1> -s <singletons>  

# merge with allSE
cat  <singletons> >> <allSE.fastq> 
```
All corresponding pair fastq files, previously treated with ``trim_GC.py`` had undergone this filtering.

## Assembly Evaluation: K-mer Sizes and Read Correction

The assemblies had been generated using the ``cond_test.sh`` script. This script comprises all SPAdes assembly commands and is run as follows:

```
<path-to-script>/cond_test.sh
```

Secondly, HiSat2 was executed on all generated assemblies to remap the original reads. For read corrected assembly runs, the corrected reads had been used. Additional had been calculated using the ``HiSat2_cond_test.sh``. Those are contig length and (un)mapped read counts (``samtools idxstats``), total remapping rates (``samtools flagstat``), and per base depth (``samtools depth``). 

```
<path-to-script>/HiSat2_cond_test.sh
```

The remapping rates can be read from the HiSat2 log file and the flagstats output.
The average depth per contig was computed from the sum of the per base depth in a contig divided by its length. The information were read out from ``samtools depth`` output using the ``seq_depth.py`` script:

```
<path-to-script>/seq_depth.py <samtools_depth_output.tsv>
```
The output file is named extending the basename of the input file with "_coverage.tsv".

Third, the contigs had been BLASTed on viral NCBI genome sequences, using the ``BLAST.sh``. For information on the previously build blastdb consider ``makeblastdb.txt``. The script can be executed like this:

```
# move to blastdb directory so that blast can assign the TaxIDs properly to scientific names
cd <targetdb_dir>/taxonomy
# run blast
<path-to-script>/BLAST.sh   <query.fastq> <blastdb-basename> <out_csv>

# only the best hit (highest bit score) for a contig was retained:
sort -u -k1,1  <raw_blastout.csv> | sort -rnk2,2 > <tmp.csv>; cat <tmp.csv> > <blastout_unique.csv>; rm <tmp.csv> 
```

The output table contains a comma-separated list in the format specified in the script. Only queries with target hits are reported. From this, the fraction of contigs with target hits (of different similarity degree, of course) was calculated divided . Note, that blastn was used in megablast modus which is designed for very similar sequences. This is fine for now as we are rather interested in comparing the hit count fractions of different assemblies. For the qualitative search this approach is not used (see below LINK) and replaced with an MMSEQS approach.


## Assembly Evaluation: Input Library and Assembler Subtype

The input file configuration and different assembly algorithms had been tested. The assembly commands are comprised in ``test_lib.sh``:

```
<path-to-script>/test_lib.sh
```

HiSat2 remapping was performed with ``HiSat2_lib.sh``:
```
<path-to-script>/HiSat2_lib.sh
```

Blast search and contig depth calculation was performed as described in the [previous chapter](#assembly-evaluation-k-mer-sizes-and-read-correction).


## Assembly Evaluation: Contiguity and GC Filtering

Assemblies from read sets treated with different GC filters had been created. The assembly commands are comprised in ``GCfilter_test.sh``:

```
<path-to-script>/GCfilter_test.sh
```

HiSat2 remapping was performed in ``HiSat2_GCfilter.sh``

```
<path-to-script>/HiSat2_GCfilter.sh
```

Blast search and contig depth calculation was performed as described [before](#assembly-evaluation-k-mer-sizes-and-read-correction).

## Final Assembly Set Up

The viral and NA read assemblies had been chosen from the [Contiguity and GC Filtering Evaluation](#assembly-evaluation-contiguity-and-gc-filtering). The missing bacterial assembly was generated with the following command:

```
INPATH=<path-to-reads_bigger25%GC>
OUT=<outpath>

metaspades.py --only-assembler -m 1500 -t 40 -1 <${INPATH}/all_bacteria_pair1_trimmed_filterGC.fastq> -2 "${INPATH}/all_bacteria_pair2_trimmed_filterGC.fastq"  --merged "${INPATH}/all_bacteria_allSE_trimmed_filterGC.fastq" -k 15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123,127 -o "${OUT}/" &> "${OUT}/bacteria/bigger25/bacteria_bigger25.log" 

# metaSPAdes Parameter:
#   --only-assembler    runs only the assembler, read correction disabled
#   -t  number of threads to use
#   -1  specify pair1 input reads
#   -2  specify pair2 input reads
#   --merged    specify merged input reads
#   -k  list of k-mers to use for assembly graph construction
#   -o  output directory
```

The read correction was disabled due to exhaustive memory consumption assinging a maximum of available RAM. Input reads were filtered on > 25% GC contents.
The remapping was performed during the comparative [remapping evaluation](#remapping-evaluation) on all three bins. Bwa mem results had been considered for downstream analyses for all three data sets. Homology searches had been performed with mmseqs2, described in the [section below](#contig-homology-search-on-comprehensive-databases) 


## Assembly Evaluation: Other Assemblers

Megahit was executed on the viral and NA read set using equivalent conditions as for metaSPAdes (PE + merged reads > 25% GC, k-mers set to list, remapping with bwa mem). The assembly was executed with ``megahit.sh``:

```
<path-to-script>/megahit.sh
```
Remapping with bwa mem and flagstats are comprised in ``megahit_remapping.sh``:

```
<path-to-script>/megahit_remapping.sh
```

Remapped read counts per contig, contig length and per base depth was calculated with ``megahit_remapping_stats.sh``:

```
<path-to-script>/megahit_remapping_stats.sh
```

Blast search and contig depth calculation was performed as described [before](#assembly-evaluation-k-mer-sizes-and-read-correction).


## Remapping Evaluation

Different alignment tools had been tested on the viral and NA read sets using the chosen conditions from the assembly evaluations. The commands are summarized in ``remapping.sh``. This script also includes samtools flagstat, idxstats and depth calculations.

```
<path-to-script>/remapping.sh
```
The remapping statistics had been read from the ``samtools flagstat`` output.
For bwa mem and segemehl, paired-end reads and singletons/merged pairs had been remapped separately due to available options. Therefore, the remapping rates had to be summarized afterwards as well as the contig depths.

The latter had been calculated as described [before](#assembly-evaluation-k-mer-sizes-and-read-correction) for PE and SE reads individuallyand merged with a command pipeline using ``join`` and ``awk``:

```
join -1 2 -2 2 <SE_coverage.tsv> <PE_coverage.tsv> | awk '{print $1 "\t"  $3+$5}' > <merged_coverage.tsv>

```
Susequently, contig lengths extracted with ``samtools idxstats`` and the average contig depths had been merged into a tab separated file using the script ``contig_deplen.py``:

```
<path-to-script>/contig_deplen.py <idxstat_output.tsv> <merged_coverage.tsv> <merged.tsv>
```

## Contig Homology Search on Comprehensive Databases

The contigs had been filtered requiring minimum 300 nt length and 10x average depth. A text file denoting the qualified contig identifiers had been generated from the merged length and depth table:

```
awk '($2>=10 && $3>=300)' <path-to-bin-assembly>/final_contigs_cov_length.tsv > <path-to-bin-assembly>/contig_names_cov10_len300.txt
```
The script ``filter_contigs.py`` writes the subset of fasta sequences to a multi fasta file:

```
<path-to-script>/filter_contigs.py <path-to-bin-assembly>/contig_names_cov10_len300.txt <path-to-bin-assembly>/contigs_len300_cov10.fasta
```
This fasta was adressed to homology search with Mmseqs2. 

Viral, bacterial, fungal and human sequences had been downloaded from different sources:

```
# Viral Sequences were obtained from a metagenomics pipeline provided by Dr. Martin Hoelzer (database still the latest version)
wget -nH ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/IMG_VR_2018-07-01_4.tar.gz && tar zxvf IMG_VR_2018-07-01_4.tar.gz && tar zxvf <viral_dir>/IMG_VR_2018-07-01_4.tar.gz

# GTDB bacterial/archaeal representatives obtained from the GTDB ftp server:
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/gtdb_rep_genomes.tar.gz
tar -zxf gtdb_rep_genomes.tar.gz && gunzip <bacteria_dir>/*.gz && gunzip <archaea_dir>/*.gz

# Download latest patch of human genome assembly (GRCh38.p13)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz -P <human_dir>/human -o <human_dir>/human_wget.log

# Download all RefSeq fungi sequences
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.*.genomic.fna.gz -P <fungi_dir>/fungi -o <fungi_dir>/fungi_wget.log

```
All sequences had been merged into one fasta file:
```
# GTDB bacteria - all representatives
for file in <bacteria_dir>/*.fna; do cat ${file} >> <targetdb_dir>/all_sequences.fna; done
# Archea
for file in <archaea_dir>/archaea/*.fna; do cat ${file} >> <targetdb_dir>/all_sequences.fna; done
# Fungi
for file in <fungi_dir>/fungi/*.fna; do cat ${file} >> <targetdb_dir>/all_sequences.fna; done
# Human
cat <human_dir>/GCA_000001405.28_GRCh38.p13_genomic.fna >> <targetdb_dir>/all_sequences.fna
# Virus 
cat <virus_dir>/IMGVR_all_nucleotides.fna >> <targetdb_dir>/all_sequences.fna
```

The homology search requires five steps: 
1. Query database creation
2. Target database creation
3. Target database indexing
4. Homology search
5. Results conversion into tsv format

Those steps are summarized in the script ``mmseqs2.sh``:

```
<path-to-script>/mmseqs2.sh <query.fasta> <query-type>

# <query-type>: "bacteria", "viral", "NA"
```

Per query only the highest target hit was retained. This function is executed with ``filter_mmseqs2.py``. It extracts all target hits per query, sorts it descending by bit score and writes the highest scoring hit per query in the blast-like tsv format.

```
<path-to-script>/filter_mmseqs2.py <mmseqs_out.tsv> <mmseqs_out_unique.tsv>
```

## Contig Clustering

Clustering with MMseqs requires several commands comprised in the script ``mmseqs_cluster.sh``. Those cover the following steps:

1. Query Database creation
2. Clustering 
3. Clustering results conversion into tab separated format.
4. Output cluster sequence information in fasta format.
5. Cluster representative extraction to a multi-fasta file.

```
<path-to-script>/mmseqs_cluster.sh <query.fasta> <querydb_dir> <prefix_querydb> <resultsdb_dir> <tmp_dir> 5 0.8 0.001 
```

To calculate the cluster sizes, a header file had been extracted from each cluster multi fasta file. 
```
awk '/>/ {print $0}' <path-to-cluster-results>/seqfiledb.fasta > <path-to-cluster-results>/headers.txt
```
The format indicates each new cluster by a single repetition of an header ID which is the representative of the cluster. Two repeated header IDs are filled with other headers belonging to one cluster. From this information, the cluster sizes was calculated with the script ``clustersizes.py``.

```
<path-to-script>/clustersizes.py <header_file.txt> <output.tsv>
```
The output is a two column tab separated file (representative, cluster size) where each line correspond to a new cluster. The file is sorted for descending clustersize.

The cluster affiliation and homology information per contig had been comprised into a summarizing table. In this format, each line shows information on a particular contig, comprising contig length, average read depth and mmseqs target hit metrices (if no hits found, indicated as 'None'). Contig entries are intermitted with representative header IDs which indicate the start of a new cluster. This representative is indicated by a leading '>' character and gives further the information of affiliated contig number in the second column.

Those files were generated with the ``contig_summary.py`` algorithm. The script takes as input the unique mmseqs hit file, the contig length, depth tsv and the cluster header txt:

```
<path-to-script>/contig_summary.py -cl <path-to-assembly-dir>final_contigs_cov_length.csv -ci <path-to-cluster-results>/header.txt -m <path-to-mmseqs-results>/contigs_len300_cov10_unique.tsv -o <outpath>/contig_cluster_hits.tsv
```

## Clustering on Target Hit IDs

MMseqs2 search results of viral, NA and bacterial read bin contigs had been clustered on hit target IDs. This was performed on the target hits with highest bit score ("unique" MMseqs2 results). The script ``cluster_taxes.py`` implements this approach going once through the list of all query records in the unique MMseqs2 output tsv. It checks for first occurence of the corresponding target ID or appends the its cluster. Thereby, the number of contigs per cluster is counted and the cumulative contig length calculated:

```
<path-to-script>/cluster_taxes.py <mmseqs_unique.tsv> <targetID_clusters.tsv>
# Output format is:
# targetID    contig_counts     target_sciname  cum_contig_len  mean_contig_len
```

This list is already sorted for descending cumulative contig lengths of the hits.

## Selection of Ancient Candidate Reference Sequences

The 1000 target ID clusters with highest cumulative contig lengths had been appended on average homology search measures of belonging contigs. To achieve this, the script ``top_hitstats.py`` calculates these mean values from the combining contig clustering and MMseqs2 hit tsv and the target ID clusters. A third parameter defines the number of clusters to process going through the list linewise (sorting was performed, therefore). The output is generated in tsv format, too:

```
<path-to-script>/top_hitstats.py -i <targetID_clusters.tsv> -m <contig_cluster_hits.tsv> -n 1000 -o <top_hit_stats.tsv>
```

The upper 1000 target cluster IDs from viral, NA and bacterial homology searches had been consulted for manual selection of promising ancient candidate reference sequences. From those, the target hits with highest cumulative contig lengths had been selected. Further criterions had been considered for candidate selection:
1. For high query coverage (from roughly 70%), the sequence identity was moderate (~ 70 - 80 %).
2. In case of the viral target ID clustering, from rare viral target clusters, the upper most three targets had been selected.

Of note, the candidate selection was done manually on the given criteria. No further systematic filtering could be applied. The candidate lists therefore are a subjective selection. 

The candidate target IDs had been written to a text file to guide the extraction of their corresponding fasta sequences from the bulk database fasta (created for MMseqs search, explained [here](#contig-homology-search-on-comprehensive-databases)). The script ``sub_target_extraction.py`` serves for this purpose:

```
<path-to-script>/sub_target_extraction.sh <candidate_IDs.txt>
```

Each identified fasta record is written to an individual fasta file. The output path was changed for viral, NA, or bacterial read bin candidate references.

## Creating a Microbacterium Consensus Sequence

To demonstrate the consensus approach, 30 microbacterium target IDs with highest cumulative contig lengths had been extracted from the bacterial target ID clustering. The following command (rip)greps for the key word 'Microbacterium' which might occur on the 'scientific name' column. Doing this, the sequence space for the consensus construction is restricted to the microbacterium genus.

```
rg 'Microbacterium' <targetID_clusters.tsv>  > <path-to-consensus-dir>/microbacteirum_headers.txt | head -n 30 | cut -f1 > <path-to-consensus-dir>/microbac_30_IDs.txt
```

Cactus requires phylogenetic pre-knowledge of the input sequences. To approximate the microbacterium sequence relation 16S rRNA sequences had been utilized instead of the entire sequence records. Those had been obtained from the GTDB taxonomy which is based on the 16S rRNAs of taxonomically incorporated database sequences.

```
# get the 16S rRNA sequences from GTDB taxonomy
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ssu.fna

# extract the 16S header sequences that correspond to the microbacterial genomes
rg -f microbac_30_IDs.txt ssu.fna > microbac_16S_headers.txt
```
Only primary hits were considered, for example:
```
# hits on the same microbacterial ID:
>RS_GCF_001314225.1~NZ_CP012697.1           -> this was considered for tree construction
>RS_GCF_001314225.1~NZ_CP012697.1-#2
```

From 30 microbacterial IDs, 10 showed a corresponding 16S rRNA within the taxonomy. The fasta records of those had been extracted using ``sub_target_extraction.py`` (bulk fasta path was changed to the ``ssu.fna``). The rRNA fasta sequences had been merged into one multi-fasta file. From this, a MSA was generated using MAFFT g-ins-i. This command is recommended for sequences of similar length. MAFFT generally aligns sequences on their fast Fourier transformation on volume and polarity range. This approach makes it faster than clustalW, for example.

```
# merge rRNA fastas
for seq in <path-to-fastas>/*.fasta; do cat ${seq} >> <path-to-fastas>/all_microbac_ssu.fasta; done

# create MSA
nice ginsi --thread 40 <path-to-fastas>/all_microbac_ssu.fasta > <path-to-fastas>/all_microbac_ssu_msa.fasta
```

From this, a phylogenetic tree was constructed with fasttree based on maximum likelihood estimations of nearest-neighbor interchanges and the minimum-evolution criterion:

```
fasttree <path-to-fastas>/all_microbac_ssu_msa.fas > <path-to-tree>/all_microbac_ssu_msa.tree
```
The full microbacterium sequences with identified rRNA records had been extracted to individual fasta files using the ``sub_target_extraction.py``. The header had been truncated to consist only of the target IDs:

```
sed -i 's/ Microbacterium.*//' <microbacterium.fasta
```

The phylogeny is structured in newick format. For cactus input, the rRNA header IDs in the pyhlogeny had to be trimmed to the exact microbacterial genome header. Beneath the altered newick tree, the modified IDs are matched to the corresponding microbaterium fasta path. The final input file looks like this:

```
(BDCY01000002.1:0.003203152,NZ_CP025422.1:0.005458557,((NZ_JHET01000001.1:0.000000005,NZ_JHET01000006.1:0.000000005)0.997:0.012483338,(NZ_PDJE01000001.1:0.052728170,(NZ_LMGU01000001.1:0.002956147,(NZ_FNRY01000001.1:0.029714099,(NZ_CP012697.1:0.020212651,(NZ_CP019892.1:0.011022607,NZ_LT629770.1:0.001496155)0.993:0.010784790)0.988:0.022530779)1.000:1.243046039)0.931:0.016352431)0.962:0.009052175)0.759:0.004973046);

BDCY01000002.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/BDCY01000002.1.fasta 
NZ_CP025422.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_CP025422.1.fasta 
NZ_JHET01000001.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_JHET01000001.1.fasta
NZ_JHET01000006.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_JHET01000006.1.fasta
NZ_PDJE01000001.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_PDJE01000001.1.fasta
NZ_LMGU01000001.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_LMGU01000001.1.fasta
NZ_FNRY01000001.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_FNRY01000001.1.fasta
NZ_CP012697.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_CP012697.1.fasta
NZ_CP019892.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_CP019892.1.fasta
NZ_LT629770.1 /data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/NZ_LT629770.1.fasta
```

Cactus was installed from precompiled binaries following the instructions in the [github repo](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v1.0.0/BIN-INSTALL.md). The [latest cactus release](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v1.0.0/BIN-INSTALL.md) was utilized.
Cactus was set up to run in a virtual environment:

```
source venv/bin/activate
cactus jobStore <path-to-consensus>/microbacteria_input.txt <path-to-consensus>/10microbacteria_msa.hal --maxCores 30 --maxMemory 500G --defaultMemory 15G 

# cactus parameters:
#   --maxCores      number of threads to use
#   --maxMemory     memory limit to be assigned by cactus
#   --defaultMemory default memory per sub-task
```

The output in `hal` format is efficient but not human readable. The ``hal2maf`` tool comes along with cactus and was used to convert the alignment in human readable multiple alignment format (MAF). 

```
hal2maf <path-to-consensus>/10microbacteria_msa.hal <path-to-consensus>/10microbacteria_msa.maf
```

With default settings, ``hal2maf`` writes the alignment referenced to root. This root sequence can be considered as consensus. The sequence record was parsed to a fasta file with ``maf_to_cons.py``. 

```
<path-to-consensus>/maf_to_cons.py <path-to-consensus>/10microbacteria_msa.maf <path-to-consensus>/consensus.fasta
```

## Damage Pattern Analysis

Extracted candidate sequences and the generated consensus served as reference sequences in read mappings using ``bwa aln`` and ``bwa samse``:

```
# Indexing References
nice bwa index -p <ref_index_prefix> <ref.fasta> &> <indexing.log>

# Mapping with bwa aln
nice bwa aln -t 30 -n 0.01 -o 2 -l 16500 <ref_index_prefix> <single_end_reads.fasta> > <alignment.sai> 2> <alignment.log>

# Converting from index to SAM
nice bwa samse <ref_index_prefix> <alignment.sai> > <alignment.sam>

# BWA Parameters:
#   -p  prefix of sequence database
#   -n  fraction of missing alignments giving a 2% base error rate
#   -t  threads to use
#   -o  number of allowed gaps
#   -l  seeding length 
```
The mapping rates had been generated with ``samtools flagstat`` as described before.

Each candidate was selected from an target ID clustering that originated from the homology search of one set of read bin contigs (viral, NA, bacteria). Those read bins had only been mapped to candidates from the corresponding clustering. 
From remappings, substitution frequencies and base damage probability estimations were calculated using the ``mapdamage2`` package:

```
# PE mapping
mapDamage -i <alignment_PE.sam> -r <ref.fasta> -d <output_dir> -l 95

# allSE mapping:
mapDamage -i <alignment_SE.sam> -r <ref.fasta> -d <output_dir> -l 43

# Mapdamage Parameters:
#   -i  input mapping file in SAM format
#   -r  used reference sequence in fasta format 
#   -d  output directory for all output files; generated if not present
#   -l  length of the reads to assess
```

The read length was selected from the quality reports of the > 25% GC conntaining fastq files from the viral, NA, and bacterial bin.



# Visualization with Matplotlib

Charts and plots had been generated with functions from the ``matplotlib`` python library. For that, the data had been parsed into suitable text-based formatting to generate a graphical output. Some figure captions had been modified using the inkscape (v0.92.4) software for better overview. The method summary flowcharts had been generated on the free ``drawio`` online platform.