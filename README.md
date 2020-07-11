# MA_Methods
This file shall supplement the Master's Thesis ``ENTER TITLE`` with detailed descriptions of all relevant commands and custom scripts. The structure is oriented on the Methods section which I recommend to conduct in case of uncertainties. Additional chapters are introduced if this is necessary to keep the pipeline character of the script explanations.
The parameter functionality of utilized tools is described within the custom scripts, if present. Custom scripts implementing command line interfaces (CLI) are explained in place. If not described further, default tool settings apply.

If not described any further, programs had been used in a conda environment for better version control and a conflict free dependency management. 

- ADD AN TABEL OF CONTENTS
- 
## Summary Workflow Description
- TO DESCRIBE: How were the workflow graphics prepared?
- 
## Data Preprocessing
The input data had been downloaded using a recursive ``wget`` command:
```
cd <Data_Directory>
wget -r -np -nd http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/unmapped_qualfail/

Parameters:
-r  recursive download
-np exclude parent directories from download
-nd do not create directories on local drive
```
This command downloads all files present in the specified directory, in this case 33 bam files. Since the read names were unconveniently long, files had been given shorter, unique names (commands not shown). A filename mapping is provided here (LINK).

### Extraction of Input Reads

The script `fastq_pairedend.sh` includes the following steps: sorting bam files by read names (``samtools sort``), calculating read mapping statistics (``samtools flagstat``) and extracting the read sequences (``samtools fastq``). The latter comprises reads flagged as pair1, pair2, singleton, or both or none of pair1 and pair2 in each bam file. This leads to 4 fastq files per input bam. In addition, it creates a fastqc report in html format assessing several quality measures for all fastq files.

```
BAM="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/"
FASTQ="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA"
QC="/data/mahlzeitlocal/projects/ma_neander_assembly/anc_virus_MA/qual_reports/raw"

<path_to_script>/fastq_pairedend.sh  "${BAM}" "${FASTQ}" "${QC}"
```

Fastq files labelled 'unflag' correspond to merged read pairs that were generated elsewhere and explained in VERWEIS AUF THEORY. For convenience, singletons and unflagged reads had been merged and resulting fastq files receive the label "allSE" (commands not shown).

 
### Adapter Content Estimation

Adapter sequences had been summarized in a multi-fasta file (with header), and additionally written to a text file (without header) for pattern matching automation. The files are provided elsewhere(LINK).
The script ``rg_adapters.sh`` prints read counts and exact adapter matches per fastq file to a tab-separated table. Columns correspond to the adapter (and read counts) and lines to the fastq file. The script is executed as follows:

```
INDIR="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA"
AD_BN="adapter_seqs"
OUT="/data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped/fastq_MA/adapter_counts.tsv"
<path_to_script>/rg_adapters.sh "${INDIR}" "${AD_BN}" "${OUT}"
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

Secondly, HiSat2 was executed on all generated assemblies to remap the original reads. For read corrected assembly runs, the corrected reads had been used. Additional had been calculated using the ``HiSat2_cond_test.sh``. Those are contig length and (un)mapped read counts (``samtools idxstats``), total remapping rates (``samtools flagstats``), and per base depth (``samtools depth``). 

```
<path-to-script>/HiSat2_cond_test.sh
```

The remapping rates can be read from the HiSat2 log file and the flagstats output.
The average depth per contig was computed from the sum of the per base depth in a contig divided by its length. The information were read out from ``samtools depth`` output using the ``seq_depth.py`` script:

```
<path-to-script>/seq_depth.py <samtools_depth_output.tsv>
```

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
The remapping was performed during the comparative [remapping evaluation](#remapping-evaluation). Bwa mem results had been considered for downstream analyses for all three data sets. Homology searches had been performed with mmseqs2, described in the [section below](#contig-homology-search-on-comprehensive-databases) 


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

Different alignment tools had been tested on the viral and NA read sets using the chosen conditions from the assembly evaluations. The commands are summarized in ``remapping.sh``:

```
<path-to-script>/remapping.sh
```

The remapping statistics had been read from the ``samtools flagstats`` output.
For bwa mem and segemehl, paired-end reads and singledtons/merged pairs had been remapped separately due to available options. The remapping rates had been summed afterwards. 


## Contig Homology Search on Comprehensive Databases

The viral, bacterial, fungal and human sequences had been downloaded from different sources:

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
1. query database creation
2. target database creation
3. target database indexing
4. hiomology search
5. results conversion into tsv format

Those steps are summarized in the script ``mmseqs2.sh``:

```
<path-to-script>/mmseqs2.sh <query.fasta> <query-type>

# <query-type>: "bacteria", "viral", "NA"
```

Per query only the highest target hit was retained. This function is executed with ``filter_mmseqs2.py``. It extracts all target hits per query, sorts it descending by bit score and writes the highest scoring hit per query in the blast-like tsv format.

```
<path-to-script>/filter_mmseqs2.py <mmseqs_out.tsv> <mmseqs_out_unique.tsv>
```

# Contig Clustering

- describe the single steps
```

<path-to-script>/mmseqs_cluster.sh <query.fasta> <querydb_dir> <prefix_querydb> <resultsdb_dir> <tmp_dir> 5 0.8 0.001 

```
- extract header file
- clustersize calculation
- cluster mmseqs assignment