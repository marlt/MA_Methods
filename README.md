# MA_Methods
This file shall supplement the Master's Thesis ``ENTER TITLE`` with detailed descriptions of all relevant commands and custom scripts. The structure is oriented on the Methods section which I recommend to conduct in case of uncertainties. Additional chapters are introduced if this is necessary to keep the pipeline character of the script explanations.
The parameter functionality of utilized tools is described within the custom scripts, if present. Custom scripts implementing command line interfaces (CLI) are explained in the text. If not described further, default tool settings apply.

The tools had been used in a conda environment for better version control. 
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

# K-mer based Read Classification

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

# Read Binning

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
Note, that the shell script executes the python script for the per file extraction prior to merging as described.

```
# all paths are specified in the file
<path-to-script>/fromdict_extract_reads.sh
```


