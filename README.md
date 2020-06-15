# MA_Methods
This file shall supplement the Master's Thesis ``ENTER TITLE`` with detailed descriptions of all relevant commands and custom scripts. It is provided in addition and structured in accordance to the Methods section which I recommend to conduct in case of uncertainties.
The parameter functionality of utilized tools is described within the custom scripts, if present. Custom scritps implementing command line interfaces (CLI) are explained in the text.

## Summary Workflow Description
- TO DESCRIBE: How were the workflow graphics prepared?
- 
## Data Preprocessing
The input data had been downloaded using a recursive wget command:
```
cd <Data_Directory>
wget -r -np -nd http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/unmapped_qualfail/

Parameters:
-r  recursive download
-np exclude parent directories from download
-nd do not create directories on local drive
```
This command downloads all files present in the specified directory, in this case 33 bam files.

### Extraction of Input Reads

The script `fastq_pairedend.sh` extracts pair1, pair2, singleton, as well as both or none of pair1 and pair2 flagged reads from each bam file. This creates 4 fastq files per bam. In addition, it creates a fastqc report in html format for all fastq files.

```
cd /data/fass2/reads/max_crass/homo_sapiens_neanderthalensis_altai_unmapped
/data/mahlzeitlocal/projects/ma_neander_assembly/anc_virus_MA/scripts/fastq_pairedend.sh 
```

