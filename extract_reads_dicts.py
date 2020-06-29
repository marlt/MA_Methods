#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3

# extract fasstq entries from a multi-fastq file on the basis of provided text-files (TaxID_sci_names.py) containing headers to look for
# the script loads all fastq records in a file to RAM in a dictionary and searches in it; depending on fastq-file size and available RAM this approach is limited, but also fast

# PACKAGES
import sys, os

# GLOBALS
# readID to phylogenetic path mappings
fastq_path = sys.argv[1] # fastq file to search in: this is completely loaded to RAM
vir_path = sys.argv[2]  # viral txt file
NA_path = sys.argv[3]   # NA header txt file
bac_path = sys.argv[4]  # bacterial header txt file
outpath = sys.argv[5]   # outpath to write to

BN = os.path.basename(fastq_path)[:-14]     # basename for the extracted reads - where are the reads from is tracable
print(BN)
fastq_dict = {} # save all the fastq entries here with header as key

# load fastq data:
print(f'Load {BN}_trimmed.fastq to dictionary...')
with open(fastq_path, 'r') as fastq_f:
    #count = 0
    for line in fastq_f:
        #count += 1
        #if count <= 10:                         # for debugging: check the first 10 lines
        header = line[1:].rstrip('\n')                                   # read every forth line as header, crop the '@'
        seq = fastq_f.readline().rstrip('\n')                             # second line is sequence
        plus = fastq_f.readline().rstrip('\n')                            # third line is additional context
        phred = fastq_f.readline().rstrip('\n')                           # forth line is quality string
        fastq_dict.update({header:[seq,plus,phred]})        # add the whole fastq entry to the fastq_dict
        # print(fastq_dict)
        #else:
            #break

# create file to write viral reads:
with open(f'{outpath}/{BN}_viral_extract.fastq', 'w') as vir_reads:
    pass
# create file to write NA reads:
with open(f'{outpath}/{BN}_NA_extract.fastq', 'w') as NA_reads:
    pass
# create file to write bacterial reads:
with open(f'{outpath}/{BN}_bacteria_extract.fastq', 'w') as bac_reads:
    pass


print('Extract viral reads...')
with open(vir_path, 'r') as vir_IDs:
    with open(f'{outpath}/{BN}_viral_extract.fastq', 'a') as vir_reads:
        for header in vir_IDs:
            try:
                header_ = header.rstrip('\n')
                entry = fastq_dict.get(header_)          # find header in O(1)
                #print(entry)                    
                vir_reads.write(f'@{header_}\n')           # write: header
                vir_reads.write(f'{entry[0]}\n')         # seq
                vir_reads.write(f'{entry[1]}\n')         # plus
                vir_reads.write(f'{entry[2]}\n')         # qual string
            except KeyError:
                print(f'{header} not found.')           # all entries should be found! Otherwise the strings are not correctly processed...

# same for NA
print('Extract NA reads...')
with open(NA_path, 'r') as NA_IDs:
    with open(f'{outpath}/{BN}_NA_extract.fastq', 'a') as NA_reads:
        for header in NA_IDs:
            try:
                header_ = header.rstrip('\n')
                entry = fastq_dict.get(header_)          
                #print(entry)                    
                NA_reads.write(f'@{header_}\n')           
                NA_reads.write(f'{entry[0]}\n')         
                NA_reads.write(f'{entry[1]}\n')         
                NA_reads.write(f'{entry[2]}\n')         
            except KeyError:
                print(f'{header} not found.')          

# same for bacteria
print('Extract bacterial reads...')
with open(bac_path, 'r') as bac_IDs:
    with open(f'{outpath}/{BN}_bacteria_extract.fastq', 'a') as bac_reads:
        for header in bac_IDs:
            try:
                header_ = header.rstrip('\n')
                entry = fastq_dict.get(header_)
                #print(entry)
                bac_reads.write(f'@{header_}\n')
                bac_reads.write(f'{entry[0]}\n')
                bac_reads.write(f'{entry[1]}\n')
                bac_reads.write(f'{entry[2]}\n')
            except KeyError:
                print(f'{header} not found.')