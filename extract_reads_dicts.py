#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3

# Script to extract fastq entries from a fastq file on the basis of provided text-files that contain the headers to check for
# PACKAGES
import sys, os

# GLOBALS
# readID to phylogenetic path mappings
fastq_path = sys.argv[1] # fastq file to extract the reads from
NA_path = sys.argv[2]   # NA header txt file
bac_path = sys.argv[3]  # bacterial header txt file
outpath = sys.argv[4]   # outpath to write to

BN = os.path.basename(fastq_path)[:-14]     # basename for the extracted reads - where are the reads from is tracable
print(BN)
fastq_dict = {} # save all the fastq entries here by header

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

# create file to write NA reads:
with open(f'{outpath}/{BN}_NA_extract.fastq', 'w') as NA_reads:
    pass
# create file to write bacterial reads:
with open(f'{outpath}/{BN}_bacteria_extract.fastq', 'w') as bac_reads:
    pass


print('Extract NA reads...')
with open(NA_path, 'r') as NA_IDs:
    with open(f'{outpath}/{BN}_NA_extract.fastq', 'a') as NA_reads:
        for header in NA_IDs:
            try:
                header_ = header.rstrip('\n')
                entry = fastq_dict.get(header_)          # find header in O(1), as Flo said
                #print(entry)                    
                NA_reads.write(f'@{header_}\n')           # write: header
                NA_reads.write(f'{entry[0]}\n')         # seq
                NA_reads.write(f'{entry[1]}\n')         # plus
                NA_reads.write(f'{entry[2]}\n')         # qual string
            except KeyError:
                print(f'{header} not found.')           # all entries should be found! Otherwise the strings are not correctly processed...

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



