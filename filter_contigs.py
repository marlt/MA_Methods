#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3

# script filters contigs based on provided textfile; this contains the contig IDs to retain


######## IMPORTS ###################
import os, sys

######### FUNCTIONS ################
def write_entry(head, sequence):
    global matches
    if matches == 0:                
        otype = 'w'                 # overwrites the file, if its the first contig to add
    else:
        otype = 'a'                 # append for all other contigs
    matches += 1
    with open(dname + "/contigs_len300_cov10.fasta", otype) as output:
        output.write(head)          # write header
        for line in sequence:       # write sequence
            output.write(line)

########## MAIN ####################

contig_IDs = []             # load contig IDs from text file here
seq = []
matches = 0
total = 0
txt_file = sys.argv[1]              # load the text file with the qualified contig headers
contig_file = sys.argv[2]           # fasta file containing all contigs

dname=os.path.dirname(contig_file)

with open(txt_file, 'r') as txt:
    for line in txt:
        # print(line)
        contig_IDs.append(">" + line)

with open(contig_file, 'r') as contigs:
    for line in contigs:
        if line.startswith('>'):            # is beginning of new contig
            # print('Encountered header')
            total += 1
            if seq != [] and line in contig_IDs:            # if the contig is in the qualified header list, write the contig to output 
                write_entry(header, seq)
            seq = []
            header = line
        else:                               # is sequence 
            seq.append(line)

print(f"Write {matches} of {total} sequences to {dname}/contigs_len300_cov10.fasta")

