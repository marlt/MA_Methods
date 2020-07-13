#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6

# Extract consensus sequence from a MAF file containing alignments that are referenced to root

import sys
from collections import OrderedDict

anchor_pos = OrderedDict()  # the alignment is splitted in denoted "chromosomes"; save the length of each of this subregions
anchor_seq = OrderedDict()  # save the sequence of each chromosome as list - alignment is broken up into several alignment blocks

maf_file = sys.argv[1]      # alignemnt if MAF format
headline = sys.argv[2]      # string to write as header of the fasta entry 
outfile = sys.argv[3]       # path for output
count = 0                   # counter for the sequences in an alignment block; the root/consensus sequence occurs for count == 1

with open(maf_file, 'r') as infile:
    for line in infile:
        count += 1                          # increment for each new line
        if line.startswith('a'):            # indicates the start of a new alignment block             
            count = 0                       # then, counter is incremented
        elif count==1 and line.startswith('s'): # if the first line in an alignment block is a sequence, it is the root/consensus
            line_ = line.rstrip('\n')
            entry = line_.split('\t')
            length = entry[5]                   # length of the alignment block
            key = anchor_pos.get(entry[1])      # the alignment is splitted
            if key == None:                     # if chromosome name is not in dictionary, initialise:

                anchor_pos.update({entry[1]:length})        # enter the chromosome and its length
                anchor_seq.update({entry[1]:[entry[6]]})    # save the sequence record of the alignment block

            else:
                anchor_seq[entry[1]].append(entry[6])       # if present, append the seuqence with the sequence from the encountered alignment block

        else:
            pass

sum_contig_lengths = 0

for key in anchor_pos.keys():
    # print(anchor_pos[key])
    sum_contig_lengths += int(anchor_pos[key])  # add up chromosomes to the entire genome length

print(f'Summed contig length columns:\n{sum_contig_lengths}')

consensus = ''          # save the consensus here
N_count = 0             # count 'N' nucleotide occurences in the consensus to estimate the quality

for contig in anchor_seq.keys():            # add up the sequence stretches to the full genome
    for subseq in anchor_seq[contig]:   
        subseq_upper = subseq.upper()       # make the m all upper case letters, no differences are recognized
        N_count += subseq_upper.count('N')  # count N 
        consensus += subseq_upper           # add up sequence stretches

with open(outfile, 'w') as out:
    out.write(f'>{headline}\n')               # write header
    out.write(consensus)                        # write sequence              

print(f'Total base counts of concat consensus:\n{len(consensus)}')              # this number shall be the same as the sum_contig_lengths
print(f'N-count in sequences:\t{N_count}\t{round(N_count*100/len(consensus),4)} %')