#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6

# Cluster mmseqs targets in unique mmseqs search output on their target IDs (Accession number, from header sequences)


############### IMPORTS ######################
import argparse
import sys, os
from collections import OrderedDict

################ CLI #######################
parser = argparse.ArgumentParser()

# input file: mmseqs convertalis output tsv, uniquely filtered for highest bit score
parser.add_argument('-i', '--input', dest='infile', metavar='STR', required = True, help="Input file is an mmseqs convertalis tsv, filtered by highest bit score per query for this purpose.")
# tab separated output statistics
parser.add_argument('-o', '--output', dest='outfile', metavar='STR', required = True, help="Outfile to target ID, clustersize and taxon name")


args = parser.parse_args()

infile = args.infile
outfile = args.outfile

################ MAIN ########################
names = {}                                              # save in dict: {targetID:sciname}
counts = OrderedDict()                                  # save in dict: {targetID:occurrences}
lengths = {}                                            # save in dict: {targetID:concatenated contig lengths}
with open(infile, 'r') as ifile:
    for line in ifile:
        line_ = line.rstrip('\n')
        entry = line_.split('\t')                       # split every line on tab delimiter
        if entry[1] != "target":                        # if not headline
            count = counts.get(entry[1])                # check if there is already an entry for the encountered target ID
            if count == None:                           # if encountered ID for the first time, initialise:
                counts.update({entry[1]:1})                 # set counts of this target ID to 1
                names.update({entry[1]:entry[3]})           # save the taxonomic name of the target
                lengths.update({entry[1]:int(entry[8])})    # save contig length per target
            else:                                       # extend hit counts, concatenated contig length
                count += 1                              # increment for every encountered hit of same name
                # print(f'New contig\t{entry[8]}')
                length = int(lengths.get(entry[1]))          # get length entry of encountered target ID  
                # print(f'Length\t{length}')
                length = length + int(entry[8])         # add length of the new contig to the existing length
                lengths.update({entry[1]:length})       # update length entry for target ID  
                # print(length)
                counts.update({entry[1]:count})         # update count entry for target ID  

counts_sorted = OrderedDict(sorted(counts.items(), key=lambda x: x[1], reverse = True))         # sort for most encountered target ID

with open(outfile, 'w') as out:                                                                 # open output write stream
    total = 0
    tax_count = 0
    contig_length = 0
    out.write(f'header_ID\tHit_Counts\tTax_Name\tConcat_ContigLen\tContigLen_per_Hit\n')        # write header line for output file
    for hID in counts_sorted.keys():                                                            # iterate over all target IDs
        count = counts_sorted.get(hID)
        total = total + int(count)                                                              # calculate the total number of mmseq hits                                                               
        tax_count += 1                                                                          # number of unqiue target IDs
        contig_length = lengths.get(hID)                                                        # get the concatenated contig length 
        taxname = names.get(hID)                                                                # get the scientific name
        av_conlen_perhit = round(contig_length/count, 0)                                        # calculate the average contig length per hit
        out.write(f'{hID}\t{count}\t{taxname}\t{contig_length}\t{av_conlen_perhit}\n')          # write to tsv
    print(f'In {total} hits found {tax_count} unique IDs')                                      # write total and uniquely encountered target numbers
        

