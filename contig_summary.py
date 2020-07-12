#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6

# Script to summarize per contig cluster affiliation, lengths, coverage and best mmseq hit (from bit score)

############## IMPORTS ####################
import argparse
import sys, os
from collections import OrderedDict

############# CLI ##################
parser = argparse.ArgumentParser()

# Coverage length information
parser.add_argument('-cl', '--covlen', dest='covlen', metavar='INT', required = True, help="Enter a tab separated file that includes 'contigID \t length \t coverage'.")
# mmseqs results
parser.add_argument('-m', '--mmseqs', dest='mmseqs', metavar='STR', required = True, help="Hand over a mmseqs2 search results file in tsv format. (can be generated with 'mmseqs convertalis')")
# cluster information, header file from result2flat header file
parser.add_argument('-ci', '--cluster', dest='cluster', metavar='STR', required = True, help="Hand over an header file, containing the representative and belonging contig header")
# cluster summary output:
parser.add_argument('-o', '--osum', dest='osum', metavar='STR', required = True, help="Output file for the cluster summary file without sequence information")

args = parser.parse_args()

covlen = args.covlen
mmseqs = args.mmseqs
cluster = args.cluster
osum = args.osum

############## MAIN ########################
# GLOBALS
covlens = {}            # contains all contigs derived from binned reads (see k-mer classification)
mmseq_results = {}      # contains all unique mmseqs search hits (filtered from max 500 hits per query)
headline = 'representativeID/contigID\tcoverage\tlength'
clusternames = {}           # contig IDs saved with representative as key
clustersizes = {}    # representatives with the clustersizes, ordered dict to sort descending

with open(covlen, 'r') as infile:
    for line in infile:
        line_ = line.rstrip('\n')
        entry = line_.split('\t')                        # split line content upon tab delimiter
        header = entry[0].lstrip('>')
        info = ''
        for elem in entry[1:]:
            info = info + '\t' + elem
        covlens.update({header:info})                   # add the coverage and length info as entry with contigID as key

with open(mmseqs, 'r') as infile:
    for line in infile:
        line_ = line.rstrip('\n')
        entry = line_.split('\t')                        # split line content upon tab delimiter
        if entry[0] == "query":                          # this is the header line of the mmseqs results
            header = ''
            for elem in entry[1:]:
                header = header + '\t' + elem
            header = header + '\n'
            headline = headline + header                # append header line about the mmseqs column header
        else:
            hit = ''
            if entry[1] == "NZ_KI911782.1":
                print(entry)
            for elem in entry:
                if elem == entry[0]:
                    contigID = elem
                else:    
                    hit = hit + '\t' + elem
            mmseq_results.update({contigID:hit})           # add the mmseqs entry with contigID as key


before = '' 
represent = ''          # 
contigs = []
size = 0
with open(cluster, 'r') as infile:
    for line in infile:
        line1 = line.rstrip('\n')
        line_ = line1.lstrip('>')
        # print(line_)
        # print(size)
        if line_ == before:                                     # new cluster encountered
            if represent != '':
                real_size = size -1
                contigs.remove(before)                              # the latest element was the new representative of the new cluster, so remove from this!
                clustersizes.update({represent:real_size})           # write the cluster representative and the cluster size
                clusternames.update({represent:contigs})            # write all
                # print(f'{represent} {real_size}')
            represent = line_                                      # this is the new representative now
            contigs = []
            contigs.append(represent)
            size = 1
        else:
            size += 1
            before = line_
            contigs.append(before)

# print(covlens)
# print("\n")
# print(mmseq_results)
# print("\n")
# print(clusternames)
# print("\n")
# print(clustersizes)

clustersizes_sorted = OrderedDict(sorted(clustersizes.items(), key= lambda x: x[1], reverse = True))
# print(clustersizes_sorted)

# with open(osum, 'w') as outfile:
#     outfile.write(f'{headline}')
#     for contigID in clustersizes_sorted.keys():
#         outfile.write(f'>{contigID}\t{clustersizes[contigID]}\n')           # representative is written with leading '>' and with clustersize
#     # for contigID in clusternames.keys():                    
#         for contig in clusternames[contigID]:                               # iterate over all cluster members, including the representative
#             try:
#                 mmseqs_text = mmseq_results.get(contig)                     # try to get mmseqs hit information
#             except KeyError:
#                 mmseqs_text = 'NoHit'     
#             # print(contig)                                  # of no homologies were written, write 'not assigned'
#             # print(covlens.get(contig))
#             covlen_text = covlens.get(contig)[1:]                               # get coverage and length information
#             covlen_text = covlen_text.rstrip('\t')
#             # print(covlen_text)
#             # if contig == contigID:
#             #     outfile.write(f'>{contig}\t{clustersizes[contigID]}\t{covlen_text}\t{mmseqs_text}\n')
#             # else:
#             if mmseqs_text != None:
#                 mmseqs_text = mmseqs_text.lstrip('\t')    
#             outfile.write(f'{contig}\t{covlen_text}\t{mmseqs_text}\n')



