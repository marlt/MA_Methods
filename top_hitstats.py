#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6

# append n target ID clusters with the corresponding average homology search measures

############ IMPORTS ###############################

import argparse, sys, os
import numpy as numpy
from collections import OrderedDict 

################ CLI #######################
parser = argparse.ArgumentParser()

# target ID cluster tsv
parser.add_argument('-i', '--inIDs', dest='inIDs', metavar='STR', required = True, help="Sorted file comprising cumulative contig length on hit ID cluster.")
# combined contig cluster and unique MMseqs2 results tsv
parser.add_argument('-m', '--mmseqs', dest='mmseqs', metavar='STR', required = True, help="MMSEQs unique hit file")
# number of target ID clusters to process
parser.add_argument('-n', '--numhits', dest='numhits', metavar='INT', required = True, default = 50, help="Number of target cluster IDs to process")
# tab separated output statistics
parser.add_argument('-o', '--output', dest='outfile', metavar='STR', required = True, help="Output is statistics for all hits on alnlen, pident, bitscore, average query coverage")

args = parser.parse_args()

# parse the CLI parameters
inIDs = args.inIDs
mmseqs = args.mmseqs
if type(args.numhits) == list:
    num = args.numhits[0]
else:
    num = args.numhits
num = int(num)   
outfile = args.outfile


################# FUNCTIONS #######################
def load_IDs(id_file, hitnum):                                          # load the target cluster IDs
    id_dict = OrderedDict()
    with open(id_file, 'r') as itfile:
        linecount = 0
        for line in itfile:
            linecount += 1
            if linecount <= hitnum:                                     # load as many as restricted by hitnum
                entry = line.split("\t")
                # print(entry)
                id_info = entry[1:]                                     # hitcounts, taxname, cum_contiglen, av_contiglen
                id_dict.update({entry[0]:id_info})                      # save entry
            else:
                # print(f'Extract from {linecount-1} IDs...')
                break   
    return id_dict

def get_dicts(id_dict, mm_cluster_list):                                # collect the homology search measures per target cluster ID (used as key); one dict used per measure
    pident = OrderedDict()                                              # percentage identity
    alnlen = OrderedDict()                                              # alignment length
    bitscore = OrderedDict()                                            # bit score
    qlen = OrderedDict()                                                # query length
    depth = OrderedDict()                                               # average read depths on contig
    cluster = OrderedDict()                                             # number of contig clusters covered (contigs in one cluster share pairwise seq identity of 80%)                            

    id_list = id_dict.keys()

    with open(mm_cluster_list, 'r') as mmf:
        count = 0
        for line in mmf:
            entry = line.split('\t')
            if entry[0].startswith('>'):                                # if new cluster is encountered:
                rep = entry[0]                                          # save as representative
            else:
                if entry[3] in id_list:                                 # test, if the encountered contigs has a hit of the loaded target cluster IDs       
                    if pident.get(entry[3]) == None:                    # if yes, but their is not yet an entry for this cluster, initialise all dictionary entries
                        # print(entry[3])
                        pident.update({entry[3]:[float(entry[8])]})     
                        alnlen.update({entry[3]:[int(entry[7])]})
                        bitscore.update({entry[3]:[int(entry[6])]})
                        qlen.update({entry[3]:[int(entry[10])]})
                        depth.update({entry[3]:[float(entry[1])]})
                        cluster.update({entry[3]:[rep]})
                    else:                                               # if encountered again, append the homology measure entries of the target hit ID in the dictionaries
                        pident[entry[3]].append(float(entry[8]))
                        alnlen[entry[3]].append(int(entry[7]))
                        bitscore[entry[3]].append(int(entry[6]))
                        qlen[entry[3]].append(int(entry[10]))
                        depth[entry[3]].append(float(entry[1]))
                        if rep not in cluster[entry[3]]:
                            cluster[entry[3]].append(rep)
    return pident, alnlen, bitscore, qlen, depth, cluster

def av_compute_avrg(id_dict, p, a, b, q, d, cluster):                  # compute the (mean) values for all selected target cluster IDs

    av_pident = OrderedDict()           # average percentage identity
    av_alnlen = OrderedDict()           # average alignment length
    av_bitscore = OrderedDict()         # average bit score
    av_qlen = OrderedDict()             # average contig length target ID cluster
    av_depth = OrderedDict()            # average read depth on contigs of a cluster
    av_coverage = OrderedDict()         # average contig coverage 
    clusters = OrderedDict()            # number of contig clusters covered in the target ID cluster
    cum_qlen = OrderedDict()            # cumulative contig length, as computed before

    for headID in id_dict.keys():           # go through the list of target cluster IDs
        # this two measures can be read from previous calculation
        cum_q = int(id_dict[headID][2])
        av_q =  float(id_dict[headID][3])
        length = len(p[headID])
        # calculate the remaining measures
        av_p = round(sum(p[headID])/length, 3)
        av_a = round(sum(a[headID])/length, 1)
        av_b = round(sum(b[headID])/length, 1)
        av_d = round(sum(d[headID])/length, 2)
        av_cov = round((sum(a[headID])/length)/(cum_q/length), 3)
        unique_clusters = len(cluster[headID])
        # write them to the dictionaries to store the average homology measures per target cluster ID
        av_pident.update({headID:av_p})
        av_alnlen.update({headID:av_a})
        av_bitscore.update({headID:av_b})
        av_qlen .update({headID:av_q})
        av_depth.update({headID:av_d})
        av_coverage.update({headID:av_cov})
        clusters.update({headID:unique_clusters})
        cum_qlen.update({headID:cum_q})
        # assert 
    return av_pident, av_alnlen, av_bitscore, av_qlen, av_depth, av_coverage, clusters, cum_qlen

def write_stats(outfile, id_dict, av_pident, av_alnlen, av_bitscore, av_qlen, av_depth, av_coverage, clusters, cum_qlen): # write the output tsv for the calculated output tsv
    # header to write first to the output file
    header = "Target ID\tTarget Name\tHits\tCumulative Contig Length\tAv. Contig Length\tAv. Depth\tAv. Alignment Length\tAv. Contig Coverage\tAv. Pident\tClusters\tav. Bit Score\n"
    with open(outfile, 'w') as out:
        # print(header)
        out.write(header)
        for headerID in id_dict.keys():
            # print(id_dict.get(headerID)[1])
            # write an target cluster ID per line, including the calculated measures
            out.write(f'{headerID}\t{id_dict.get(headerID)[1]}\t{id_dict.get(headerID)[0]}\t{cum_qlen.get(headerID)}\t{av_qlen.get(headerID)}\t{av_depth.get(headerID)}\t{av_alnlen.get(headerID)}\t{av_cov.get(headerID)}\t{av_pident.get(headerID)}\t{clusters.get(headerID)}\t{av_bitscore.get(headerID)}\n')
            # print(f'{headerID}\t{id_dict.get(headerID)[1]}\t{id_dict.get(headerID)[0]}\t{cum_qlen.get(headerID)}\t{av_qlen.get(headerID)}\t{av_depth.get(headerID)}\t{av_alnlen.get(headerID)}\t{av_coverage.get(headerID)}\t{av_pident.get(headerID)}\t{clusters.get(headerID)}\t{av_bitscore.get(headerID)}')

##################### MAIN #######################

ids = load_IDs(inIDs, num)

pi, al, bs, ql, dep, clust = get_dicts(ids, mmseqs)
av_pi, av_al, av_bs, av_ql, av_dep, av_cov, clusts, cumlen = av_compute_avrg(ids, pi, al, bs, ql, dep, clust)


write_stats(outfile, ids, av_pi, av_al, av_bs, av_ql, av_dep, av_cov, clusts, cumlen)

