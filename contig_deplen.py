#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6


# The script merges contig lengths and depths on equal contig identifiers
# IMPORTS
import matplotlib.pyplot as plt
import matplotlib.axes as axes
import numpy as np
import sys, os



# Specify input files
len_file = sys.argv[1]
depth_file = sys.argv[2]
outpath = sys.argv[3]


depth = {}
length = {}

# average read depth of all contigs is loaded to a dictionary
with open(depth_file, 'r') as depthf:
    for line in depthf:
        elem = line.split('\t')
        if elem[0] != "qseqid":
            depth.update({elem[0]:elem[1].rstrip('\n')})    # 1.col: contig ID; 2.col: average depth of the contig

# length of all contigs is loaded to a dictionary
with open(len_file, 'r') as lenf:
    for line in lenf:
        elem = line.split('\t')
        length.update({elem[0]:elem[1]})        # 1.col: contig ID; 2.col: length of the contig

with open(f'{outpath}/final_contigs_cov_length.csv', 'w') as out:
    for key in depth.keys():
        out.write(f'{key}\t{depth[key]}\t{length[key]}\n')