#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6

# Calculate the clustersize from a mmseqs result2flat header file

############ IMPORTS #######################
import sys, os

############### MAIN ########################
# GLOBALS
clusters = {}       # save clustersize as value of the cluster's representative

infile = sys.argv[1]        # input file: headers from results2flat output: can be generated with 'awk '/>/ {print $0}' <file.fasta>
outfile = sys.argv[2]       # file to save: tab-separated representatives and clustersize

before = ''                 # the header of the line before
represent = ''              # the representative of the cluster
size = 0                    # size of the cluster

with open(infile, 'r') as fa:
    for line in fa:
        if line == before:          # new cluster encountered
            if represent != '':
                real_size = size -1
                represent = represent.rstrip("\n")
                clusters.update({represent:real_size})   # write the cluster representative and the cluster size
                print(f'{represent} {real_size}')
            represent = line
            size = 1
        else:
            size += 1
            before = line
    # writing last element
    represent = represent.rstrip("\n")
    clusters.update({represent:size})   # write the cluster representative and the cluster size
    print(f'{represent} {size}')

with open(outfile, 'w') as out:
    for key in clusters.keys():
        out.write(f'{key}\t{clusters[key]}\n')
