#!/usr/bin/python3


# python script to compute average read depth on contigs
# Script works for samtools depth output files

# IMPORTS
import sys, os
import pandas as pd

# infile is samtools depth file
infile = str(sys.argv[1])

# get path
path = os.path.dirname(infile) + "/"
# get basename
name = os.path.basename(infile)
bname = name[:-10]



out_depth = pd.DataFrame(columns = ['qseqid','ntcov'])
prev_name = ""
ident = ""
nt_pos = 0
nt_sum = 0
pos = 0
counter = 0
depth = 0

zaehle = 0
with open(infile, 'r') as depthfile:
    for line in depthfile:
        # if (zaehle <= 800): 
        line_content = line.split('\t')
        ident = str(line_content[0])            # the first column is the contig name
        pos = int(line_content[1])              # second column is position in contig
        # print("pos: \t " + str(pos))
        nt_pos = int(line_content[2])           # third column is depth for base
        zaehle += 1
        # print(line + '\t' + str(zaehle) + '\t' )
        
        
        if (prev_name == ""):   # initialise: for line 1
            prev_name = ident
            counter = pos
            nt_sum = int(nt_pos) 
        elif (prev_name == ident):              # if the contig name stays the same
            nt_sum+=nt_pos                      # add up the depth
            counter = pos                       
        elif (prev_name != ident and prev_name != "" or (zaehle+1) == len(depthfile)):  # if new contig loaded, meaning the ID in the first column changes
            depth = nt_sum/counter                                                      # contig is complete, caclulate the depth
            data = {"qseqid":[prev_name], "ntcov":[depth]}                          # save the contig ID and the average read depth
            data_pd = pd.DataFrame(data, columns = ['qseqid','ntcov'])              # convert into pandas dataframe                   
            out_depth = out_depth.append(data_pd, ignore_index = True)              # append the summary output dataframe by the new contig information
            counter = pos                                                           # set position counter to the first position in the new contig
            nt_sum = int(nt_pos)                                                    # set the initial depth
            prev_name = str(ident)                                                  # set the new contig name
            depth = 0                                                               # reset average depth


# Adding the last element
depth = nt_sum/counter
data = {"qseqid":[prev_name], "ntcov":[depth]}
data_pd = pd.DataFrame(data, columns = ['qseqid','ntcov'])                        
out_depth = out_depth.append(data_pd, ignore_index = True)



counter = pos
nt_sum = int(nt_pos)
prev_name = str(ident)
depth = 0


print(out_depth)
print(path)
print(bname)

# # Write a tsv to same directory as input: 
out_depth.to_csv(path + bname + "_coverage.tsv", sep="\t")
#   1. column: contig name
#   2. column:  average depth
