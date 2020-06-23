#!/usr/bin/python3
###07.02.20

# script visualizes all adapter counts per fastq files in a heatmap

######################### IMPORTS #################################
import matplotlib.pyplot as plt
import numpy as np
import sys, os

####################### FUNCTIONS #################################
def load_data(table):
    '''Load tab separated output table generated with rg_adapters.sh'''
    file_labels = []                                                # save fastq file name here
    adap_labels = []                                                # save the adapter identifiers here here
    adap_counts = []                                                # 2D array including adapter counts per fastq (line: fastq file; column: adapter content)
    
    with open(table, 'r') as tsv:                                   
        for line in tsv:                                            # iterate over each line in tsv file
            norm_counts = []
            line_ = line.replace(" ", "\t")                         # rg_adapters did not correctly placed delimiters: replace with tab
            line_ = line_.rstrip("\n")                              # get rid of terminal newline
            entry = line_.split('\t')                               # split line on tab
            if entry[0] == "Filename":                              # headline of file: 
                adap_labels.extend(entry[2:])                       # extract the adapter names from here
            elif "trimmed" not in entry[0] and "allSE" in entry[0]: # exclude the merged (unflag+single) files from before filtering, they are displayed separately; this does not cover the trimmed allSE from trimmed fastq files!
                pass
            elif "NIOBE" in entry[0]:                               # exclude the ds library files, all starting with "NIOBE" in file
                pass
            else:
                label = trim_filename(entry[0])                     # trim teh filename to unqiue but shorter names for heatmap labeling
                file_labels.append(label)
                readnum = entry[1]                                  # store read number of a fastq file
                for count in entry[2:]:                             # for all adapter counts in fastq file
                    norm_count = int(count)/int(readnum)*1000       # calculate the relative adapter count per 1000 reads
                    norm_counts.append(norm_count)                  # collect them in norm_counts array
                adap_counts.extend([norm_counts])                   # finally write the adapter counts to next line of the 2D adapter counts array for visualization
    adap_counts = np.array(adap_counts)                             # convert to numpy array

    # print(adap_counts)
    # print(file_labels)
    # print(adap_labels)

    return adap_counts, file_labels, adap_labels


def trim_filename(fname):
    ''' Modify filenames to unique and meaningful but shorter versions'''
    if "trimmed" in fname:                                          
        label = fname[:-8]                                         
    else:
        label = fname
    # print(label)
    return label

def heatmap(counts, xlabel, ylabel, outpath):
    ''' Calculate the heatmap for a 2D numpy array with labeling data and save the figure'''
    # split counts and fastq file labels into half
    counts1 = counts[:int(len(counts)/2)]                           
    counts2 = counts[int(len(counts)/2):]

    ylabel1 = ylabel[:int(len(ylabel)/2)]
    ylabel2 = ylabel[int(len(ylabel)/2):]

    # plot in subplots, one row, two columns
    fig, (ax1, ax2) = plt.subplots(1,2)

    # plot the generated halfs separately
    hm1 = ax1.imshow(counts1, aspect = 'auto') #, aspect = 'auto'
    hm2 = ax2.imshow(counts2, aspect = 'auto') #, aspect = 'auto'

    # set x-ticks (for orientation of adapter labels)
    ax1.set_xticks(np.arange(len(xlabel)))
    ax2.set_xticks(np.arange(len(xlabel)))

    # set y-ticks (for orientation of adapter labels)
    ax1.set_yticks(np.arange(len(ylabel1)))
    ax2.set_yticks(np.arange(len(ylabel2)))

    # write adapter identifiers as x labels
    ax1.set_xticklabels(xlabel, fontsize = 6)
    ax2.set_xticklabels(xlabel, fontsize = 6)

    # write fastq names as y labels (in halfs)
    ax1.set_yticklabels(ylabel1, fontsize = 5.5)
    ax2.set_yticklabels(ylabel2, fontsize = 5.5)

    # remove spines between individual pixels ("heatmap boxes")
    for edge, spine in ax1.spines.items():
        spine.set_visible(False)
    for edge, spine in ax2.spines.items():
        spine.set_visible(False)

  
    # write a colorbar for assignment of color to exact adapter counts
    cbar = ax2.figure.colorbar(hm2, ax=ax2)
    cbar.ax.set_ylabel("Rel. Matches per 1000 Reads", rotation=-90, va="bottom", fontsize = 10)
    cbar.ax.tick_params(labelsize = 6)
    
    # rotate x-labels by 45°
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # add supertitle
    fig.suptitle("Exact Adapter Matches on Trimmed Reads")
    fig.tight_layout()
    # save figure to outpath
    plt.savefig(outpath)
    plt.show()

####################### MAIN #####################################

# Input file:
infile=str(sys.argv[1])             # specify path to tab-separated table containing the adapter counts (generated with "rg_adapters.sh")

counts, flabel, alabel = load_data(infile)

heatmap(counts, alabel, flabel, '/data/mahlzeitlocal/projects/ma_neander_assembly/anc_virus_MA/visualization/preprocessing/s.pdf')

