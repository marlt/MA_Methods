#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3

# This script filters reads from a fastq file on a GC content threshold defined by the user
 
# PACKAGES
import argparse
import sys, os 

# CLI
parser = argparse.ArgumentParser()

# Define a GC threshold
parser.add_argument('-gc', '--gc-content', dest='gc', action='store', metavar='INT', nargs=1, default=20, required=True, help='Define NUM GC threshold for filtering.')
# Input multi-fastq file
parser.add_argument('-i', '--input', dest='infile', action='store', metavar='STR', nargs=1, required=True, help='Define input fastq.')
# Ooutput path for the filtered reads
parser.add_argument('-o', '--output', dest='outfile', action='store', metavar='STR', nargs=1, required=True, help='Define output fastq.')
# Specify a path for reads that dont meet the threshold
parser.add_argument('-f', '--failed', dest='e_file', action='store', metavar='STR', nargs=1, required=True, help='Save the filtered reads fastq.')

args=parser.parse_args()

# GLOBALS

infile = args.infile[0]
e_file = args.e_file[0]
outfile = args.outfile[0]
print(infile)
print(e_file)
print(outfile)
gc = int(args.gc[0])

# Open an output stream to the qualified and failed output files in parallel
with open(outfile, 'w') as ofile:
    with open(e_file, 'w') as efile:
        with open(infile, 'r') as ifile:
            for line in ifile:                                  # get every fourth line; this is the header
                header = line.rstrip('\n')                      
                seq = ifile.readline().rstrip('\n')                # get the sequence
                plus = ifile.readline().rstrip('\n')                # get the additional comment line
                phred = ifile.readline().rstrip('\n')               # get the phred score
                gc_count = seq.count('G') + seq.count('C')          # count occurences of bases G and C in the read
                part = (gc_count * 100) / int(len(seq))             # calculate the GC fraction in the read
                #print(seq)
                #print(f'gc_count:\t{gc_count}')
                #print(f'length:\t{len(seq)}')
                #print(f'part:\t{part}')
                if part > gc:               # if threshold is passed, write to qualified output
                    ofile.write(f'{header}\n')
                    ofile.write(f'{seq}\n')
                    ofile.write(f'{plus}\n')
                    ofile.write(f'{phred}\n')
                else:                        # if below threshold, write to failed fastq file
                    efile.write(f'{header}\n')
                    efile.write(f'{seq}\n') 
                    efile.write(f'{plus}\n')
                    efile.write(f'{phred}\n')


# Parameters: consider CLI Description