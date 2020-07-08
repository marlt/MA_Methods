#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3

# DESCRIPTION
# Filter singletons from paired end fastq files
 
# PACKAGES
import argparse
import sys, os 

# CLI

parser = argparse.ArgumentParser()

# pair1
parser.add_argument('-i1', '--input1', dest='input1', action='store', nargs=1, metavar='PATH', required=True, help='Input file 1')
# pair2
parser.add_argument('-i2', '--input2', dest='input2', action='store', nargs=1, metavar='PATH', required=True, help='Input file 2')
# out pair1
parser.add_argument('-o1', '--output1', nargs=1, action='store', metavar='PATH', required=True, help='Output file 1')
# out pair2
parser.add_argument('-o2', '--output2', nargs=1, action='store', metavar='PATH', required=True, help='Output file 2')
# singleton
parser.add_argument('-s', '--single', nargs=1, action='store', metavar='PATH', required=True, help='Singleton files')

args = parser.parse_args()

# GLOBALS

input1 = args.input1[0]
input2 = args.input2[0]
output1 = args.output1[0]
output2 = args.output2[0]
single = args.single[0]


# FUNCTIONS

# load fastq file
def load_fastq(infile):
    '''Load the paired fastq files'''
    with open(infile, 'r') as file:
        count = 0
        container = {}
        for line in file:               # read the header, every fourth line
            count += 1
            header = line.rstrip('\n')
            seq = file.readline().rstrip('\n')      # get sequence
            plus = file.readline().rstrip('\n')     # get additional options line
            phred = file.readline().rstrip('\n')    # get the phred score line
            container.update({header:[seq,plus,phred]})
        assert count == len(container.keys())
    return container


# MAIN
# load pair1
pair1_fastq = load_fastq(input1)                # sequences are saved with header as key in dictionary
# load pair2
pair2_fastq = load_fastq(input2)

# get number of entries in pair1 fastq
length1 = len(pair1_fastq.keys())
# get number of entries in pair2 fastq
length2 = len(pair2_fastq.keys())

# print(length1)
# print(length2)

count_pair = 0
count_single1 = 0
count_single2 = 0

# open three outout channels for the pair output files and the singleton
with open(output1, 'w') as out1:
    with open(output2, 'w') as out2:
        with open(single, 'w') as singleton:
            for header in pair1_fastq.keys():       # conduct all header IDs in pair1
                con1 = pair1_fastq.get(header)      # the h1 header must be in the dict 
                try:
                    con2 = pair2_fastq.get(header)   # if the header is found, write the entries in the pair1 and pair2 file
                    pair2_fastq.pop(header)          # remove the paired read from pair2 dict
                    out1.write(header + '\n')        # write read1 to outfile1
                    out2.write(header + '\n')        # write read2 to outfile2
                    for line in con1:
                        out1.write(line + '\n')
                    for line in con2:
                        out2.write(line + '\n')
                    count_pair += 1
                except KeyError:                # if the header is not in the pair2 file, it is a singleton and written to the single file
                    singleton.write(header + '\n')  
                    count_single1 += 1
                    for line in con1:
                        singleton.write(line + '\n')
            for header in pair2_fastq.keys():   # the remaining reads in pair2 are all singletons since the pairs had been removed
                con = pair2_fastq.get(header)
                singleton.write(header + '\n')  # write them to the singleton file
                count_single2 += 1
                for line in con:
                    singleton.write(line + '\n')

print(f'paired reads:\t {count_pair}')
print(f'singletons 1: \t {count_single1}')
print(f'singletons 2: \t {count_single2}')
print(f'single TOTAL: \t {count_single1 + count_single2}')
print(f'files in 1:\t {count_pair + count_single1} == {length1}')
print(f'files in 2:\t {count_pair + count_single2} == {length2}')

