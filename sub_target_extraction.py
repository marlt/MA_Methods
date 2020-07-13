#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3.6

# Extract the top 10 accIDs, get the fa_file sequences and write all to a multi-fa_file file
# requires a sorted accession ID and scientific name mapping

import sys, os

sci_keys = str(sys.argv[1])
accIDs = []
print(len(accIDs))

fastas = {}         # save the fasta records of the given target IDs

dbfasta = "/data/fass1/database/blast_bak_virus_fungi_human/GTDB_TMGVR_REFSEQ/all_sequences.fna"   # bulk fasta file (MMseqs2 search database) for candidate selection
# dbfasta = "/data/fass1/database/blast_bak_virus_fungi_human/GTDB_TMGVR_REFSEQ/taxonomy/ssu.fna"  # bulk fasta file (GTDB taxonomy 16S rRNAs) for tree construction (required from cactus)
with open(sci_keys, 'r') as ids:            # write the target IDs from text file to array
    for line in ids:
        line_ = line.rstrip("\n")
        accIDs.append(line_)

print(f'Found {len(accIDs)} accessionIDs:')

count = 0
header = ""
seq = ""

with open(dbfasta, 'r') as fa_file:                 # go through the bulk fasta file
    for line in fa_file:
            count += 1
            if count % 1000000 == 0:
                print(count)
            if line.startswith('>'):                # check header
                if header != "" and seq != "":      # if header and sequence were recorded, save them to the fastas dictionary
                    fastas.update({header:seq})
                    seq = ""
                    header = ""
                for accID in accIDs:                # go through the list of target IDs
                    if accID in line:               # if a target ID is found in the list
                        print("Found")
                        header = line               # set header only if target ID was found in line = one of the desired hits; otherwise stays empty and sequence will not be written
                        accIDs.remove(accID)
            elif not line.startswith('>') and header != "":     # collect fasta sequence for an encountered candidate
                seq = seq + line                    
         

print(len(fastas.keys()))

for headid in fastas.keys():            # go through the target IDs of candidates
    accID = headid.split(' ')
    accID_ = accID[0].lstrip('>')
    # write each to an individual fasta file, name set to the target ID
    with open(f'/data/fass2/projects/ma_neander_assembly/MA_data/consensus/microbacterium/genomes/{accID_}.fasta', 'w') as out:
        out.write(headid)               # write header
        out.write(fastas[headid])       # write sequence
        print(f'{headid} length:\n{len(fastas[headid])}')

