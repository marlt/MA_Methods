#!/home/ha74mit/bin/miniconda3/envs/anc_virus/bin/python3

# translate TaxIDs into their corresponding phylogenetic path using the NCBI taxonomy files "names.dmp" and "nodes.dmp"
# write one file each containing the corresponding read header ID, the TaxID and the phylogenetic path

import sys
import os
import pandas as pd
from collections import defaultdict



file = str(sys.argv[1]) # input: csv table containing classification results of CLARK classify_metagenome.sh

#directory for all output is the inputdirectory you define via the input file
path = os.path.dirname(file)			# save input file path
name = os.path.basename(file)			# save filename from that read names will be extracted
basename = name[:-20]					# cut '_trimmed.results.csv' from input filenames

# format: {readID:[TaxID, phyl_path]}
viral_reads = {}				# dict for viral reads
bacterial_reads = {}			# dict for bacterial reads
unclass_reads = {}				# dict for NA reads






def get_dicts(data_dir="/data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/db/taxonomy/"):
	"""get ncbi data and parse to dicts
	parent_dict: dict returning parent of input id
	scientific_names: returning clade name of input id
	phylo_names: returning taxonomical rank of id
	"""
	col_delimiter = '\t|\t'
	row_delimiter = '\t|\n'
	parent_dict = {}
	phylo_names = {}
	scientific_names = {}
	common_names = defaultdict(set)
	genbank_common_name = defaultdict(set)
	with open(os.path.join(data_dir, 'names.dmp')) as names_file:			# open names.dmp 
		for line in names_file:
			line = line.rstrip(row_delimiter)
			values = line.split(col_delimiter)
			tax_id, name_txt, _, name_type = values[:4]						# save taxID, taxon name and name type
			# print(tax_id)
			# print(name_type)
			# print(name_type=='scientific name')
			if name_type=='scientific name':								# handle different types of names assigned in the NCBI database
				# print("add to scientific names")
				scientific_names[int(tax_id)] = name_txt
			elif name_type == 'common name':
				common_names[int(tax_id)].add(name_txt)
			elif name_type == 'genbank common name':
				genbank_common_name[int(tax_id)].add(name_txt)
	with open(os.path.join(data_dir, 'nodes.dmp')) as nodes_file:			# load the nodes.dmp
		for line in nodes_file:
			line = line.rstrip(row_delimiter)
			values = line.split(col_delimiter)
			tax_id, parent_id, phylo_rank = values[:3]
			tax_id = int(tax_id)
			parent_id = int(parent_id)
			parent_dict.update({tax_id: parent_id})                 # save taxID and its parent ID in the phylogenetic tree, one level up
			phylo_names.update({tax_id: str(phylo_rank)})			# save the rank of each taxID
	return parent_dict, scientific_names, common_names, phylo_names, genbank_common_name


def get_tax_path(id_val, parent_dict, phylo_names):		
	"""generate alist of taxIDs that represent one phylogenetic path"""
	current_id = int(id_val)
	taxonomy = [id_val]			# initialize the list with the "leaf" taxID
	# print(current_id)
	# print(taxonomy)
	while current_id != 1 and phylo_names[current_id] != "superkingdom":			# add parent IDs as long as superkingdom is reached
		current_id = parent_dict[current_id]
		taxonomy.append(current_id)
	# print(taxonomy)
	return taxonomy

# fill the dictionaries
parent_dict, scientific_names, common_names, phylo_names, genbank_common_name = get_dicts(data_dir="/data/prostlocal/programs/clark_1.2.6.1/CLARKSCV1.2.6.1/db/taxonomy/")



# Open clark classification csv, get the first hit TaxID + scientific name
with open(file, "r") as infile:
	count = 0
	for line in infile:				# iterate over all entries in the file
		count+=1
		# if count < 10:
		line_ = line.split(",")
		# print(line_)
		readID = line_[0]			# save the read header ID
		taxID = line_[3]			# track the first assignemtn, second is discarded
		if taxID == 'NA':			# save the NA reads in dictionary
			unclass_reads.update({readID:taxID})
		# print(readID)
		# print(taxID)					
		taxID = int(taxID)
		taxonomy = get_tax_path(taxID, parent_dict, phylo_names)				# gets the taxID path 
		# print(taxonomy)
		pfad = ""
		for number in taxonomy:				
		  	# print(scientific_names[number])
			pfad =  pfad + " > " + str(scientific_names[number])				# converts taxID into scientific name
			# print(pfad)
		pfad = pfad.lstrip(" > ")
		# distinguish between viral and bacteria assignments:
		if "Bacteria" in pfad:										
			bacterial_reads.update({readID:[taxID,pfad]})
		elif "Viruses" in pfad:
			viral_reads.update({readID:[taxID,pfad]})
		# taxID_names[taxID] = pfad



# # output table:
# #	col1: TaxID
# #	col2: Scientific name (taxonpath, NCBI)
# #	>col3: read counts among all files
# output_name = "taxonnames_allcounts.csv"


# write NA, viral and bacterial files 
print("Write viral reads")
with open(f'{path}/{basename}_viral_hit1.csv', 'w') as vfile:
for key in viral_reads.keys():
		vfile.write(f'{key},{viral_reads[key][0]},{viral_reads[key][1]}\n')

print("Write bacterial reads")
with open(f'{path}/{basename}_bacterial_hit1.csv', 'w') as bfile:
	for key in bacterial_reads.keys():
		bfile.write(f'{key},{bacterial_reads[key][0]},{bacterial_reads[key][1]}\n')	

print("Write NA reads")
with open(f'{path}/{basename}_NA_hit1.csv', 'w') as nfile:
	for key in unclass_reads.keys():
		nfile.write(f'{key},{unclass_reads[key]}\n')	


# print(f'Unclassified: \t {len(unclass_reads.keys())}')
# print(f'Bacteria: \t {len(bacterial_reads.keys())}')
# print(f'Virus: \t {len(viral_reads.keys())}')





