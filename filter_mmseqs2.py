#!/usr/bin/python3

# Filter the highest bit scoring hit from potentially all 500 hits 


##### IMPORTS #########
import sys, os



############ FUNCTIONS ############################

def sort_cluster(cluster):              # sort the hit cluster descending sorted on score
    print("Sort list...")
    cluster.sort(key= lambda x: x[4], reverse = True)  # sort by score descending
    print("Done.")

def append_out(outfile, cluster):                # write highest score entry of the cluster to file
    with open(out, 'a') as outf:
        string = ""
        # print(hits_cluster[0])
        for elem in hits_cluster[0]:
            string = string + str(elem) + delim
        string = string.rstrip(delim)
        outf.write(string)       # write hit with highest score to output csv

######## MAIN #########################

csv=str(sys.argv[1])    # hand over a blast output csv file
delim="\t"              # specify the set delimiter in csv file
out=str(sys.argv[2])    # outfile
                            

clusternum=1
with open(out, 'w') as outf:
    # retain the blast-like tsv format
    outf.write("query" +delim+	"target"+delim+	"taxid"+delim+	"taxname"+delim+	"bits"+delim+"alnlen"+delim+	"pident"+delim+	"evalue"+delim+	"qlen"	+delim+"qstart"+delim+	"qend"+delim+	"qcov"+delim+	"tlen"	+delim+"tstart"+delim+	"tend"+delim+	"tcov"+"\n")

with open(csv, 'r') as file:
    line1 = "$" # initialize line1
    queryA = "$"
    hits_cluster = []
    while line1 != "":          # as long as no EOF
        line1 = file.readline()
        if line1 == "":            # if readline encounters end of file, return is empty string
            print("End of file")
            break
        elif line1 != "":
            line1_ = line1.split(delim)
            line1_[4] = int(line1_[4])
            queryB = line1_[0]
            if queryA == queryB or queryA == "$":       # if first entry or hit on same query as before
                hits_cluster.append(line1_)
                queryA = queryB
            elif queryA != queryB and queryA != "$":    # if hit of new query occurs
                clusternum+=1
                sort_cluster(hits_cluster)              # sort the cluster       
                append_out(out, hits_cluster)           # write the highest score hit to output file
                queryA = queryB                         # hand over the new query name for comparing the next line
                hits_cluster = []                       # remove the previous cluster entries
                hits_cluster.append(line1_)             # put in the first element of the new hit cluster
                print(clusternum)
            else:
                print("Unclassified hit. Please clearify.")





