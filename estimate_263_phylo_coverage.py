# This script writes a file containing the sequences matching the first input keyword found from each file contained in the input directory. To call the script, the user must provide:
# 1. a first argument which is the (full) path of the directory containing the fasta files. Note: must write the '/' symbol at the end
# 2. a list of keywords to be found in these files. Note: the list must be ordered as if a match is found, the next keywords will not be searched
# Example: python phylo_protein_extract.py PATH_OF_MY_DIRECTORY Dinophyceae Gyrodinium Psammosa

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import glob
import os, os.path

DIR = str(sys.argv[1]) # directory name with / at the end
keywords = sys.argv[2:] # search for these keywords

new_filename = DIR + 'lineage_phylogenomic_one_query_per_263_protein.fasta'

files_without_keywords = []
proteins = []

for file_ in os.listdir(DIR): # loop through each file within a folder
    if file_ != 'lineage_phylogenomic_one_query_per_263_protein.fasta':
        kw_match = False # kw_match assets if a keyword has been found in the records
        # print "===== Open the file " + file_ + " with keyword match set to " + str(kw_match)
        for kw in keywords:
            # print "          Start the search for the keyword " + kw + " with keyword match set to " + str(kw_match) + " in " + file_
            for record in SeqIO.parse(DIR + file_, "fasta"): # a keyword corresponds to a protein
                if kw in record.id: # if one of the given keywords is in the name of a sequence
                    # print "          Match obtained (keyword match was set to " + str(kw_match) +")"
                    kw_match = True
                    new_record = SeqRecord(record.seq, file_ + "_" + record.id, "", "")
                    proteins.append(new_record)
                    break
            if kw_match:
                # print "          Match successful in " + file_ + " with keyword " + kw
                break
            # else:
                # print "          ~~~~~ FAILURE TO FIND " + kw + " in " + file_
        if not kw_match:
            files_without_keywords.append(file_)

print "The following files do not contain any of the input keywords:"
print files_without_keywords

SeqIO.write(proteins, new_filename, "fasta") # creates a new file if it does not already exist
