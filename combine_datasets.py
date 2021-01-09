# COMBINE_DATASETS.PY
#
# Given a list of folders, combine the datasets contained in each folder without duplicating sequences.
#
# Input: two or more folders containing gene files. The folders must all be within one folder whose path must be given in the variable inPATH below.
#
# Output: one out file for each gene file; each out file contains all the datasets (from the input folders) with no duplicates. The out files are created in the folder whose path must be given in the variable outPATH below.
#
# NOTE: the code assumes that for a given gene, two repeated organism's sequence ID (or sequence name) have the exact same protein sequence.
#
# EXAMPLE: python combine_datasets.py dataset1/ dataset2/



import sys
import glob
import os, os.path




inPATH = os.getcwd() #sets the output file path to your current working directory
outPATH = inPATH + "out_folder/"
if not os.path.exists(os.path.dirname(out_folder)):
    try:
        os.makedirs(os.path.dirname(out_folder))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

table_protein_seq = {} # creates a dictionnary according to the relation: file ↦ list of protein sequence

for i in range(1,len(sys.argv)): # loop through the folders given as arguments
    DIR = inPATH+str(sys.argv[i])
    for file in os.listdir(DIR): # loop through each file within a folder
        open(outPATH+str(file)+'.out', 'a').close() # creates a new file if it does not already exist
        table_protein_seq.append(file : [''])


for i in range(1,len(sys.argv)): # loop through the folders given as arguments
    DIR = inPATH+str(sys.argv[i])
    for file in os.listdir(DIR): # loop through each file within a folder
        with open(outPATH+str(file)+'.out', 'a') as out: # append sequence name and sequence protein in the corresponding out file
            sequences = open(DIR+'/'+str(file), 'r').read()
            lines = sequences.split('>')
            n = 1
            while n < len(lines):
                parts = lines[n].split('\n')
                seq_name = parts[0].split(' ')[0]
                protein_seq = ''.join(parts[1:-1])

                if protein_seq not in table_protein_seq[file]:
                    out.write('>'+seq_name+'\n') # write sequence name
                    out.write(protein_seq+'\n') # write protein sequence
                    table_protein_seq[file].append(protein_seq) # add the new protein sequence to the list corresponding to the in-file

                n += 1

            # out.close() # should be superfluous -> with already close it at the end
