# rename_fasta_for_phylogenomics
#
# Renames the output of Trinity and TransDecoder v5.1.0 (peptide prediction) to work with Scafos ex. >Genus_species@ID
#
# Usage: python rename_fasta_for_phylogenomics.py [fasta_file_to_rename]
#
# You can also specify many files using '*'
#
# Ex: python rename_fasta_for_phylogenomics.py gregarine_transdecoder_pep.fa 
# gregarine_transdecoder_pep.fa
# Enter species name: Monocystis_agilis
#

import sys
from glob import glob

fasta_files = sys.argv[1]

for fname in glob(fasta_files):
    print fname
    species_name = raw_input('Enter species name: ')
    lines = open(fname, 'r').readlines()
    output = open(species_name+'.fasta.pep', 'w')
    for line in lines:
       if line.startswith('>') == True:
            output.write('>'+species_name+'@'+line.split(' ')[0].strip('>TRINITY_')+'\n')
       else:
           output.write(line)
    output.close()

