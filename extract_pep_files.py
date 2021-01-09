# first argument: peptide.fasta
# second argument: blast output associated with the peptide.fasta
# execute this script in the directory of the blast output

from Bio import SeqIO
import sys
import os
import errno


with open(sys.argv[2]) as f:
    names = []
    names_ext = []
    for x in f.readlines():
        col = x.split('\t')
        names.append(col[1])
        names_ext.append(col[0].split('_')[0]+'_'+col[1])

records = SeqIO.parse(sys.argv[1], "fasta")

plastid_filename = os.getcwd()+ '/plastid_genes/plastid_genes_queries.fasta'
if not os.path.exists(os.path.dirname(plastid_filename)):
    try:
        os.makedirs(os.path.dirname(plastid_filename))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

with open(plastid_filename, 'a') as out:
    for record in records:
        for i, name in enumerate(names):
            if record.id == name:
                SeqIO.write(record, out, 'fasta')
                SeqIO.write(record, os.getcwd() + '/plastid_genes/' + names_ext[i] + ".fasta", 'fasta')
