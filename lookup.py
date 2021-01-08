# LOOKUP.PY
#
# Remove the sequences defined in a list from a fasta file.
#

import sys

n = 1
sequences = open(sys.argv[1], 'r').read()
lines = sequences.split('>')
out = open(sys.argv[1] + '.clean', 'w')
remove = open(sys.argv[2],'r').readlines()

while n < len(lines):
    parts = lines[n].split('\n')
    seq_name = parts[0].split(' ')[0]+'\n'
    if seq_name not in remove:
        out.write('>'+seq_name)
        for part in parts[1:-1]:
            out.write(part)
        out.write('\n') 
    n = n + 1
out.close()
