# single_line_fasta.py
#
# Convert a multiline fasta into a single line fasta.
# Outputs *.fasta.sl
#
# Usage: python single_line_fasta.py [multiline.fasta]
# Example: python single_line_fasta.py multiline.fasta
#

import sys

n = 0

file = open(sys.argv[1], 'r')
lines = file.readlines()
file.close()
output = open(sys.argv[1] + '.sl', 'w')

for i in lines:
	if i.startswith('>') == True:
		if n == 0:
			output.write(i)
			n += 1
		else:
			output.write('\n' + i)
			n += 1
	else:
		output.write(i.strip('\n'))
		n += 1
output.write('\n')
output.close()
