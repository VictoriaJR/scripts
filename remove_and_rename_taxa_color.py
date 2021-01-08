import os
import sys
import glob

def read_fasta(filename):
  infile = open(filename)
  line = infile.read()
  infile.close()
  seqs = line[1:].split('\n>')

  seq_d = {}
  for seq in seqs:
    seq_d[seq.split('\n')[0]] = ''.join(seq.split('\n')[1:])

  return seq_d


try:
  treefile = sys.argv[1]
  color_delete = sys.argv[2]
  color_rename = sys.argv[5]
  aln = sys.argv[3]
  outfile = sys.argv[4]



except IndexError:
  print "missing input parameter"

out = open(outfile,'w')


infile = open(treefile)
line = infile.read()
infile.close()

list = line.split('taxlabels\r\n')[1]
list = list.split(';\nend;\n')[0]

list = list.split('\r\n')

taxa = []
for i in list:
  if i.strip() != '':
    taxa.append(i.strip())

#print taxa

new_names = {}

for i in taxa:
#  print i.split('=#')[-1][:-1]
#  print color

  if i.split('=#')[-1][:-1] == color_rename:
    name = i.split('[')[0]
    name = name.replace("'",'')
    print name
    new_name = '%s-1@%s' % (name.split('@')[0], name.split('@')[1])
    new_names[name] = new_name

  elif i.split('=#')[-1][:-1] != color_delete:
    name = i.split('[')[0]
    name = name.replace("'",'')
    new_names[name] = name
  else:
    print '%s fucking deleted' % (i)
   
seq_d = read_fasta(aln)

out = open(sys.argv[4],'w')

for name in new_names:
  out.write('>%s\n%s\n' % (new_names[name], seq_d[name]))

out.close()

