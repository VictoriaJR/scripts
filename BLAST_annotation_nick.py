# BLAST ANNOTATION
#
# Make an annotation file for your phylogeny for FigTree based on a BLAST against SWISS-PROT
#
# Requires a BLAST output file against SWISS-PROT with -max_target_seqs 1 -max_hsps 1 -outfmt 6
#
# IMPORTANT: the blast needs to be done against a properly formatted swissprot database
#    To make this:
#    1. Download swiss-prot and unzip
#        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#        gunzip uniprot_sprot.fasta.gz
#    2. Replace spaces with '@'
#        sed -i 's/ /@/g' uniprot_sprot.fasta
#    3. Make a blast database
#        makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot.fasta.db -dbtype prot
#    4. Run the blast
#        blastp -query proteinA.fasta -db uniprot_sprot.fasta.db -outfmt 6 -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 -out ProteinA.fasta.blastout
#    5. Make the annotation file
#        python blast_annotation.py proteinA.fasta proteinA.fasta.blastout
#    6. Apply it to your phylogeny - open your phylogeny in FigTree, go to File>Import Annotations, and select your ProteinA.fasta.annotation file.
#
# Usage: python blast_annotation.py [fasta_file] [blast_output]
# Example: python blast_annotation.py proteinA.fasta proteinA.fasta.blastout
#

import sys

# 1. load in fasta file
fasta = open(sys.argv[1],'r').read()
seq_d = {}
for line in fasta.split('>')[1:]:
    seq_d[line.split('\n')[0].strip().strip('>')] = []

# 2. load in blast file
blast = open(sys.argv[2],'r').readlines()
blast_d = {}
evalue_d = {}
for line in blast:
    blast_d[line.split('\t')[0].strip()] = line.split('\t')[1].split('OS=')[0].split('@',1)[1].strip('@').replace('@','_')
    evalue_d[line.split('\t')[0].strip()] = line.split('\t')[-2].strip()

# 3. Output annotation
out = open(sys.argv[1]+'.annotation','w')
out.write('seq\tprotein\tevalue\n')
for s in seq_d.keys():
    try:
        out.write(s + '\t' + blast_d[s] + '\t' + evalue_d[s] + '\n')
    except:
        out.write(s + '\tno_annotation\tNA\n')
out.close()

