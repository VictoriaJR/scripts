# For each gene in a given directory, extract the original taxa and coloured OTUs from tree.
# This script assumes that:
# - each [.fasta] file has the same prefix as a [.fasta.new.sl] file and a [.fasta.new.linsi.trimal.clean.treefile_coloured] file.
#   The single line format is important (.sl extension).
#   Use single_line_fasta.py on the appropriate [.fasta.new] file to create the corresponding [.fasta.new.sl] file.
# - each desired OTUs have been tagged with the same colour.
# For each gene, the output file name is the prefix of the [.fasta] file with the extension [.fasta.new.sl.cleaned].
#
# Usage: python extract_original_and_coloured_taxa.py [directory_path] [colour_code]
# Example: python extract_original_and_coloured_taxa.py /My/Path/ ff0000
# Note: the input [dirctory_path] must end with the character /

import os, sys

directory = sys.argv[1]
colour = sys.argv[2]

for file_ in os.listdir(directory):
    if file_.endswith(".fasta"):
        gene = file_[:-6]
        prefix = directory + gene
        new_taxa_fasta_file = prefix + ".fasta.new.sl"
        original_taxa_fasta_file = prefix + ".fasta"
        coloured_tree_file = prefix + ".fasta.new.linsi.trimal.clean.treefile_coloured"

        output_file = new_taxa_fasta_file + ".cleaned"
	
	with open(new_taxa_fasta_file) as new_taxa:
            new_taxa_lines = new_taxa.readlines()

        #

        coloured_taxa = []
        with open(coloured_tree_file) as tree:
            for line in tree:
                if colour in line:
                    coloured_taxa.append(line.split("[")[0].strip("\t").strip("'"))

	#

        with open(original_taxa_fasta_file) as original, open(output_file, "w") as output:
	    i = 0
            while i < len(new_taxa_lines):
                if new_taxa_lines[i].startswith(">"):
                    if new_taxa_lines[i].split(" ")[0].strip(">").strip("\n").replace("@", "_") in coloured_taxa or new_taxa_lines[i] in original:
                        output.write(new_taxa_lines[i].strip("\n") + "\n" + new_taxa_lines[i+1].strip("\n") + "\n")
                    i += 2
                else:
                    i += 1

