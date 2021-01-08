# extract nucleotide sequences from a clean nucleotide file using a query

import sys
import os
import errno
import subprocess

path = sys.argv[1] # path of the cleaned nucleotide sequences
query_file = sys.argv[2] # query used to blast against nucleotides

for file_ in os.listdir(path):
    if file_.endswith(".fasta"):
        organism_name = file_[:-6] # name of the organism under study
        nucleotide_file = path + organism_name + ".fasta"

        database_dir = path + "nucl_databases/"
        if not os.path.exists(os.path.dirname(database_dir)):
            try:
                os.makedirs(os.path.dirname(database_dir))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        blast_db_file = database_dir + organism_name + ".fasta.DB"

#        subprocess.call(["makeblastdb",
#                            "-in", nucleotide_file,
#                            "-dbtype", "nucl",
#                            "-out", blast_db_file])

        # BLAST 18S QUERY AGAINST NEW DB (max target sequences = 15; e-value=25)

        blastn_dir = path + "blastn/"
        if not os.path.exists(os.path.dirname(blastn_dir)):
            try:
                os.makedirs(os.path.dirname(blastn_dir))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        blastn_output_file = blastn_dir + organism_name + "_18S_queries_e-25.blastnout"

#        subprocess.call(["blastn",
#                            "-query", query_file,
#                            "-db", blast_db_file,
#                            "-outfmt", "6",
#                            "-num_threads", "1",
#                            "-evalue", "1e-25",
#                            "-max_hsps", "15",
#                            "-out", blastn_output_file])

        #

        sequence_dir = path + "gene_sequences/"
        if not os.path.exists(os.path.dirname(sequence_dir)):
            try:
                os.makedirs(os.path.dirname(sequence_dir))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        gene_sequence_file = sequence_dir + organism_name + "_extracted_18S_queries.fasta"

        with open(nucleotide_file) as nucl_file:
            lines_nucl = nucl_file.readlines()
            len_lines_nucl = len(lines_nucl)

        with open(blastn_output_file) as blastn_file: # extract second column of blastn file
            sequences_names = []
            for line_blastn in blastn_file:
                sequence_name = line_blastn.split("\t")[1]
                if sequence_name not in sequences_names: # prevents duplicate
                    sequences_names.append(sequence_name)

        with open(gene_sequence_file, "w") as gene_file: # search the nucleotide file for the sequence name
            for sequence_name in sequences_names:
                for i, line_nucl in enumerate(lines_nucl):
                    if sequence_name in line_nucl:
                        gene_file.write(">" + organism_name + line_nucl.replace(">", "_"))
                        for j in range(i+1, len_lines_nucl):
                            next_line_nucl = lines_nucl[j]
                            if ">" not in next_line_nucl:
                                gene_file.write(next_line_nucl)
                                j += 1
                            else:
                                break
                        break
