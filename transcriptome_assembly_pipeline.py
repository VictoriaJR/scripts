import sys
import os
import errno
import subprocess





## Layout

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'





## Initial Setup

path = sys.argv[1] # path of the two raw sequence reads file
organism_name = sys.argv[2] # name of the organism under study
lineage = sys.argv[3] # e.g. eukaryota_odb10, alveolata_odb10
prey_name = sys.argv[4] # e.g. Procryptobia, Spumella





## 1. TRIM ADAPTERS/PRIMERS FROM RAW READS

print(bcolors.HEADER + "\n1. TRIM ADAPTERS/PRIMERS FROM RAW READS" + bcolors.ENDC)

# use FASTQC for quality analysis

print(bcolors.OKBLUE + "     - starting FASTQC for quality analysis ..." + bcolors.ENDC)

fastqc_dir = path + "fastqc_raw_reads/"
if not os.path.exists(os.path.dirname(fastqc_dir)):
    try:
        os.makedirs(os.path.dirname(fastqc_dir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

raw_seq_reads_1 = path + organism_name + "_R1_001.fastq.gz"
raw_seq_reads_2 = path + organism_name + "_R2_001.fastq.gz"

#subprocess.call(["fastqc", "-o", fastqc_dir, raw_seq_reads_1, raw_seq_reads_2])

print(bcolors.OKGREEN + "          ... done! The output have been stored in " + fastqc_dir + bcolors.ENDC)

# use CUTADAPT to remove adapters from paired-end reads

print(bcolors.OKBLUE + "     - starting CUTADAPT to remove adapters from paired-end reads ..." + bcolors.ENDC)

cutadapt_dir = path + "cutadapt/"
if not os.path.exists(os.path.dirname(cutadapt_dir)):
    try:
        os.makedirs(os.path.dirname(cutadapt_dir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

cutadapt_output_file = cutadapt_dir + organism_name + "_R1_001.cutadapt.fastq"
cutadapt_paired_output_file = cutadapt_dir + organism_name + "_R2_001.cutadapt.fastq"

with open(cutadapt_dir + organism_name + "_report.text", "w") as report:
    subprocess.call(["cutadapt",
                        "-a", "AGATGTGTATAAGAGACAG",
                        "-a", "AAGCAGTGGTATCAACGCAGAGT",
                        "-a", "TACTCTGCGTTGATACCACTGCTT",
                        "-a", "ACTCTGCGTTGATACCACTGCTT",
                        "-a", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
                        "-a", "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
                        "-a", "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                        "-a", "AGATGTGTATAAGAGACAG",
                        "-A", "AAGCAGTGGTATCAACGCAGAGT",
                        "-A", "TACTCTGCGTTGATACCACTGCTT",
                        "-A", "ACTCTGCGTTGATACCACTGCTT",
                        "-A", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
                        "-A", "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
                        "-A", "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                        "-o", cutadapt_output_file,
                        "-p", cutadapt_paired_output_file,
                        raw_seq_reads_1, raw_seq_reads_2], stdout = report)

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + cutadapt_dir + bcolors.ENDC)





## 2. TRANSCRIPTOME ASSEMBLY

print(bcolors.HEADER + "\n2. TRANSCRIPTOME ASSEMBLY" + bcolors.ENDC)

# use SPADES to assemble

print(bcolors.OKBLUE + "     - starting SPADES to assemble ..." + bcolors.ENDC)

rnaspades_dir = path + "../" + organism_name + "_rnaspades/"
if not os.path.exists(os.path.dirname(rnaspades_dir)):
    try:
        os.makedirs(os.path.dirname(rnaspades_dir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

with open(rnaspades_dir + organism_name + "_spades_report.text", "w") as report:
    subprocess.call(["/opt/SPAdes-3.13.0-Linux/bin/spades.py",
                         "--rna",
                         "--pe1-1", cutadapt_output_file,
                         "--pe1-2", cutadapt_paired_output_file,
                         "--threads", "24",
                         "--memory", "100",
                         "-o", rnaspades_dir], stdout = report)

transcripts_file = rnaspades_dir + "transcripts.fasta"

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + rnaspades_dir + bcolors.ENDC)

# use BUSCO to estimate transcriptome coverage

print(bcolors.OKBLUE + "     - starting BUSCO to estimate transcriptome coverage ..." + bcolors.ENDC)

subprocess.call(["busco",
                    "--in", transcripts_file,
                    "--out", organism_name + "_transcripts.fasta_BUSCO",
                    "--out_path", rnaspades_dir,
                    "--lineage", lineage,
                    "--mode", "transcriptome"])

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + rnaspades_dir + organism_name + "_transcripts.fasta_BUSCO/" + bcolors.ENDC)





## 3. LOOK FOR CONTAMINATION

print(bcolors.HEADER + "\n3. LOOK FOR CONTAMINATION" + bcolors.ENDC)

contamination_dir = rnaspades_dir + organism_name + "_contamination_removal/"
if not os.path.exists(os.path.dirname(contamination_dir)):
    try:
        os.makedirs(os.path.dirname(contamination_dir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

print(bcolors.OKBLUE + "     - starting BLASTN to perform a megablast with nt database ..." + bcolors.ENDC)

blastn_output_file = contamination_dir + "transcripts_vs_nt.blastn"

with open(contamination_dir + organism_name + "_blastn_megablast_report.text", "w") as report:
    subprocess.call(["blastn",
                        "-task", "megablast",
                        "-query", transcripts_file,
                        "-db", "/Data/NCBI_NT/nt",
                        "-outfmt", "6 qseqid staxids bitscore std ssciorganism_names sskingdoms stitle",
                        "-culling_limit", "5",
                        "-num_threads", "24",
                        "-evalue", "1e-25",
                        "-max_target_seqs", "5",
                        "-out", blastn_output_file], stdout = report)

print(bcolors.WARNING + "     NOTE: disregard warning error suggesting to download a database to associate taxonomy names to the hits (not relevant for this step) - see " + contamination_dir + organism_name + "_blastn_megablast_report.text" + bcolors.ENDC)

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + contamination_dir + bcolors.ENDC)

print(bcolors.OKBLUE + "     - starting DIAMOND BLASTX ..." + bcolors.ENDC)

blastx_output_file = contamination_dir + "transcripts.fasta.vs.uniprot_ref.mts1.1e25.out"

with open(contamination_dir + "blastx_report.text", "w") as report:
    subprocess.call(["/Data/victoria/bin/diamond", "blastx",
                        "--query", transcripts_file,
                        "--max-target-seqs", "1",
                        "--sensitive",
                        "--threads", "24",
                        "--db", "/Data/databases/uniprot_diamond/uniprot_ref_proteomes.diamond.dmnd",
                        "--evalue", "1e-25",
                        "--outfmt", "6",
                        "--out", blastx_output_file], stdout = report)

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + contamination_dir + bcolors.ENDC)

# use BLOBTOOLS to taxify diamond results

print(bcolors.OKBLUE + "     - starting BLOBTOOLS to taxify diamond results ..." + bcolors.ENDC)

os.chdir(contamination_dir)

subprocess.call(["/opt/blobtools/blobtools", "taxify",
                    "-f", blastx_output_file,
                    "-m", "/Data/databases/uniprot_diamond/uniprot_ref_proteomes.taxids",
                    "-s", "0",
                    "-t", "2"])

blobtools_taxify_output_file = contamination_dir + "transcripts.fasta.vs.uniprot_ref.mts1.1e25.taxified.out"

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + contamination_dir + bcolors.ENDC)

# use BOWTIE2-BUILD

print(bcolors.OKBLUE + "     - starting BOWTIE2-BUILD ..." + bcolors.ENDC)

os.chdir(rnaspades_dir)

with open(rnaspades_dir + "bowtie2-build_report.text", "w") as report:
    subprocess.call(["bowtie2-build", transcripts_file, organism_name + "_transcripts.fasta"], stdout = report)

print(bcolors.OKGREEN + "          ... done! The results have been stored in " + rnaspades_dir + bcolors.ENDC)

# use BOWTIE2

print(bcolors.OKBLUE + "     - starting BOWTIE2 ..." + bcolors.ENDC)

bowtie2_output_file = contamination_dir + organism_name + "_rnaspades.sam"

with open(contamination_dir + "bowtie2_report.text", "w") as report:
    subprocess.call(["bowtie2",
                        "-x", organism_name + "_transcripts.fasta",
                        "-U", cutadapt_output_file + "," + cutadapt_paired_output_file,
                        "-S", bowtie2_output_file,
                        "-p", "24"], stdout = report)

print(bcolors.OKGREEN + "          ... done! The result has been stored in " + contamination_dir + bcolors.ENDC)

# use BLOBTOOLS to parse sam file for coverage

os.chdir(contamination_dir)

print(bcolors.OKBLUE + "     - starting BLOBTOOLS to parse sam file for coverage ..." + bcolors.ENDC)

subprocess.call(["/opt/blobtools/blobtools", "map2cov",
                    "-i", transcripts_file,
                    "-s", bowtie2_output_file])

blobtools_map2cov_output_file = contamination_dir + organism_name + "_rnaspades.sam.cov"

print(bcolors.OKGREEN + "          ... done! The result has been stored in " + contamination_dir + bcolors.ENDC)

# use BLOBTOOLS to create and plot database

print(bcolors.OKBLUE + "     - starting BLOBTOOLS to create and plot database ..." + bcolors.ENDC)

subprocess.call(["/opt/blobtools/blobtools", "create",
                    "-i", transcripts_file,
                    "-t", blastn_output_file,
                    "-t", blobtools_taxify_output_file,
                    "-c", blobtools_map2cov_output_file,
                    "-o", organism_name + "_rnaspades"])

blobtools_create_output_file = contamination_dir + organism_name + "_rnaspades.blobDB.json"

subprocess.call(["/opt/blobtools/blobtools", "blobplot",
                    "--infile", blobtools_create_output_file,
                    "--rank", "family"])

subprocess.call(["/opt/blobtools/blobtools", "blobplot",
                    "--infile", blobtools_create_output_file,
                    "--rank", "phylum"])

subprocess.call(["/opt/blobtools/blobtools", "blobplot",
                    "--infile", blobtools_create_output_file,
                    "--rank", "superkingdom"])

subprocess.call(["/opt/blobtools/blobtools", "view",
                    "--input", blobtools_create_output_file,
                    "--out", "taxonomy",
                    "--rank", "all",
                    "--taxrule", "bestsum"])

blobtools_view_output_file = contamination_dir + "taxonomy." + organism_name + "_rnaspades.blobDB.bestsum.table.txt"

print(bcolors.OKGREEN + "          ... done! The result has been stored in " + contamination_dir + bcolors.ENDC)

os.chdir(rnaspades_dir)
# Remove Chordate hits - if its nucleotide its likely that all of these contigs found are homo derived.
# if these were protiens than it could somehow be a paralog of some sort
no_chor_contigs_list = contamination_dir + "no_chor_contigs.list"

with open(no_chor_contigs_list, "w") as out_list:
    grep = subprocess.Popen(["grep", "-Ev", "Chordata", blobtools_view_output_file], stdout = subprocess.PIPE)
    subprocess.Popen(["cut", "-f", "1"], stdin = grep.stdout, stdout = out_list)
    grep.stdout.close()

subprocess.call(["perl", "/Data/victoria/scripts/lookup.pl", transcripts_file, no_chor_contigs_list])

no_chor_transcripts_fasta_file = rnaspades_dir + organism_name + "_no_chor.fasta"

subprocess.call(["mv", rnaspades_dir + "lookup_out.fasta", no_chor_transcripts_fasta_file])

no_bac_contigs_list = contamination_dir + "no_bac_contigs.list"

with open(no_bac_contigs_list, "w") as out_file:
    grep = subprocess.Popen(["grep", "-Ev", "Bacteria", blobtools_view_output_file], stdout = subprocess.PIPE)
    subprocess.Popen(["cut", "-f", "1"], stdin = grep.stdout, stdout = out_file)
    grep.stdout.close()

subprocess.call(["perl", "/Data/victoria/scripts/lookup.pl", no_chor_transcripts_fasta_file, no_bac_contigs_list])

no_chor_no_bac_transcripts_fasta_file = rnaspades_dir + organism_name + "_no_chor_no_bac.fasta"

subprocess.call(["mv", rnaspades_dir + "lookup_out.fasta", no_chor_no_bac_transcripts_fasta_file])

no_arthro_contigs_list = contamination_dir + "no_arthro_contigs.list"

with open(no_arthro_contigs_list, "w") as out_file:
    grep = subprocess.Popen(["grep", "-Ev", "Arthropoda", blobtools_view_output_file], stdout = subprocess.PIPE)
    subprocess.Popen(["cut", "-f", "1"], stdin = grep.stdout, stdout = out_file)
    grep.stdout.close()

subprocess.call(["perl", "/Data/victoria/scripts/lookup.pl", no_chor_no_bac_transcripts_fasta_file, no_arthro_contigs_list])

no_chor_no_bac_no_arthro_transcripts_fasta_file = rnaspades_dir + organism_name + "_no_chor_no_bac_no_arthro.fasta"

subprocess.call(["mv", rnaspades_dir + "lookup_out.fasta", no_chor_no_bac_no_arthro_transcripts_fasta_file])

clean_nucleotides_fasta_file = no_chor_no_bac_no_arthro_transcripts_fasta_file

if prey_name == "Procryptobia": # Colp34 data need to remove prey transcripts Procryptobia soronki

    print(bcolors.OKBLUE + "     - starting BLASTN of cleaned fasta file to remove Procryptobia soronki prey..." + bcolors.ENDC)

    Procryptobia_blastn_output_file = contamination_dir + organism_name + "_no_chor_no_bac_vs_Procryptobia.blastnout"

    with open(contamination_dir + organism_name + "_vs_Procryptobia_blastn_report.txt", "w") as report:
        subprocess.call(["blastn",
                            "-task", "megablast",
                            "-query", no_chor_no_bac_transcripts_fasta_file,
                            "-db", "/Data/victoria/transcriptomes/Procryptobia/Procryptobia_both.fa.DB",
                            "-outfmt", "6",
                            "-num_threads", "24",
                            "-evalue", "1e-25",
                            "-max_target_seqs", "1",
                            "-out", Procryptobia_blastn_output_file], stdout = report)

    Kinetoplastid_contigs_list = contamination_dir + "Kinetoplastid_contigs.list"

    with open(Kinetoplastid_contigs_list, "w") as out_file:
        subprocess.call(["cut", "-f", "1", Procryptobia_blastn_output_file], stdout = out_file)

    subprocess.call(["perl", "/Data/victoria/scripts/lookup_reverse.pl", no_chor_no_bac_transcripts_fasta_file, Kinetoplastid_contigs_list])

    clean_nucleotides_fasta_file = rnaspades_dir + organism_name + "_no_chor_no_bac_no_Procryptobia.fasta"

    subprocess.call(["mv", rnaspades_dir + "lookup_out.fasta", clean_nucleotides_fasta_file])

    print(bcolors.OKGREEN + "          ... done! The result has been stored in " + contamination_dir + bcolors.ENDC)

elif prey_name == "Spumella": # Psammosa pacifica data need to remove prey transcripts Spumella elongata

    print(bcolors.OKBLUE + "     - starting BLASTN of cleaned fasta file to remove Spumella elongata prey..." + bcolors.ENDC)

    Spumella_blastn_output_file = contamination_dir + organism_name + "_rnaspades_no_chor_no_bac_vs_Spumella.blastnout"

    with open(contamination_dir + organism_name + "_vs_Spumella_blastn_report.txt", "w") as report:
        subprocess.call(["blastn",
                            "-task", "megablast",
                            "-query", no_chor_no_bac_transcripts_fasta_file,
                            "-db", "/Data/victoria/transcriptomes/Spumella_elongata_MMETSP1098/Spumella_elongata_CCAP955_1_MMETSP1098_cds.fa.db",
                            "-outfmt", "6",
                            "-num_threads", "24",
                            "-evalue", "1e-25",
                            "-max_target_seqs", "1",
                            "-out", Spumella_blastn_output_file], stdout = report)

    Spumella_contigs_list = contamination_dir + "Spumella_contigs.list"

    with open(Spumella_contigs_list, "w") as out_file:
        subprocess.call(["cut", "-f", "1", Spumella_blastn_output_file], stdout = out_file)

    subprocess.call(["perl", "/Data/victoria/scripts/lookup_reverse.pl", no_chor_no_bac_transcripts_fasta_file, Spumella_contigs_list])

    clean_nucleotides_fasta_file = rnaspades_dir + organism_name + "_no_chor_no_bac_no_Spumella.fasta"

    subprocess.call(["mv", rnaspades_dir + "lookup_out.fasta", clean_nucleotides_fasta_file])

    print(bcolors.OKGREEN + "          ... done! The result has been stored in " + contamination_dir + bcolors.ENDC)

else:

    print(bcolors.FAIL + "     - unknown input prey name: " + prey_name + ". Could not remove prey transcripts." + bcolors.ENDC)





## 4. TRANSLATE TO PEPTIDES

print(bcolors.HEADER + "\n4. TRANSLATE TO PEPTIDES" + bcolors.ENDC)

# Create a Blast database using Swissprot proteins to blast cleaned ORFs against.

print(bcolors.OKBLUE + "     - starting MAKEBLASTDB ..." + bcolors.ENDC)

blast_db_dir = rnaspades_dir + "blast_db/"
if not os.path.exists(os.path.dirname(blast_db_dir)):
    try:
        os.makedirs(os.path.dirname(blast_db_dir))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

subprocess.call(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", "-P", blast_db_dir])
subprocess.call(["gunzip", blast_db_dir + "uniprot_sprot.fasta.gz"])
subprocess.call(["makeblastdb",
                    "-in", blast_db_dir + "uniprot_sprot.fasta",
                    "-dbtype", "prot"])

print(bcolors.OKGREEN + "          ... done! The result has been stored in " + rnaspades_dir + bcolors.ENDC)

# Identify the ORFs in the cleaned fasta assembly using ###TransDecoder:

subprocess.call(["perl5.18.2", "/opt/TransDecoder-TransDecoder-v5.3.0/TransDecoder.LongOrfs",
                    "-t", clean_nucleotides_fasta_file])

###Outputs peptides in a new TransDecoder directory as a file named "long_orfs.pep".

# blast the ORFs against the Swissprot database to find ORFs that have possible matches to known proteins:

print(bcolors.OKBLUE + "     - starting BLASTP ..." + bcolors.ENDC)

blastp_output = rnaspades_dir + organism_name + "_longest_orfs.blastout"

subprocess.call(["blastp",
                    "-query", clean_nucleotides_fasta_file + ".transdecoder_dir/longest_orfs.pep",
                    "-db", blast_db_dir + "uniprot_sprot.fasta",
                    "-evalue", "1e-5",
                    "-max_target_seqs", "1",
                    "-outfmt", "6",
                    "-out", blastp_output,
                    "-num_threads", "24"])

print(bcolors.OKGREEN + "          ... done! The result has been stored in " + rnaspades_dir + bcolors.ENDC)

# Run the final ORF annotation using TransDecoder and tell it to include the ORFs that were identified in the Blast search as being things that matched to lessen the chance that they are lost in TransDecoder:

print(bcolors.OKBLUE + "     - starting TransDecoder ..." + bcolors.ENDC)

subprocess.call(["perl5.18.2", "/opt/TransDecoder-TransDecoder-v5.3.0/TransDecoder.Predict",
                    "-t", clean_nucleotides_fasta_file,
                    "--retain_blastp_hits", blastp_output])

for file_ in os.listdir(rnaspades_dir):
    if file_.endswith(".transdecoder.pep"):
        with open(rnaspades_dir + file_) as old_file, open(rnaspades_dir + "renamed_" + file_, "w") as new_file:
            i = 0
            for line in old_file:
                if line.startswith(">"):
                    i += 1
                    new_file.write(">" + organism_name + "_" + str(i) + "\n")
                elif line == "*":
                    continue
                else:
                    new_file.write(line.replace("*", ""))

print(bcolors.OKGREEN + "          ... done! The result has been stored in " + rnaspades_dir + bcolors.ENDC)
