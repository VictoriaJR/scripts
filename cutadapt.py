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


