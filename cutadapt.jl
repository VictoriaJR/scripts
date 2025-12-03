function cutadapt(dir_path::AbstractString, organism::AbstractString)

    # Normalize directory path
    dir_path = endswith(dir_path, "/") ? dir_path : dir_path * "/"

    isdir(dir_path) ||
        throw(ArgumentError("Input path $dir_path is not a directory"))

    # Input read files
    raw_seq_reads_1 = dir_path * organism * "_R1_001.fastq.gz"
    raw_seq_reads_2 = dir_path * organism * "_R2_001.fastq.gz"

    # Output directory
    cutadapt_dir = dir_path * "cutadapt/"
    cutadapt_output_file = cutadapt_dir * organism * "_R1_001.cutadapt.fastq"
    cutadapt_paired_output_file = cutadapt_dir * organism * "_R2_001.cutadapt.fastq"

    # Ensure output directory exists
    isdir(cutadapt_dir) || mkdir(cutadapt_dir)

    # Adapter sets for R1 and R2
    adapters_R1 = [
        "AGATGTGTATAAGAGACAG",
        "AAGCAGTGGTATCAACGCAGAGT",
        "TGGTATCAACGCAGAGT",
        "TACTCTGCGTTGATACCACTGCTT",
        "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
        "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
        "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
        "ACGAGCATCAGCAGCATACGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN",
        "AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG",
        "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG",
        "ACGAGCATCAGCAGCATACGA"
    ]

    adapters_R2 = [
        "NAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCGTATGCTGCTGATGCTCGT",
        "AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG",
        "CATTGCGCAATCTGTCTCTTATACACATCTGACGCTGCCGACGA",
        "TCGTATGCTGCTGATGCTCGT",
        "AGATGTGTATAAGAGACAG",
        "AAGCAGTGGTATCAACGCAGAGT",
        "TGGTATCAACGCAGAGT",
        "TACTCTGCGTTGATACCACTGCTT",
        "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
        "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
        "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    ]

    # Build the cutadapt command
    cmd = `cutadapt`

    for a in adapters_R1
        cmd = `$cmd -a $a`
    end
    for a in adapters_R2
        cmd = `$cmd -A $a`
    end

    cmd = `$cmd -o $cutadapt_output_file -p $cutadapt_paired_output_file $raw_seq_reads_1 $raw_seq_reads_2`

    # Run cutadapt
    run(cmd)

    return nothing
end
