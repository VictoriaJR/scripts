"""
    transcriptome_assembly_db(dir_path, organism, steps)

Execute the transcriptome assembly pipeline consisting of the steps: "cutadapt", "flash2".
Inputs:
- `dir_path` = dir_path of the two raw sequence reads files
- `organism` = name of the organism under study. This will be used as a prefix throughout
"""
function transcriptome_assembly_db(dir_path::AbstractString, organism::AbstractString)
    if dir_path[end] != "/"
        dir_path *= "/"
    end
    isdir(dir_path) || return throw(ArgumentError(string("Input path ", dir_path, " is not a directory")))

    ## 1. TRIM ADAPTERS/PRIMERS FROM RAW READS

    raw_seq_reads_1 = dir_path * organism * "_1.fastq"
    raw_seq_reads_2 = dir_path * organism * "_2.fastq"


    # use CUTADAPT to remove adapters from paired-end reads

    cutadapt_dir = dir_path * "cutadapt/"
    cutadapt_output_file = cutadapt_dir * organism * "_1.cutadapt.fastq"
    cutadapt_paired_output_file = cutadapt_dir * organism * "_2.cutadapt.fastq"
    if check_cutadapt
        if !isdir(cutadapt_dir)
            mkdir(cutadapt_dir)
        end
        run(`cutadapt
            -a AGATGTGTATAAGAGACAG
            -a AAGCAGTGGTATCAACGCAGAGT
            -a AGATGTGTATAAGAGACAG
            -a AAGCAGTGGTATCAACGCAGAGT
            -a TACTCTGCGTTGATACCACTGCTT
            -a ACTCTGCGTTGATACCACTGCTT
            -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
            -a CTGTCTCTTATACACATCTGACGCTGCCGACGA
            -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
            -a AGATGTGTATAAGAGACAG
            -A AAGCAGTGGTATCAACGCAGAGT
            -A TACTCTGCGTTGATACCACTGCTT
            -A ACTCTGCGTTGATACCACTGCTT
            -A TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
            -A CTGTCTCTTATACACATCTGACGCTGCCGACGA
            -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
            -o $cutadapt_output_file
            -p $cutadapt_paired_output_file
            $raw_seq_reads_1 $raw_seq_reads_2`)
    end



    ## 2. TRANSCRIPTOME ASSEMBLY

    # use FLASH to assemble
        run(`/Data/victoria/software/FLASH2-2.2.00/flash2 flash
            -d $dir_path
            -t 12
            $cutadapt_output_file
            $cutadapt_paired_output_file`)
end
