## Initial Setup
# Example: julia PATH/ organism_name lineage_dataset NO_PREY [optional: steps...]

path = ARGS[1] # path of the two raw sequence reads file
if path[end] != "/"
    path *= "/"
end

organism_name = ARGS[2] # name of the organism under study

lineage_dataset = ARGS[3] # e.g. eukaryota_odb10, alveolata_odb10

prey = ARGS[4] # e.g. Procryptobia, Spumella

what_to_run = ARGS[5:end]

function remove_name_hits(dir_path, name, blobtools_view_file, transcripts_file)
    contigs_list = dir_path * "no_" * name * "_contigs.list"
    run(pipeline(`grep -Ev $name $blobtools_view_file`, `cut -f 1`, contigs_list))

    # the following imitates lookup.pl

    open(contigs_list, "r") do f
        lines_contigs = readlines(f)
    end

    open(transcripts_file, "r") do f
        lines_transcripts = split(read(f), ">"; keepempty = false)
    end

    clean_transcripts_file = transcripts_file[1:end-length(".fasta")] * "_no_" * name * ".fasta"

    open(clean_transcripts_file, "w") do f
        for line_transcripts in lines_transcripts
            if any(line_contigs -> occursin(line_contigs, line_transcripts), lines_contigs)
                write(f, ">" * line_transcripts)
            end
        end
    end

    return clean_transcripts_file
end

function rename_transdecoder(dir_path, organism_name)
    for f in readdir(dir_path)
        if endswith(f, ".transdecoder.pep")

            open(dir_path * f, "r") do f
                lines = split(read(f), ">"; keepempty = false)
            end

            open(dir_path * "renamed_" * f, "w") do f
                for (i, line) in enumerate(lines)
                    write(f, ">" * organism_name * "_" * str(i) * "\n" * replace(split(line, "\n"; limit = 2)[2], "*" => ""))
                end
            end

        end
    end
end

function main(path, organism_name, lineage_dataset, prey, what_to_run)

    check_fastqc = false
    check_cutadapt = false
    check_rnaspades = false
    check_busco = false
    check_blastn_megablast = false
    check_diamond_blastx = false
    check_bowtie2 = false
    check_blobtools = false
    check_contamination_removal = false
    check_prey_removal = false
    check_transdecoder = false

    if isempty(what_to_run)
        check_fastqc = true
        check_cutadapt = true
        check_rnaspades = true
        check_busco = true
        check_blastn_megablast = true
        check_diamond_blastx = true
        check_bowtie2 = true
        check_blobtools = true
        check_contamination_removal = true
        check_prey_removal = true
        check_transdecoder = true
    else
        for step in what_to_run
            if step == "fastqc"
                check_fastqc = true
            elseif step == "cutadapt"
                check_cutadapt = true
            elseif step == "rnaspades"
                check_rnaspades = true
            elseif step == "busco"
                check_busco = true
            elseif step == "blastn_megablast"
                check_blastn_megablast = true
            elseif step == "diamond_blastx"
                check_diamond_blastx = true
            elseif step == "bowtie2"
                check_bowtie2 = true
            elseif step == "blobtools"
                check_blobtools = true
            elseif step == "contamination_removal"
                check_contamination_removal = true
            elseif step == "prey_removal"
                check_prey_removal = true
            elseif step == "transdecoder"
                check_transdecoder = true
            else
                error("Step ", step, " not recognize.")
            end
        end
    end





    ## 1. TRIM ADAPTERS/PRIMERS FROM RAW READS

    # use FASTQC for quality analysis

    fastqc_dir = path * "fastqc_raw_reads/"
    raw_seq_reads_1 = path * organism_name * "_R1_001.fastq.gz"
    raw_seq_reads_2 = path * organism_name * "_R2_001.fastq.gz"
    if check_fastqc
        if !isdir(fastqc_dir)
            mkdir(fastqc_dir)
        end
        run(`fastqc -o $fastqc_dir $raw_seq_reads_1 $raw_seq_reads_2`)
    end

    # use CUTADAPT to remove adapters from paired-end reads

    cutadapt_dir = path * "cutadapt/"
    cutadapt_output_file = cutadapt_dir * organism_name * "_R1_001.cutadapt.fastq"
    cutadapt_paired_output_file = cutadapt_dir * organism_name * "_R2_001.cutadapt.fastq"
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

    # use SPADES to assemble

    rnaspades_dir = path * "../" * organism_name * "_rnaspades/"
    if check_rnaspades
        if !isdir(rnaspades_dir)
            mkdir(rnaspades_dir)
        end
        run(`/opt/SPAdes-3.13.0-Linux/bin/spades.py
            --rna
            --pe1-1 $cutadapt_output_file
            --pe1-2 $cutadapt_paired_output_file
            --threads 24
            --memory 100
            -o $rnaspades_dir`)
    end
    transcripts_file = rnaspades_dir * "transcripts.fasta"

    # use BUSCO to estimate transcriptome coverage

    if check_busco
        path_ = pwd()
        cd("/Data/victoria/db/")
        run(`busco
            --in $transcripts_file
            --out $(organism_name * "_transcripts.fasta_BUSCO")
            --out_path $rnaspades_dir
            --lineage_dataset $lineage_dataset
            --mode transcriptome`)
        cd(path_)
    end





    ## 3. LOOK FOR CONTAMINATION

    contamination_dir = rnaspades_dir * organism_name * "_contamination_removal/"
    if !isdir(contamination_dir)
        mkdir(contamination_dir)
    end

    # use BLASTN to perform a megablast with nt database

    blastn_output_file = contamination_dir * "transcripts_vs_nt.blastn"
    if check_blastn_megablast
        run(`blastn
            -task megablast
            -query $transcripts_file
            -db /Data/databases/NCBI_NT/nt
            -outfmt 6 qseqid staxids bitscore std ssciorganism_names sskingdoms stitle
            -culling_limit 5
            -num_threads 24
            -evalue 1e-25
            -max_target_seqs 5
            -out $blastn_output_file`)
    end

    # use DIAMOND BLASTX

    blastx_output_file = contamination_dir * "transcripts.fasta.vs.uniprot_ref.mts1.1e25.out"
    if check_diamond_blastx
        run(`/Data/victoria/bin/diamond blastx
            --query $transcripts_file
            --max-target-seqs 1
            --sensitive
            --threads 24
            --db /Data/databases/uniprot_ref_v2020.06_diamond/uniprot_ref_proteomes.diamond.dmnd
            --evalue 1e-25
            --outfmt 6
            --out $blastx_output_file`)
    end

    # use BOWTIE2

    bowtie2_prefix = organism_name * "_transcripts.fasta"
    bowtie2_output_file = contamination_dir * organism_name * "_rnaspades.sam"
    if check_bowtie2
        path_ = pwd()
        cd(rnaspades_dir)
        run(`bowtie2-build $transcripts_file $bowtie2_prefix`)
        run(`bowtie2
            -x $bowtie2_prefix
            -U $cutadapt_output_file, $cutadapt_paired_output_file
            -S $bowtie2_output_file
            -p 24`)
        cd(path_)
    end

    # use BLOBTOOLS

    blobtools_taxify_output_file = contamination_dir * "transcripts.fasta.vs.uniprot_ref.mts1.1e25.taxified.out"
    blobtools_prefix = organism_name * "_rnaspades"
    blobtools_map2cov_output_file = contamination_dir * blobtools_prefix * ".sam.cov"
    blobtools_create_output_file = contamination_dir * blobtools_prefix * ".blobDB.json"
    if check_blobtools
        path_ = pwd()
        cd(contamination_dir)
        run(`/opt/blobtools/blobtools taxify
            -f $blobtools_taxify_output_file
            -m /Data/databases/uniprot_ref_v2020.06_diamond/uniprot_ref_proteomes.taxids
            -s 0
            -t 2`)
        run(`/opt/blobtools/blobtools map2cov -i $transcripts_file -s $bowtie2_output_file`)
        run(`/opt/blobtools/blobtools create
            -i $transcripts_file
            -t $blastn_output_file
            -t $blobtools_taxify_output_file
            -c $blobtools_map2cov_output_file
            -o $blobtools_prefix`)
        run(`/opt/blobtools/blobtools blobplot
            --infile $blobtools_create_output_file
            --rank family`)
        run(`/opt/blobtools/blobtools blobplot
            --infile $blobtools_create_output_file
            --rank phylum`)
        run(`/opt/blobtools/blobtools blobplot
            --infile $blobtools_create_output_file
            --rank superkingdom`)
        run(`/opt/blobtools/blobtools view
            --input $blobtools_create_output_file
            --out taxonomy
            --rank all
            --taxrule bestsum`)
        cd(path_)
    end
    blobtools_view_output_file = contamination_dir * "taxonomy." * blobtools_prefix * ".blobDB.bestsum.table.txt"

    # clean transcripts

    clean_transcripts_file = transcripts_file

    if check_contamination_removal
        clean_transcripts_file = remove_name_hits(contamination_dir, "Chordata", blobtools_view_output_file, clean_transcripts_file)
        clean_transcripts_file = remove_name_hits(contamination_dir, "Bacteria", blobtools_view_output_file, clean_transcripts_file)
        clean_transcripts_file = remove_name_hits(contamination_dir, "Arthropoda", blobtools_view_output_file, clean_transcripts_file)
    end

    if check_prey_removal
        if prey == "Procryptobia" # Colp34 data need to remove prey transcripts Procryptobia soronki
            prey_db = "/Data/victoria/transcriptomes/Procryptobia/Procryptobia_both.fa.DB"
            isvalidprey = true
        elseif prey == "Spumella" # Psammosa pacifica data need to remove prey transcripts Spumella elongata
            prey_db = "/Data/victoria/transcriptomes/Spumella_elongata_MMETSP1098/Spumella_elongata_CCAP955_1_MMETSP1098_cds.fa.db"
            isvalidprey = true
        else
            isvalidprey = false
        end

        if isvalidprey
            # NOTE: change "_no_Chordata_no_Bacteria_no_Arthropoda_vs_" accordingly to the above
            prey_blastn_output_file = contamination_dir * organism_name * "_no_Chordata_no_Bacteria_no_Arthropoda_vs_" * prey * ".blastnout"

            run(`blastn
                -task megablast
                -query $clean_transcripts_file
                -db $prey_db
                -outfmt 6
                -num_threads 24
                -evalue 1e-25
                -max_target_seqs 1
                -out $prey_blastn_output_file`)

            prey_contigs_list = contamination_dir * prey * "_contigs.list"
            run(pipeline(`cut -f 1`, prey_blastn_output_file, prey_contigs_list))

            # the following imitates lookup_reverse.pl

            open(prey_contigs_list, "r") do f
                lines_prey_contigs = readlines(f)
            end

            open(clean_transcripts_file, "r") do f
                lines_clean_transcripts = split(read(f), ">"; keepempty = false)
            end

            clean_transcripts_file = clean_transcripts_file[1:end-length(".fasta")] * "_no_" * prey * ".fasta"

            open(clean_transcripts_file, "w") do f
                for line_clean_transcripts in lines_clean_transcripts
                    if !any(line_prey_contigs -> occursin(line_prey_contigs, line_clean_transcripts), lines_prey_contigs)
                        write(f, ">" * line_clean_transcripts)
                    end
                end
            end

        end
    end





    ## 4. TRANSLATE TO PEPTIDES

    if check_transdecoder
        # Identify the ORFs in the cleaned fasta assembly using ###TransDecoder:

        run(`perl5.18.2 /opt/TransDecoder-TransDecoder-v5.3.0/TransDecoder.LongOrfs -t $clean_transcripts_file`)
        # Outputs peptides in a new TransDecoder directory as a file named "long_orfs.pep".

        # blast the ORFs against the Swissprot database to find ORFs that have possible matches to known proteins:

        for f_ in readdir("/Data/databases/")
            if "uniprot_sprot_" in f_
                database = "/Data/databases/" * f_ * "/uniprot_sprot.fasta"
                break
            end
        end
        blastp_output = rnaspades_dir * organism_name * "_longest_orfs.blastout"
        run(`blastp
            -query $(clean_transcripts_file * ".transdecoder_dir/longest_orfs.pep")
            -db $database
            -evalue 1e-5
            -max_target_seqs 1
            -outfmt 6
            -out $blastp_output
            -num_threads 24`)

        # Run the final ORF annotation using TransDecoder and tell it to include the ORFs that were identified in the Blast search as being things that matched to lessen the chance that they are lost in TransDecoder:

        run(`perl5.18.2 /opt/TransDecoder-TransDecoder-v5.3.0/TransDecoder.Predict -t $clean_transcripts_file --retain_blastp_hits $blastp_output`)

        rename_transdecoder(rnaspades_dir, organism_name)
    end
end

main(path, organism_name, lineage_dataset, prey, what_to_run)
