"""
    transcriptome_assembly_paired_2_jezero(dir_path, organism, lineage_dataset, prey, steps)
Execute the transcriptome assembly pipeline consisting of the steps: "fastqc", "cutadapt", "rnaspades", "busco", "blastn_megablast", "diamond_blastx", "bowtie2", "blobtools", "contamination_removal", "busco_clean", "prey_removal", "transdecoder".
Inputs:
- `dir_path` = dir_path of the two raw sequence reads files
- `organism` = name of the organism under study. This will be used as a prefix throughout
- `lineage_dataset` = e.g. eukaryota_odb10, alveolata_odb10
- `contaminations` = e.g. ["Chordata", "Bacteria",  "Arthropoda"]
- `prey` = name of the prey, e.g. "Procryptobia" for Colp34, "Spumella" for Psammosa pacifica
- `steps` = steps of the pipeline to run. Default to all.
"""
function transcriptome_assembly_paired_2_jezero(dir_path::AbstractString, organism::AbstractString, lineage_dataset::AbstractString, contaminations, prey::AbstractString, steps=String[])
    if dir_path[end] != "/"
        dir_path *= "/"
    end
    isdir(dir_path) || return throw(ArgumentError(string("Input path ", dir_path, " is not a directory")))

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
    check_busco_clean = false
    check_transdecoder = false

    if isempty(steps)
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
        check_busco_clean = true
        check_transdecoder = true
    else
        for step in steps
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
                if prey == "Procryptobia"
                    prey_db = "/Data/victoria/transcriptomes/Procryptobia/Procryptobia_both.fa.DB"
                elseif prey == "Spumella"
                    prey_db = "/Data/victoria/transcriptomes/Spumella_elongata_MMETSP1098/Spumella_elongata_CCAP955_1_MMETSP1098_cds.fa.db"
                else
                    return throw(ArgumentError(string("Prey ", prey, " is not valid")))
                end
                check_prey_removal = true
            elseif step == "busco_clean"
                check_busco = true
            elseif step == "transdecoder"
                check_transdecoder = true
            else
                return throw(ArgumentError(string("Step ", step, " is not valid")))
            end
        end
    end





    ## 1. TRIM ADAPTERS/PRIMERS FROM RAW READS

    # use FASTQC for quality analysis

    fastqc_dir = dir_path * "fastqc_raw_reads/"
    raw_seq_reads_1 = dir_path * organism * "_R1_001.fastq.gz"
    raw_seq_reads_2 = dir_path * organism * "_R2_001.fastq.gz"
    if check_fastqc
        if !isdir(fastqc_dir)
            mkdir(fastqc_dir)
        end
        run(`fastqc -o $fastqc_dir $raw_seq_reads_1 $raw_seq_reads_2`)
    end

    # use CUTADAPT to remove adapters from paired-end reads
    # for this script i use the previously made cutadapt from the original assemblies
    cutadapt_dir = dir_path * "cutadapt/"
    cutadapt_output_file1 = cutadapt_dir * organism * "_S1_R1_001.cutadapt.fastq"
    cutadapt_paired_output_file1 = cutadapt_dir * organism * "_S1_R2_001.cutadapt.fastq"
    cutadapt_output_file2 = cutadapt_dir * organism * "_S2_R1_001.cutadapt.fastq"
    cutadapt_paired_output_file2 = cutadapt_dir * organism * "_S2_R2_001.cutadapt.fastq"
    if check_cutadapt
        if !isdir(cutadapt_dir)
            mkdir(cutadapt_dir)
        end
        run(`cutadapt
            -a AGATGTGTATAAGAGACAG
            -a AAGCAGTGGTATCAACGCAGAGT
            -a TGGTATCAACGCAGAGT
            -a TACTCTGCGTTGATACCACTGCTT
            -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
            -a CTGTCTCTTATACACATCTGACGCTGCCGACGA
            -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
            -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
            -A AGATGTGTATAAGAGACAG
            -A AAGCAGTGGTATCAACGCAGAGT
            -A TGGTATCAACGCAGAGT
            -A TACTCTGCGTTGATACCACTGCTT
            -A TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
            -A CTGTCTCTTATACACATCTGACGCTGCCGACGA
            -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
            -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
            -o $cutadapt_output_file
            -p $cutadapt_paired_output_file
            $raw_seq_reads_1 $raw_seq_reads_2`)
    end





    ## 2. TRANSCRIPTOME ASSEMBLY

    # use SPADES to assemble

    rnaspades_dir = dir_path * "../" * organism * "_rnaspades/"
    if check_rnaspades
        if !isdir(rnaspades_dir)
            mkdir(rnaspades_dir)
        end
        run(`/opt/spades_3.15.1/bin/rnaspades.py
            --pe1-1 $cutadapt_output_file1
            --pe1-2 $cutadapt_paired_output_file1
            --pe2-1 $cutadapt_output_file2
            --pe2-2 $cutadapt_paired_output_file2
            --threads 24
            -o $rnaspades_dir`)
    end
    transcripts_file = rnaspades_dir * "soft_filtered_transcripts.fasta"

    # use BUSCO to estimate transcriptome coverage

    if check_busco
        path_ = pwd()
        cd("/Data/victoria/db/")
        run(`busco
            --in $transcripts_file
            --out $(organism * "_transcripts.fasta_BUSCO")
            --out_path $rnaspades_dir
            --lineage_dataset $lineage_dataset
            -f
            --mode transcriptome`)
        cd(path_)
    end





    ## 3. LOOK FOR CONTAMINATION

    contamination_dir = rnaspades_dir * organism * "_contamination_removal/"
    if !isdir(contamination_dir)
        mkdir(contamination_dir)
    end

    # use BLASTN to perform a megablast with nt database

    blastn_output_file = contamination_dir * "transcripts_vs_nt.blastn"
    if check_blastn_megablast
        run(`/usr/bin/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin/blastn
            -task megablast
            -query $transcripts_file
            -db /Data/databases/NT_blast/nt
            -outfmt '6 qseqid staxids bitscore std ssciorganisms sskingdoms stitle'
            -culling_limit 5
            -num_threads 24
            -evalue 1e-25
            -max_target_seqs 5
            -out $blastn_output_file`)
    end

    # use DIAMOND BLASTX

    blastx_output_file = contamination_dir * "transcripts.fasta.vs.uniprot_ref.mts1.1e25.out"
    if check_diamond_blastx
        run(`diamond blastx
            --query $transcripts_file
            --max-target-seqs 1
            --sensitive
            --threads 24
            --db /Data/databases/uniprot_ref_diamond/uniprot_ref_proteomes.dmnd
            --evalue 1e-25
            --outfmt 6
            --out $blastx_output_file`)
    end

    # use BOWTIE2

    bowtie2_prefix = organism * "_transcripts.fasta"
    bowtie2_output_file = contamination_dir * organism * "_rnaspades.sam"
    if check_bowtie2
        path_ = pwd()
        cd(rnaspades_dir)
        run(`bowtie2-build $transcripts_file $bowtie2_prefix`)
        run(`bowtie2
            -x $bowtie2_prefix
            -U $cutadapt_output_file1,$cutadapt_paired_output_file1,$cutadapt_output_file2,$cutadapt_paired_output_file2
            -S $bowtie2_output_file
            -p 24`)
        cd(path_)
    end

    # use BLOBTOOLS

    blobtools_taxify_output_file = contamination_dir * "transcripts.fasta.vs.uniprot_ref.mts1.1e25.taxified.out"
    blobtools_prefix = organism * "_rnaspades"
    output_bam = contamination_dir * blobtools_prefix * ".bam"
    output_sorted_bam = contamination_dir * blobtools_prefix * ".sorted.bam"
    blobtools_map2cov_output_file = contamination_dir * blobtools_prefix * ".sorted.bam.cov"
    blobtools_create_output_file = contamination_dir * blobtools_prefix * ".blobDB.json"
    if check_blobtools
        path_ = pwd()
        cd(contamination_dir)
        run(`samtools view -@ 6 -b -S $bowtie2_output_file -o $output_bam`)
        run(`samtools sort -@ 6 $output_bam -o $output_sorted_bam`)
        run(`samtools index $output_sorted_bam`)
        run(`blobtools taxify
            -f $blastx_output_file
            -m /Data/databases/uniprot_ref_diamond/uniprot_ref_proteomes.taxids
            -s 0
            -t 2`)
        run(`blobtools map2cov -i $transcripts_file -b $output_sorted_bam`)
        run(`blobtools create
            -i $transcripts_file
            -t $blastn_output_file
            -t $blobtools_taxify_output_file
            -c $blobtools_map2cov_output_file
            -o $blobtools_prefix`)
        run(`blobtools plot
            --infile $blobtools_create_output_file
            --rank family`)
        run(`blobtools plot
            --infile $blobtools_create_output_file
            --rank phylum`)
        run(`blobtools plot
            --infile $blobtools_create_output_file
            --rank superkingdom`)
        run(`blobtools view
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
        for name in contaminations
            contigs_list = contamination_dir * "no_" * name * "_contigs.list"
            run(pipeline(`grep -Ev $name $blobtools_view_output_file`, `cut -f 1`, contigs_list))
            clean_transcripts_file = lookup_match(clean_transcripts_file, name, contigs_list)
        end
    end

    if check_prey_removal
        removed_contaminations = ""
        for name in contaminations
            removed_contaminations *= "no_" * name * "_" # e.g. "no_Chordata_no_Bacteria_no_Arthropoda_"
        end
        prey_blastn_output_file = contamination_dir * organism * "_" * removed_contaminations * "vs_" * prey * ".blastnout"
        run(`/usr/bin/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin/blastn
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
        clean_transcripts_file = lookup_mistmatch(clean_transcripts_file, prey, prey_contigs_list)
    end


    # use BUSCO to estimate transcriptome coverage

    if check_busco_clean
        path_ = pwd()
        cd("/Data/victoria/db/")
        run(`busco
            --in $clean_transcripts_file
            --out $(organism * "_transcripts.fasta_clean_BUSCO")
            --out_path $rnaspades_dir
            --lineage_dataset $lineage_dataset
            -f
            --mode transcriptome`)
        cd(path_)
    end



    ## 4. TRANSLATE TO PEPTIDES

    if check_transdecoder
        path_ = pwd()
        cd(rnaspades_dir)
        # Identify the ORFs in the cleaned fasta assembly using ###TransDecoder:
        clean_transcripts_file = /Data/victoria/psammosa/C34_2017_2020_rnaspades/soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Procryptobia.fasta
        run(`/opt/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $clean_transcripts_file`)
        # Outputs peptides in a new TransDecoder directory as a file named "long_orfs.pep".

        # blast the ORFs against the Swissprot database to find ORFs that have possible matches to known proteins:
        database = "/Data/databases/uniprot_sprot_blast/uniprot_sprot.fasta.blastDB"
        for f in readdir("/Data/databases/uniprot_sprot_blast/"; join = false)
            if startswith(f, "uniprot_sprot_")
                database *= f
                break
            end
        end
        blastp_output = rnaspades_dir * organism * "_longest_orfs.blastout"
        run(`/usr/bin/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin/blastp
            -query $(clean_transcripts_file * ".transdecoder_dir/longest_orfs.pep")
            -db $database
            -evalue 1e-5
            -max_target_seqs 1
            -outfmt 6
            -out $blastp_output
            -num_threads 24`)

        # Run the final ORF annotation using TransDecoder and tell it to include the ORFs that were identified in the Blast search as being things that matched to lessen the chance that they are lost in TransDecoder:

        run(`/opt/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict
            -t $clean_transcripts_file
            --retain_blastp_hits $blastp_output`)
        cd(path_)

        dir_rename_headers(rnaspades_dir, ".transdecoder.pep",  organism)
    end
end
