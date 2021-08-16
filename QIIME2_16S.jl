"""
    QIIME2_16S(dir_path, sequencing_technology, sequencing_type, primers, steps)
Execute the qiime2 pipeline consisting of the steps: "download_fastq", "import_data", "trim_data_1", "trim_data_2", "join_reads", "derep_vsearch", "denovo_clus_OTU97", "taxonomy_OTU", "mito_chloro_filter", "rename_OTU_tax_bioproj" 
Inputs:
`dir_path` = path to folder containing manifest file and sra_acc.txt file)
`bio_project` = e.g. "Massana_2015_PRJEB9133"
`sequencing_technology` = e.g. "illumina" or "pyrosequencing"
`sequencing_type` = e.g. "single_end" or "paired_end"
`primer_fwd` = e.g. ["forward_primer"]
`primer_rev` = e.g. ["reverse_primer"]
`steps` = steps of the pipeline to run. Default to all.


"""


function QIIME2_16S(dir_path::AbstractString, bio_project::AbstractString, sequencing_technology::AbstractString, sequencing_type::AbstractString, primer_fwd::AbstractString, primer_rev::AbstractString, steps=String[])
    if dir_path[end] != "/"
        dir_path *= "/"
    end
    isdir(dir_path) || return throw(ArgumentError(string("Input path ", dir_path, " is not a directory")))

    check_download_fastq = false
    check_import_data = false
    check_trim_data_1 = false
    check_trim_data_2 = false
    check_join_reads = false
    check_derep_vsearch = false
    check_denovo_clus_OTU97 = false
    check_taxonomy_OTU = false
    check_mito_chloro_filter = false
    check_rename_OTU_tax_bioproj = false

    if isempty(steps)
        check_download_fastq = true
        check_import_data = true
        check_trim_data_1 = true
        check_trim_data_2 = true
        check_derep_vsearch = true
        check_denovo_clus_OTU97 = true
        check_taxonomy_OTU = true
        check_mito_chloro_filter = true
        check_rename_OTU_tax_bioproj = true

    else
        for step in steps
            if step == "download_fastq"
                check_download_fastq = true
            elseif step == "import_data"
                check_import_data = true
            elseif step == "trim_data_1"
                check_trim_data_1 = true
            elseif step == "trim_data_2"
                check_trim_data_2 = true
            elseif step == "join_reads"
                check_join_reads = true
            elseif step == "derep_vsearch"
                check_derep_vsearch = true
            elseif step == "denovo_clus_OTU97"
                check_denovo_clus_OTU97 = true
            elseif step == "taxonomy_OTU"
                check_taxonomy_OTU = true
            elseif step == "16S_taxonomy_filter"
                check_mito_chloro_filter = true
            elseif step == "rename_OTU_tax_bioproj"
                check_rename_OTU_tax_bioproj = true
            else
                return throw(ArgumentError(string("Step ", step, " is not valid")))
            end
        end
    end




    ## 1. DOWNLOAD SRA

    fastq_dir = dir_path * "fastq/"
    sra_acc_file = dir_path * "sra_acc.txt"
    if check_download_fastq
        if !isdir(fastq_dir)
            mkdir(fastq_dir)
        end
        cd(fastq_dir)
        if sequencing_type == "single_end"
            for line in eachline(sra_acc_file)
                run(`/opt/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump $line --gzip`)
            end
        elseif sequencing_type == "paired_end"
            for line in eachline(sra_acc_file)
                run(`/opt/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump $line --split-files --gzip`)
            end
        end
        cd("../")
    end    
        

 ## 2. Import data

    manifest_file = joinpath(dir_path, bio_project * "_manifest.tsv")
    demux_file = joinpath(dir_path, "demux.qza")
    if check_import_data
        if sequencing_type == "single_end"
            run(`qiime tools import
            --type 'SampleData[SequencesWithQuality]'
            --input-path $manifest_file
            --output-path $demux_file
            --input-format SingleEndFastqManifestPhred33V2`)
            run(`qiime demux summarize
            --i-data $demux_file
            --o-visualization $(joinpath(dir_path, "single-end-demux.qzv"))`)
        elseif sequencing_type == "paired_end"
            run(`qiime tools import
            --type 'SampleData[SequencesWithQuality]'
            --input-path $manifest_file
            --output-path $demux_file
            --input-format PairedEndFastqManifestPhred64V2`)
            run(`qiime demux summarize
            --i-data $demux_file
            --o-visualization $(joinpath(dir_path, "paired-end-demux.qzv"))`)
        end
    end


## 3. Remove primer

    demux_file_trim = joinpath(dir_path, "demux_trim.qza")
    if check_trim_data_1
        if sequencing_type == "single_end"
                run(`qiime cutadapt trim-single
                --i-demultiplexed-sequences $demux_file
                --p-front $primer_fwd
                --p-match-adapter-wildcards
                --p-match-read-wildcards
                --o-trimmed-sequences $demux_file_trim`)
        elseif sequencing_type == "paired_end"
            run(`qiime cutadapt trim-paired
            --i-demultiplexed-sequences $demux_file
            --p-front-f $primer_fwd
            --p-match-adapter-wildcards
            --p-match-read-wildcards
            --o-trimmed-sequences $demux_file_trim`)
        else
            throw(DomainError)
        end
    elseif check_trim_data_2
        if sequencing_type == "single_end"
            run(`qiime cutadapt trim-single
            --i-demultiplexed-sequences $demux_file
            --p-front $primer_fwd
            --p-adapter $primer_rev
            --p-match-adapter-wildcards
            --p-match-read-wildcards
            --o-trimmed-sequences $demux_file_trim`)
        elseif sequencing_type == "paired_end"
            run(`qiime cutadapt trim-paired
            --i-demultiplexed-sequences $demux_file
            --p-front-f $primer_fwd
            --p-adapter-r $primer_rev
            --p-match-adapter-wildcards
            --p-match-read-wildcards
            --o-trimmed-sequences $demux_file_trim`)
        end
    end


## 4. Join reads together
## warning: replacing demux_file_trim with demux_file_join ; must run this to obtain proper demux file for downstream processes

    if check_join_reads
    demux_trim_join = joinpath(dir_path, "demux-joined.qza")
    demux_file_trim = demux_trim_join
        if sequencing_type == "paired_end"
            run(`qiime vsearch join-pairs
            --i-demultiplexed-seqs $demux_file
            --o-joined-sequences $demux_trim_join`)
        end
    end



 ## 5. Dereplicate using vsearch 
 derep_table = joinpath(dir_path, "table.qza")
 derep_seqs = joinpath(dir_path, "rep-seqs.qza")
    if check_derep_vsearch
        run(`qiime vsearch dereplicate-sequences
        --i-sequences $demux_file_trim
        --o-dereplicated-table $derep_table
        --o-dereplicated-sequences $derep_seqs`)
        run(`qiime metadata tabulate
        --m-input-file $derep_table
        --o-visualization $(joinpath(dir_path, "derep_table.qzv"))`)
    end


## 6. Cluster into OTUS 

clus_table = joinpath(dir_path, "table-dn-97.qza")
clus_seqs = joinpath(dir_path, "rep-seqs-dn-97.qza")
feature_table_out = joinpath(dir_path, "exported-feature-table/")
feature_table_biom = feature_table_out * "feature-table.biom"
feature_table_tsv = feature_table_out * "feature-table.tsv"
clus_seqs_out = joinpath(dir_path, "rep-seqs-dn-97/")
clus_seqs_fasta = clus_seqs_out * "dna-sequences.fasta"

    if check_denovo_clus_OTU97
        run(`qiime vsearch cluster-features-de-novo
        --i-table $derep_table
        --i-sequences $derep_seqs
        --p-perc-identity 0.97
        --o-clustered-table $clus_table 
        --o-clustered-sequences $clus_seqs
        --p-threads 20`)
        run(`qiime tools export
        --input-path $clus_table
        --output-path $feature_table_out`)
        cd(feature_table_out)
        run(`biom convert -i $feature_table_biom -o $feature_table_tsv --to-tsv`)
        run(`qiime tools export
        --input-path $clus_seqs
        --output-path $clus_seqs_out`)
        println(run(`grep ">" -c $clus_seqs_fasta`)) # * "OTUs"
        cd("../")
    end 


## 7. Annotate OTUs 
    silva_ref_seqs = joinpath(dir_path, "../ref_db/silva-138-99-seqs.qza")
    silva_ref_tax = joinpath(dir_path, "../ref_db/silva-138-99-tax.qza")
    tax_OTU = joinpath(dir_path,  "taxonomy.qza")
    tax_OTU_out = joinpath(dir_path, "taxonomy/")
    tax_OTU_tsv = tax_OTU_out * "taxonomy.tsv"
    if check_taxonomy_OTU
        run(`qiime feature-classifier classify-consensus-blast
        --i-query $clus_seqs
        --i-reference-reads $silva_ref_seqs
        --i-reference-taxonomy $silva_ref_tax
        --o-classification $tax_OTU`)
        run(`qiime metadata tabulate
        --m-input-file $tax_OTU
        --o-visualization $(joinpath(dir_path, "taxonomy.qzv"))`)
        run(`qiime taxa barplot
        --i-table $clus_table
        --i-taxonomy $tax_OTU
        --m-metadata-file $manifest_file
        --o-visualization $(joinpath(dir_path, "taxa-bar-plots.qzv"))`)
        run(`qiime tools export
        --input-path $tax_OTU
        --output-path $tax_OTU_out`) # file used for annotation of OTUs 
        replace_OTU_header_taxonomy(clus_seqs_fasta, tax_OTU_tsv, bio_project)
    end


    
## 8. filter out mitochondria and chloroplast data from 16S
    clean_table = joinpath(dir_path, "table-no-mito-no-chloro.qza")
    clean_table_out = joinpath(dir_path, "table-no-mito-no-chloro/")
    clean_table_out_biom = clean_table_out * "feature-table.biom"
    clean_table_out_tsv = clean_table_out * "feature-table.tsv"
    clean_seqs = joinpath(dir_path, "sequences-with-phyla-no-mito-no-chloro.qza")
    clean_tax_table = joinpath(dir_path, "taxonomy-no-mito-no-chloro.qza")
    clean_tax_table_out = joinpath(dir_path, "taxonomy-no-mito-no-chloro/")
    clean_tax_table_tsv = clean_tax_table_out * "taxonomy.tsv"
    clean_seqs_out = joinpath(dir_path, "sequences-with-phyla-no-mito-no-chloro/")
    clean_seqs_fasta = clean_seqs_out * "dna-sequences.fasta"
    if check_mito_chloro_filter
        run(`qiime taxa filter-table
        --i-table $clus_table
        --i-taxonomy $tax_OTU
        --p-exclude mitochondria,chloroplast
        --o-filtered-table $clean_table`)
        run(`qiime taxa filter-seqs
        --i-sequences $clus_seqs
        --i-taxonomy $tax_OTU
        --p-exclude mitochondria,chloroplast
        --o-filtered-sequences $clean_seqs`)
        run(`qiime feature-classifier classify-consensus-blast
        --i-query $clean_seqs
        --i-reference-reads $silva_ref_seqs
        --i-reference-taxonomy $silva_ref_tax
        --o-classification $clean_tax_table`)
        run(`qiime tools export
        --input-path $clean_tax_table
        --output-path $clean_tax_table_out`)
        run(`qiime tools export
        --input-path $clean_seqs
        --output-path $clean_seqs_out`)
        println(run(`grep ">" -c $clean_seqs_fasta`))
        run(`qiime tools export
        --input-path $clean_table
        --output-path $clean_table_out`)
        run(`biom convert -i $clean_table_out_biom -o $clean_table_out_tsv --to-tsv`)
    end

## 9. Rename headers and count number of OTUs per SRA experiment

    if check_rename_OTU_tax_bioproj
        replace_OTU_header_taxonomy(clean_seqs_fasta, clean_tax_table_tsv, bio_project)
        OTUs_per_SRA_experiment(clean_table_out_tsv)
    else
        return throw(ArgumentError(string("Rename and count OTU error")))
    end


end #end function

 
        
    


