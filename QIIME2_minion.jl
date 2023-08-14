"""
    QIIME2_minion(dir_path, experiment_name, gene_amplicon_name, primer_set, cluster_perc, steps)
Execute the qiime2 pipeline consisting of the steps: "create_manifest", "import_data", "trim_data", "derep_vsearch", "denovo_clus_OTU", "taxonomy_OTU", "16S_taxonomy_filter", "rename_OTU_tax" "phylogenetic_tree"
Inputs:
`dir_path` = path to folder containing fastq files)
`experiment_name` = e.g. David_may20_soil
`gene_amplicon_name` = e.g. RbcL
`primer_set` = Primer set name e.g. RbcL_1.
`cluster_perc = Cluster percent to use e.g. 97
` add info here for primer sets
`steps` = steps of the pipeline to run. Default to all.


"""


function QIIME2_minion(dir_path::AbstractString, experiment_name::AbstractString, gene_amplicon_name::AbstractString, primer_set::AbstractString, cluster_perc::AbstractString, steps=String[])
    if dir_path[end] != "/"
        dir_path *= "/"
    end
    isdir(dir_path) || return throw(ArgumentError(string("Input path ", dir_path, " is not a directory")))

    check_create_manifest = false
    check_import_data = false
    check_trim_data = false
    check_derep_vsearch = false
    check_denovo_clus_OTU = false
    check_taxonomy_OTU = false
    check_16S_taxonomy_filter = false
    check_taxonomy_filter = false
    check_rename_OTU_tax = false
    check_phylogenetic_tree = false


    if isempty(steps)
        check_create_manifest = true
        check_import_data = true
        check_trim_data = true
        check_derep_vsearch = true
        check_denovo_clus_OTU = true
        check_taxonomy_OTU = true
        check_16S_taxonomy_filter = true
        check_taxonomy_filter = true
        check_rename_OTU_tax = true
        check_phylogenetic_tree = true

    else
        for step in steps
            if step == "create_manifest"
                check_create_manifest = true
            elseif step == "import_data"
                check_import_data = true
            elseif step == "trim_data"
                check_trim_data = true
            elseif step == "derep_vsearch"
                check_derep_vsearch = true
            elseif step == "denovo_clus_OTU"
                check_denovo_clus_OTU = true
            elseif step == "taxonomy_OTU"
                if gene_amplicon_name == "RbcL"
                    ref_db_seqs = "/Data/victoria/Qiime2/RbcL_project/database/rbcL/rbcl_seqs_diat.barcode.qza"
                    ref_db_tax = "/Data/victoria/Qiime2/RbcL_project/database/rbcL/rbcl_seqs_diat.barcode.tax.qza"
                elseif gene_amplicon_name == "18S"
                    ref_db_seqs = joinpath(dir_path, "../ref_db/silva-138-99-tax.qza")
                    ref_db_tax = joinpath(dir_path, "../ref_db/silva-138-99-tax.qza")
                elseif gene_amplicon_name == "16S"
                    ref_db_seqs = ""
                    ref_db_tax = ""
                else
                    return throw(ArgumentError(string("gene_amplicon_name", gene_amplicon_name, " is not valid")))
                end
                check_taxonomy_OTU = true
            elseif step == "16S_taxonomy_filter"
                check_16S_taxonomy_filter = true
            elseif step == "rename_OTU_tax"
                check_rename_OTU_tax = true
            elseif step == "phylogenetic_tree"
                check_phylogenetic_tree = true
            else
                return throw(ArgumentError(string("Step ", step, " is not valid")))
            end
        end
    end

## Identify primer sets
    # RbcL set 1: 1500bp diatom specific amplicon yield
    if primer_set == "RbcL_1"
        primer_fwd = "AAGGAGAAATHAATGTCT" #DPrbcL1
        primer_rev = "AARCAACCTTGTGTAAGTCTC" #DPrbcL7
        
    # RbcL set 2: 1200bp diatom specific amplicon yield
    elseif primer_set == "RbcL_2"
        primer_fwd = "AAAAGTGACCGTTATGAATC" #NDrbcL2
        primer_rev = "CCAATAGTACCACCACCAAAT" #NDrbcL8       
    
    # 18S set (EukA & EukB primer set
    elseif primer_set == "18S"
        primer_fwd = "AACCTGGTTGATCCTGCCAGT"
        primer_rev = "TGATCCTTCTGCAGGTTCACCTAC"
    
    # 16S set 
    elseif primer_set == "16S"
        primer_fwd = ""
        primer_rev = ""
    else
        return throw(ArgumentError("primer_set not recognized"))
    end
    
## 1. Create manifest file
    
    if check_create_manifest
    generate_manifest(dir_path, experiment_name, primer_fwd, primer_rev)
    end

 ## 2. Import data

    fastq_dir = dir_path * "fastq/"
    manifest_file = joinpath(dir_path, experiment_name * "_manifest.tsv")
    demux_file = joinpath(dir_path, "demux.qza")
    if check_import_data
            run(`qiime tools import
            --type 'SampleData[SequencesWithQuality]'
            --input-path $manifest_file
            --output-path $demux_file
            --input-format SingleEndFastqManifestPhred33V2`)
            run(`qiime demux summarize
            --i-data $demux_file
            --o-visualization $(joinpath(dir_path, "single-end-demux.qzv"))`)
    end


## 3. Remove primer

    demux_file_trim = joinpath(dir_path, "demux_trim.qza")
    if check_trim_data
        run(`qiime cutadapt trim-single
        --i-demultiplexed-sequences $demux_file
        --p-front $primer_fwd
        --p-adapter $primer_rev
        --p-match-adapter-wildcards
        --p-match-read-wildcards
        --o-trimmed-sequences $demux_file_trim`)
    end


 ## 4. Dereplicate using vsearch 
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

clus_table = joinpath(dir_path, "table-dn-" * cluster_perc * ".qza")
clus_seqs = joinpath(dir_path, "rep-seqs-dn-" * cluster_perc * ".qza")
feature_table_out = joinpath(dir_path, "exported-feature-table/")
feature_table_biom = feature_table_out * "feature-table.biom"
feature_table_tsv = feature_table_out * "feature-table.tsv"
clus_seqs_out = joinpath(dir_path, "rep-seqs-dn-" * cluster_perc * "/")
clus_seqs_fasta = clus_seqs_out * "dna-sequences.fasta"
clust_perc = "0." * cluster_perc 
                                

    if check_denovo_clus_OTU
        run(`qiime vsearch cluster-features-de-novo
        --i-table $derep_table
        --i-sequences $derep_seqs
        --p-perc-identity $clust_perc
        --o-clustered-table $clus_table 
        --o-clustered-sequences $clus_seqs
        --p-threads 4`)
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
    tax_OTU = joinpath(dir_path,  "taxonomy.qza")
    tax_OTU_out = joinpath(dir_path, "taxonomy/")
    tax_OTU_tsv = tax_OTU_out * "taxonomy.tsv"
    ref_db_seqs = "/Data/victoria/Qiime2/RbcL_project/database/rbcL/rbcl_seqs_diat.barcode.qza"
    ref_db_tax = "/Data/victoria/Qiime2/RbcL_project/database/rbcL/rbcl_seqs_diat.barcode.tax.qza"
    if check_taxonomy_OTU
        run(`qiime feature-classifier classify-consensus-blast
        --i-query $clus_seqs
        --i-reference-reads $ref_db_seqs
        --i-reference-taxonomy $ref_db_tax
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
        if !check_16S_taxonomy_filter
            replace_OTU_header_taxonomy(clus_seqs_fasta, tax_OTU_tsv, experiment_name)
        end
    end


    
## 8. filter 16S data 
    clean_table = joinpath(dir_path, "table-no-bacteria-no-archaea.qza")
    clean_seqs = joinpath(dir_path, "sequences-with-phyla-no-bacteria-no-archaea.qza")
    clean_tax_table = joinpath(dir_path, "taxonomy-no-bacteria-no-archaea.qza")
    clean_seqs_out = joinpath(dir_path, "sequences-with-phyla-no-bacteria-no-archaea/")
    clean_seqs_fasta = clean_seqs_out * "dna-sequences.fasta"
    if check_16S_taxonomy_filter && gene_amplicon_name != "16S"
        run(`qiime taxa filter-table
        --i-table $clus_table
        --i-taxonomy $tax_OTU
        --p-exclude bacteria,archaea
        --o-filtered-table $clean_table`)
        run(`qiime taxa filter-seqs
        --i-sequences $clus_seqs
        --i-taxonomy $tax_OTU
        --p-include p__
        --p-exclude bacteria,archaea
        --o-filtered-sequences $clean_seqs`)
        run(`qiime feature-classifier classify-consensus-blast
        --i-query $clean_seqs
        --i-reference-reads $silva_ref_seqs
        --i-reference-taxonomy $silva_ref_tax
        --o-classification $clean_tax_table`)
        run(`qiime tools export
        --input-path $clean_seqs
        --output-path $clean_seqs_out`)
        println(run(`grep ">" -c $clean_seqs_fasta`))
    end

## 9. Rename headers 

    if check_rename_OTU_tax
        final_fasta = joinpath(clean_seqs_out, ".taxonomy")
        replace_OTU_header_taxonomy(clean_seqs_fasta, tax_OTU_tsv, experiment_name)
    end



## 10. Make phylogenetic tree
    final_fasta = clean_seqs_fasta
    fasta_file_import = final_fasta * ".qza"
    fasta_file_import = joinpath(dir_path, "final_fasta_file/")
    aligned_seqs = joinpath(dir_path, "aligned-rep-seqs.qza")
    masked_alignment = joinpath(dir_path, "masked-aligned-rep-seqs.qza")
    unrooted_tree = joinpath(dir_path, "unrooted-tree/")
    unrooted_tree_out = joinpath(dir_path, "")
    rooted_tree = joinpath(dir_path, "rooted-tree.qza")
    rooted_tree_out = joinpath(dir_path, "rooted-tree/")
    if check_phylogenetic_tree
        run(`qiime tools import
        --input-path $final_fasta
        --output-path $fasta_file_import
        --type 'FeatureData[Sequence]'`)
        run(`qiime phylogeny align-to-tree-mafft-fasttree
        --i-sequences $fasta_file_import 
        --p-n-threads 4
        --o-alignment $aligned_seqs
        --o-masked-alignment $masked_alignment
        --o-tree $unrooted_tree
        --o-rooted-tree $rooted_tree`)
        run(`qiime tools export \
        --input-path $unrooted_tree \
        --output-path $unrooted_tree_out`)
        run(`qiime tools export \
        --input-path $rooted_tree \
        --output-path $rooted_tree_out`)
    end


end #end function

 
        
    


