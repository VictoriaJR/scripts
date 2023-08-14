module QIIME2_pipeline
using CSV
using DataFrames

include("QIIME2_minion.jl")
    export QIIME2_minion

# utils 

# create tsv formatted manifest file
function generate_manifest(dir_path, experiment_name, primer_fwd, primer_rev)
    manifest = dir_path * "/" * experiment_name * "_" * "manifest.tsv"
    open(manifest, "w") do io 
        title = "sample-id" * "\t" * "absolute-filepath" * "\t" *  "forward-primer-5-3" * "\t" * "reverse-primer-5-3" *  "\r" * "\n"
        write(io, title)
        for f in readdir(dir_path; join = false)
            if isfile(f) && endswith(f, ".fastq")
                suffix = length(".fastq")
                f_name = f[1:end-suffix]
                filepath = dir_path * "/" * f 
                write(io, f_name * "\t" * filepath * "\t" * primer_fwd * "\t" * primer_rev * "\r" * "\n")
            end
        end
    end
end

# replace OTU headers with taxonomic annotation
function replace_OTU_header_taxonomy(fasta_file, taxonomy_file, experiment_name)
    d = Dict{String, Tuple{String, String}}()  # Use a tuple to store DNA sequence and taxonomy
    for line in eachline(taxonomy_file)
        tmp = split(strip(line), "\t"; limit = 3, keepempty = false)
        if length(tmp) ≥ 2
            taxonomy_id = tmp[1]
            taxonomy_affiliation = tmp[2]
            d[taxonomy_id] = (taxonomy_affiliation, "")
        end
    end
    
    fasta_lines = split(read(fasta_file, String), ">"; keepempty = false)
    modified_seqs = Vector{String}()
    
    for (i, line) in enumerate(fasta_lines)
        id_seq = split(strip(line), "\n"; limit = 2)
        taxonomy_id = strip(id_seq[1])
        
        if haskey(d, taxonomy_id)
            taxonomy_affiliation, sequence = d[taxonomy_id]
            modified_header = ">" * experiment_name * "_OTU_" * string(i) * "_" * taxonomy_affiliation
            modified_seq = modified_header * "\n" * id_seq[2]
            push!(modified_seqs, modified_seq)
        else
            println("Warning: Taxonomy ID not found for sequence ", i)
        end
    end
    
    new_file = fasta_file * ".taxonomy"
    open(new_file, "w") do io
        write(io, join(modified_seqs, "\n"))
    end
    
    return new_file
end

function OTUs_per_SRA_experiment(feature_table)
    df = DataFrame(CSV.File(feature_table; delim = "\t", header = 2))
    n = size(df, 2)
    sample = view(names(df), 2:n)
    results = Vector{Pair{String,Int}}(undef, n-1)
    
    for i in 2:n
        x = 0
        for g in df[!, i]
            if !iszero(g)
                x = x + 1
            end
        end
        results[i-1] = sample[i-1] => x
    end       
end


end # end module
