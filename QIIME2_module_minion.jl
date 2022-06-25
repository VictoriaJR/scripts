module QIIME2_pipeline
using CSV
using DataFrames

include("QIIME2_minion.jl")
    export QIIME2_minion

# utils 

# create tsv formatted manifest file
function generate_manifest(dir_path, forward_primer, reverse_primer)
    list_file_name = String[]
    list_file_fastq = String[]
    manifest = dir_path * "/" * "manifest.tsv"
    open(manifest, "w") do io 
        title = "sample-id" * "\t" * "absolute-filepath" * "\t" *  "forward-primer-5-3" * "\t" * "reverse-primer-5-3" *  "\r" * "\n"
        write(io, title)
        c = forward_primer
        d = reverse_primer
        for f in readdir(dir_path; join = false)
            if isfile(f) && endswith(f, ".fastq")
                suffix = length(".fastq")
                f_name = f[1:end-suffix]
                filepath = dir_path * "/" * f 
                write(io, f_name * "\t" * filepath * "\t" * forward_primer * "\t" * reverse_primer * "\r" * "\n")
            end
        end
    end
end

# replace OTU headers with taxonomic annotation
function replace_OTU_header_taxonomy(fasta_file, taxonomy_file, experiment_name)
    d = Dict{String, String}()
    for line in eachline(taxonomy_file)
        tmp = split(line, "\t"; limit = 3, keepempty = false)
        if length(tmp) ≥ 2
            taxonomy_id_1 = tmp[1]
            taxonomy_id_2 = tmp[2]
            push!(d, taxonomy_id_1 => taxonomy_id_2)
        end
    end
    fasta_lines = split(read(fasta_file, String), ">"; keepempty = false)
    new_file = fasta_file * ".taxonomy"
    open(new_file, "w") do io
        for (i,line) in enumerate(fasta_lines)
            id_seq = split(line, "\n"; limit = 2)
            write(io, ">" * bioproject * "_OTU_" * string(i) * "_" * d[strip(id_seq[1], '\r')] * "\n" * id_seq[2])
        end
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
