module Transcriptome

include("assembly_pipeline.jl")
    export transcriptome_assembly
include("assembly_pipeline_unpaired.jl")
    export transcriptome_assembly_unpaired
include("assembly_pipeline_paired_2_transcriptomes.jl")
    export transcriptome_assembly_paired_2_transcriptomes

## utils

"""
    lookup_match(transcripts, name, list)

Create a transcripts where `name` has been removed by keeping from `transcripts` everything in `list`.
"""
function lookup_match(transcripts, name, list)
	endswith(transcripts, ".fasta") || return throw(ArgumentError(string(transcripts, " must be a file ending with .fasta")))
	lines_list = readlines(list)
	lines_transcripts = split(read(transcripts, String), ">"; keepempty = false)
	clean_transcripts = transcripts[1:end-6] * "_no_" * name * ".fasta"
	open(clean_transcripts, "w") do f
		for line_transcripts in lines_transcripts
			header = split(line_transcripts, "\n"; limit = 2)[1]
			if any(line_list -> line_list == header, lines_list)
				write(f, ">" * line_transcripts)
			end
		end
	end
	return clean_transcripts
end

"""
    lookup_mismatch(transcripts, name, list)

Create a transcripts where `name` has been removed by keeping from `transcripts` everything not in `list`.
"""
function lookup_mismatch(transcripts, name, list)
	endswith(transcripts, ".fasta") || return throw(ArgumentError(string(transcripts, " must be a file ending with .fasta")))
	lines_list = readlines(list)
	lines_transcripts = split(read(transcripts, String), ">"; keepempty = false)
	clean_transcripts = transcripts[1:end-6] * "_no_" * name * ".fasta"
	open(clean_transcripts, "w") do f
		for line_transcripts in lines_transcripts
			header = split(line_transcripts, "\n"; limit = 2)[1]
			if !any(line_list -> line_list == header, lines_list)
				write(f, ">" * line_clean_transcripts)
			end
		end
	end
	return clean_transcripts
end

"""
    dir_rename_headers(dir, suffix, organism)

Given an organism name `organism`, rename the headers inside every files in the directory `dir` which end with `suffix`.
"""
function dir_rename_headers(dir, suffix, organism)
	renamed = String[]
	for f in readdir(dir; join = true)
		if endswith(f, suffix)
			push!(renamed, file_rename_headers(f, organism))
		end
	end
	return renamed
end

"""
    file_rename_headers(file, organism)

Given an organism name `organism`, rename the headers inside the file `file`.
"""
function file_rename_headers(file, organism)
	lines = split(read(file, String), ">"; keepempty = false)
	renamed = file * ".renamed"
	open(renamed, "w") do f
		for (i, line) in enumerate(lines)
			seq = replace(split(line, "\n"; limit = 2)[2], "*" => "")
			write(f, ">" * organism * "_" * string(i) * "\n" * seq)
		end
	end
	return renamed
end

function dir_rename_headers_phylo(dir, suffix, organism)
	renamed = String[]
	for f in readdir(dir; join = true)
		if endswith(f, suffix)
			push!(renamed, file_rename_headers_phylo(f, organism))
		end
	end
	return renamed
end

function file_rename_headers_phylo(file, organism)
	lines = split(read(file, String), ">"; keepempty = false)
	renamed = file * ".renamed"
	open(renamed, "w") do f
		for line in lines
			seq = replace(split(line, "\n"; limit = 2)[2], "*" => "")
			write(f, ">" * organism * "@" * split(line, " "; limit = 2)[1] * "\n" * seq)
		end
	end
	return renamed
end


"""
For each gene in a given directory, extract the original taxa and coloured OTUs from tree.
This script assumes that:
- each [.fasta] file has the same prefix as a [.fasta.new.sl] file and a [.fasta.new.linsi.trimal.clean.treefile_coloured] file.
  The single line format is important (.sl extension).
  Use single_line_fasta.py on the appropriate [.fasta.new] file to create the corresponding [.fasta.new.sl] file.
- each desired OTUs have been tagged with the same colour.
For each gene, the output file name is the prefix of the [.fasta] file with the extension [.fasta.new.sl.cleaned].
"""
function dir_extract_original_and_coloured_taxa(dir, colour)
	extracted = String[]
	for f in readdir(dir; join = true)
		if endswith(f, ".fasta")
			prefix = f[1:end-6]
			new_taxa = prefix * ".fasta.new.sl"
			coloured_tree = prefix * ".fasta.new.linsi.trimal.clean.treefile_coloured"
			push!(extracted, file_extract_original_and_coloured_taxa(f, new_taxa, coloured_tree, colour))
		end
	end
	return extracted
end

function file_extract_original_and_coloured_taxa(original, new, coloured_tree, colour)
	read_original = read(original, String)
	lines_new = split(read(new, String), ">"; keepempty = false)
	coloured = String[]
	open(coloured_tree, "r") do f
		for line in eachline(f)
			if occursin(colour, line)
				push!(coloured, strip(split(line, "["; limit = 2)[1], ['\t', '\'']))
			end
		end
	end
	clean_new = new * ".cleaned"
	open(clean_new, "w") do f
		for line_new in lines_new
			if any(x -> occursin(x, line_new), coloured) || occursin(line_new, read_original)
				write(f, ">" * line_new)
			end
		end
	end
	return clean_new
end

function single_line_fasta(file)
	lines = split(read(file, String), ">"; keepempty = false)
	out_file = file * ".sl"
	open(out_file, "w") do f
		for line in lines
			parts = split(line, "\n"; limit = 2)
			write(f, ">" * parts[1] * "\n" * replace(parts[2], "\n" => "") * "\n")
		end
	end
	return out_file
end

function combine_datasets(out_dir, in_dirs...)

	dict_seq = Dict{String,Vector{String}}() # creates a dictionnary according to the relation: file => list of protein sequence

	for dir in in_dirs # loop through the directories given as arguments
		for f in readdir(dir; join = false) # loop through each file within a directory
			f_out = joinpath(out_dir, f * ".out")
			close(open(f_out, "a")) # creates a new file if it does not already exist
			push!(dict_seq, f_out => [""])
		end
	end

	for dir in in_dirs # loop through the directories given as arguments
		for f in readdir(dir; join = false) # loop through each file within a directory
			lines = split(read(joinpath(dir, f), String), ">"; keepempty = false)
			f_out = joinpath(out_dir, f * ".out")
			open(f_out, "a") do io # append sequence name and sequence protein in the corresponding out file
				for line in lines
					parts = split(line, "\n"; limit = 2)
					seq_name = split(parts[1], " "; limit = 2)[1]
					seq = parts[2]
					seq_f_out = dict_seq[f_out]
					if !(seq in seq_f_out)
						write(io, ">" * seq_name * "\n") # write sequence name
						write(io, seq) # write protein sequence
						push!(seq_f_out, seq) # add the new protein sequence to the list corresponding to the in-file
					end
				end
			end
		end
	end

	return collect(keys(dict_seq)) # turns KeySet into Vector
end

	export lookup_match, lookup_mismatch, dir_rename_headers, file_rename_headers,
		dir_extract_original_and_coloured_taxa, file_rename_headers_phylo,
		dir_rename_headers_phylo, file_extract_original_and_coloured_taxa,
		single_line_fasta, combine_datasets

end  # module Transcriptome
