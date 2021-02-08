module Transcriptome2

	using FASTX

##

function print_file_tree(root, maxdepth)
    @assert maxdepth ≥ -1 # -1 for all depth
    println(root)
    dirscount, filescount = print_file_tree(root, 0, [true], 0, 0, maxdepth)
    println(string("\n", dirscount, " directories, ", filescount, " files"))
    return nothing
end

function print_file_tree(root, depth, opendirs, dirscount, filescount, maxdepth)
    files = readdir(root)
    len = length(files)
    lastdepth = (depth == maxdepth)
    for (i, f) in enumerate(files)
        startswith(f, ".") && continue
        lastitem = (i == len)
        lastitem && (opendirs[end] = false)
        for p in opendirs[1:end-1]
            print(p ? "│   " : "    ")
        end
        println("$(lastitem ? "└" : "├")── " * f)
        path = joinpath(root, f)
        if isdir(path)
            dirscount += 1
            lastdepth && continue
            push!(opendirs, true)
            dirscount, filescount = print_file_tree(path, depth+1, opendirs, dirscount, filescount, maxdepth)
            pop!(opendirs)
        else
            filescount += 1
        end
    end
    return dirscount, filescount
end

	export print_file_tree

##

function convert_old_workflow_transcriptome_assembly(dir_path, organism)
	isdir(dir_path) || return throw(ArgumentError(string("Input path ", dir_path, " is not a directory")))
    dir_path = normpath(dir_path * "/")
	fastqc_dir = dir_path * "fastqc/"
	mv(dir_path * "fastqc_raw_reads/", fastqc_dir)
	cutadapt_dir = dir_path * "cutadapt/"
	transcriptome_dir = dir_path * "transcriptome/"
	mv(dir_path * "../" * organism * "_rnaspades/", transcriptome_dir)
	contamination_dir = transcriptome_dir * "contamination_removal/"
	mv(dir_path * "../" * organism * "_rnaspades/" * organism * "_contamination_removal/", contamination_dir)
	return fastqc_dir, cutadapt_dir, transcriptome_dir, contamination_dir
end

function setup_workflow_transcriptome_assembly(dir_path)
    isdir(dir_path) || return throw(ArgumentError(string("Input path ", dir_path, " is not a directory")))
    dir_path = normpath(dir_path * "/")
    fastqc_dir = dir_path * "fastqc/"
    if !isdir(fastqc_dir)
        mkdir(fastqc_dir)
    end
    cutadapt_dir = dir_path * "cutadapt/"
    if !isdir(cutadapt_dir)
        mkdir(cutadapt_dir)
    end
    transcriptome_dir = dir_path * "transcriptome/"
    if !isdir(transcriptome_dir)
        mkdir(transcriptome_dir)
    end
    contamination_dir = transcriptome_dir * "contamination_removal/"
    if !isdir(contamination_dir)
        mkdir(contamination_dir)
    end
    print_file_tree(dir_path, 0)
    return fastqc_dir, cutadapt_dir, transcriptome_dir, contamination_dir
end

	export convert_old_workflow_transcriptome_assembly, setup_workflow_transcriptome_assembly

##

"""
	single_line(outfile_name, transcripts_name)

Returns `outfile_name` which contains `transcripts_name` with sequences on a single line.
"""
function single_line(outfile_name, transcripts_name)
	transcripts = open(FASTA.Reader, transcripts_name)
	outfile = open(FASTA.Writer, outfile_name)

	record = FASTA.Record()
	while !eof(transcripts)
		read!(transcripts, record)
		write(outfile, record)
	end

	close(transcripts)
	close(outfile)

	return outfile_name
end

"""
    lookup(outfile_name, transcripts, list, method)

Returns a transcripts `outfile_name` from a transcripts `transcripts_name` and `list` according to a `method`.
"""
function lookup(outfile_name, transcripts_name, list, method)
	if method == :include
		return _lookup_include(outfile_name, transcripts_name, list)
	elseif method == :exclude
		return _lookup_exclude(outfile_name, transcripts_name, list)
	else
		return throw(ArgumentError("Input method not defined."))
	end
end

function _lookup_include(outfile_name, transcripts_name, list)
	lines_list = eachline(list)
	transcripts = open(FASTA.Reader, transcripts_name)
	outfile = open(FASTA.Writer, outfile_name)

	record = FASTA.Record()
	while !eof(transcripts)
	    read!(transcripts, record)
		id = identifier(record)
		if id in lines_list
			write(outfile, record)
		end
	end

	close(transcripts)
	close(outfile)

	return outfile_name
end

function _lookup_exclude(outfile_name, transcripts_name, list)
	lines_list = eachline(list)
	transcripts = open(FASTA.Reader, transcripts_name)
	outfile = open(FASTA.Writer, outfile_name)

	record = FASTA.Record()
	while !eof(transcripts)
	    read!(transcripts, record)
		id = identifier(record)
		if !(id in lines_list)
			write(outfile, record)
		end
	end

	close(transcripts)
	close(outfile)

	return outfile_name
end

"""
	extract_original_and_coloured_taxa(outfile, original, new, coloured_tree, colour)

Extract the original taxa and coloured OTUs from tree.
"""
function extract_original_and_coloured_taxa(outfile, original, new, coloured_tree, colour)
	coloured_ids = String[]
	for line in eachline(coloured_tree)
		if occursin(colour, line)
			push!(coloured_ids, strip(split(line, "["; limit = 2)[1], ['\t', '\'']))
		end
	end

	reader = open(FASTA.Reader, new)
	writer = open(FASTA.Writer, outfile)

	record = FASTA.Record()
	while !eof(reader)
		read!(reader, record)
		id = identifier(record)
		if any(coloured_id -> occursin(coloured_id, id), coloured_ids)
			write(writer, record)
		else
			reader_ = open(FASTA.Reader, original)
			if any(record_ -> id == identifier(record_), reader_)
				write(writer, record)
			end
			close(reader_)
		end
	end

	close(reader)
	close(writer)

	return outfile
end

"""
	rename_ids(outfile_name, organism, transcripts_name, format)
"""
function rename_ids(outfile_name, organism, transcripts_name, format)
	if format == :iterative
		return _rename_ids_iterative(outfile_name, organism, transcripts_name)
	elseif format == :phylo
		return _rename_ids_phylo(outfile_name, organism, transcripts_name)
	else
		return throw(ArgumentError("Input format not defined."))
	end
end

function _rename_ids_iterative(outfile_name, organism, transcripts_name)
	transcripts = open(FASTA.Reader, transcripts_name)
	outfile = open(FASTA.Writer, outfile_name)

	i = 1
	record = FASTA.Record()
	while !eof(transcripts)
		read!(transcripts, record)
		write(outfile, FASTA.Record(string(organism * "_", i), sequence(record)))
		i += 1
	end

	close(transcripts)
	close(outfile)

	return outfile_name
end

function _rename_ids_phylo(outfile_name, organism, transcripts_name)
	transcripts = open(FASTA.Reader, transcripts_name)
	outfile = open(FASTA.Writer, outfile_name)

	record = FASTA.Record()
	while !eof(transcripts)
		read!(transcripts, record)
		write(outfile, FASTA.Record(organism * "@" * split(identifier(record), " "; limit = 2)[1], sequence(record)))
	end

	close(transcripts)
	close(outfile)

	return outfile_name
end

"""
	combine_datasets(out_dir, in_dirs...)
"""
function combine_datasets(out_dir, in_dirs...)

	dict_seq = Dict{String,Vector{String}}() # creates a dictionnary according to the relation: file => list of protein sequence

	record = FASTA.Record()

	for dir in in_dirs # loop through the directories given as arguments
		for f in readdir(dir; join = false) # loop through each file within a directory
			f_out = joinpath(out_dir, f * ".out")
			push!(dict_seq, f_out => String[])
			seq_f_out = dict_seq[f_out]
			reader_f_out = FASTA.Reader(open(f_out, "a")) # creates a new file if it does not already exist
			while !eof(reader_f_out)
				read!(reader_f_out, record)
				push!(seq_f_out, sequence(record))
			end
			close(reader_f_out)
		end
	end

	for dir in in_dirs # loop through the directories given as arguments
		for f in readdir(dir; join = false) # loop through each file within a directory
			reader = open(FASTA.Reader, joinpath(dir, f))
			f_out = joinpath(out_dir, f * ".out")
			seq_f_out = dict_seq[f_out]
			writer_f_out = FASTA.Writer(open(f_out, "a"))
			# append sequence name and sequence protein in the corresponding out file
			while !eof(reader)
				read!(reader, record)
				seq = sequence(record)
				if !(seq in seq_f_out)
					write(writer_f_out, record) # write record
					push!(seq_f_out, seq) # add the new protein sequence to the list corresponding to the in-file
				end
			end
			close(reader)
			close(writer_f_out)
		end
	end

	return collect(keys(dict_seq)) # turns KeySet into Vector
end

	export single_line, lookup, extract_original_and_coloured_taxa,
		rename_ids, combine_datasets

end  # module Transcriptome2
