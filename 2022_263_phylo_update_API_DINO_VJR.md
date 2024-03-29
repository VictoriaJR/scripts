## DOCUMENTATION FOR THE 263 PHYLOGENOMIC DATASET COMBINED APIS AND DINOS - ADDITION OF NEW PARASITES, CORALLICOLIDS, BLASTODINIALES (CURACAO & BC), UNKNOWN COPEPOD INFESTING DINOFLAGELLATE, AMYLOODINIUM (SRA SRR8776921), AND PSAMMOSA PACIFICA + C34

#### Jan 14 2022
last update March 1 2022


### Amyloodinium ocellatum
#### Download amyloodinium transcriptome on rosetta
```bash
/Data/victoria/software/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump SRR8776921 --split-files
```
#### rename fastq files
```bash
mv SRR8776921_1.fastq amyloodinium_ocellatum_R1_001.fastq
mv SRR8776921_2.fastq amyloodinium_ocellatum_R2_001.fastq
```
#### Assemble, clean and translate transcriptome on jezero
```julia
  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
  dir = "/Data/victoria/parasites_proj/Dec_2021_parasites/assemblies/"
for d in readdir(dir; join=true)
        d_name = relpath(d, dir)
        transcriptome_assembly_paired_jezero(
                d, d_name,
                "alveolata_odb10", ["Chordata", "Bacteria", "Malassezia"], "NoPrey", ["fastqc", "cutadapt", "rnaspades", "blastn_megablast", "diamond_blastx", "bowtie2"])
end
```

### Look into older transcriptomes that failed transdecoder:
- take from rosetta and move to jezero
- adding in --no_refine_starts on transdecoder and re-run
- These libraries contain dinoflagellate hits, however the reads are short and full of k-mers. Only transdecoder works when translating without looking for start codon. These libraries may hold interesting information as the cob gene from blastodinium are present. However it is so significantly bad, they will not be used for multi-protein phylogenomics.

```bash
W_blasto1_S7 
W_blast3_S15
```
#### clean data and translate into peptides
```julia
  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
  dir = "/Data/victoria/parasites_proj/Dec_2021_parasites/redo_transdecoder/"
for d in readdir(dir; join=true)
        d_name = relpath(d, dir)
        transcriptome_assembly_paired_jezero(
                d, d_name,
                "alveolata_odb10", ["Chordata", "Bacteria", "Arthropoda"], "NoPrey", ["contamination_removal", "transdecoder"])
end
```

### Combined Psammosa pacifica & Colp-34 transcriptomes

#### 5 psammosa transcriptomes: Psp_2017_2020
- Combine all Psammosa transcriptomes (sequenced using miseq or nextseq)
- these raw data were taken from my scratch folder in soyouz
- I need to rename these files (thats ok because these are not original fastqs)

##### Victoria (2019) 5 individual cells
```bash
mv Ppac1_S1_R1_001.fastq.gz Psp_2017_2020_S1_R1_001.fastq.gz
mv Ppac1_S1_R2_001.fastq.gz Psp_2017_2020_S1_R2_001.fastq.gz
mv Undetermined_S0_R1_001.fastq.gz Psp_2017_2020_S2_R1_001.fastq.gz
mv Undetermined_S0_R2_001.fastq.gz Psp_2017_2020_S2_R2_001.fastq.gz
```
##### Denis (2017) 20 cells
```bash
mv PsammosaS-C_S1_L001_R1_001.fastq.gz Psp_2017_2020_S3_R1_001.fastq.gz
mv PsammosaS-C_S1_L001_R2_001.fastq.gz Psp_2017_2020_S3_R2_001.fastq.gz
```
##### Victoria (2020) 3 transcriptomes
```bash
mv Psp2020_S14_R1_001.cutadapt.fastq Psp_2017_2020_S4_R1_001.cutadapt.fastq
mv Psp2020_S14_R2_001.cutadapt.fastq Psp_2017_2020_S4_R2_001.cutadapt.fastq
mv Psp2020_S18_R1_001.cutadapt.fastq Psp_2017_2020_S5_R1_001.cutadapt.fastq
mv Psp2020_S18_R2_001.cutadapt.fastq Psp_2017_2020_S5_R2_001.cutadapt.fastq
mv Psp2020_S19_R1_001.cutadapt.fastq Psp_2017_2020_S6_R1_001.cutadapt.fastq
mv Psp2020_S19_R2_001.cutadapt.fastq Psp_2017_2020_S6_R2_001.cutadapt.fastq
```

#### Combine all C-34 transcriptomes: C34_2017_2020

##### Denis (2017) whole culture transcriptome
```bash
mv Colp-34RNA_S2_L001_R1_001.fastq.gz C34_2017_2020_S1_R1_001.fastq.gz
mv Colp-34RNA_S2_L001_R2_001.fastq.gz C34_2017_2020_S1_R2_001.fastq.gz
```
##### Victoria (2020) single cell transcriptome
```bash
mv C34-3_S2_R1_001.cutadapt.fastq C34_2017_2020_S2_R1_001.cutadapt.fastq
mv C34-3_S2_R2_001.cutadapt.fastq C34_2017_2020_S2_R2_001.cutadapt.fastq
```

#### Clean up Psammosa data 

##### Assembly of Psammosa pacifica
```julia
  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
  transcriptome_assembly_paired_3_jezero("/Data/victoria/psammosa/Psp_2017_2020/", "Psp_2017_2020", "alveolata_odb10", ["Chordata", "Bacteria"], "Spumella", ["rnaspades", "blastn_megablast", "diamond_blastx", "bowtie2"])
```

#### Assembly of Colp-34
```julia
  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
transcriptome_assembly_paired_2_jezero("/Data/victoria/psammosa/C34_2017_2020/", "C34_2017_2020", "alveolata_odb10", ["Chordata", "Bacteria"], "Procryptobia", ["rnaspades", "blastn_megablast", "diamond_blastx", "bowtie2"])
```


##### Cleaning step
```bash
makeblastdb -in both_Procyrptobia.fa -dbtype nucl -out both_Procyrptobia.DB

blastn -task megablast -query soft_filtered_transcripts_no_Chordata_no_Bacteria.fasta -db /Data/victoria/transcriptomes/Procryptobia/Procryptobia_both.fa.DB -outfmt 6 -num_threads 24 -evalue 1e-25 -max_target_seqs 1 -out c34_clean_vs_Procryptobia.blastnout

#get script from soyouz
cut -f 1 c34_clean_vs_Procryptobia.blastnout > Kinetoplastid_contigs.list
996 Kinetoplastid_contigs.list

perl /Data/victoria/scripts/lookup_reverse.pl soft_filtered_transcripts_no_Chordata_no_Bacteria.fasta Kinetoplastid_contigs.list

mv lookup_out.fasta Trinity_noBac_noPrey.fasta
```


#### Collect peptides that will be added this round in one folder and have them renamed to the appropriate name:
my parasites:
GB_who_decH9
Cur_GB1_H2803
Cur_GB2_H2806
Cur_GB3_H2810
Cur_GB4_H2813
Cur_GB5_H2815
GB_gonspo_decH7
GB_spo_decH7
mysid_parasite_3273_H1
mysid_parasite_3537_H2
copepod_infesting_alv_3147

downloaded parasites (cleaned, removed malassezia):
amyloodinium_ocellatum

Corallicolids (MiSeq and NextSeq combined):
GZC_10102021
GZC_10142021
M6

Final combined Psammosa transcriptomes:
Psp_2017_2020
C34_2017_2020


```julia
include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "GB_who_decH9")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed GB_who_decH9.fasta

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "Cur_GB1_H2803")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed Cur_GB1_H2803.pep.fasta
mv Cur_GB1_H2803.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "Cur_GB2_H2806")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed Cur_GB2_H2806.pep.fasta
mv Cur_GB2_H2806.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "Cur_GB3_H2810")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed Cur_GB3_H2810.pep.fasta
cp Cur_GB3_H2810.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "Cur_GB4_H2813")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed Cur_GB4_H2813.pep.fasta
cp Cur_GB4_H2813.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "Cur_GB5_H2815")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed Cur_GB5_H2815.pep.fasta
cp Cur_GB5_H2815.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "GB_gonspo_decH7")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed GB_gonspo_decH7.pep.fasta
cp GB_gonspo_decH7.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Malassezia.fasta.transdecoder.pep", "A_ocellatum")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Malassezia.fasta.transdecoder.pep.renamed A_ocellatum.pep.fasta
cp A_ocellatum.pep.fasta ../../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "GB_spo_decH7")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed GB_spo_decH7.pep.fasta
cp GB_spo_decH7.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "mysid_parasite_3273_H1")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed mysid_parasite_3273_H1.pep.fasta
cp mysid_parasite_3273_H1.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "mysid_parasite_3537_H2")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed mysid_parasite_3537_H2.pep.fasta
cp mysid_parasite_3537_H2.pep.fasta ../transcriptomes_to_add


include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep", "copepod_infesting_alv_3147")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Arthropoda.fasta.transdecoder.pep.renamed copepod_infesting_alv_3147.pep.fasta
cp copepod_infesting_alv_3147.pep.fasta ../transcriptomes_to_add

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Spumella.fasta.transdecoder.pep", "Psp_2017_2020")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Spumella.fasta.transdecoder.pep.renamed Psp_2017_2020.pep.fasta
cp Psp_2017_2020.pep.fasta /Data/victoria/parasites_proj/Dec_2021_parasites/transcriptomes_to_add/

include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Procryptobia.fasta.transdecoder.pep", "C34_2017_2020")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria_no_Procryptobia.fasta.transdecoder.pep.renamed C34_2017_2020.pep.fasta
cp C34_2017_2020.pep.fasta /Data/victoria/parasites_proj/Dec_2021_parasites/transcriptomes_to_add/

```

#### All cleaned, renamed peptides:
location on jezero
moved from `/Data/victoria/parasites_proj/Dec_2021_parasites/transcriptomes_to_add/` to `/Data/victoria/263_api_dino_2022/transcriptomes_to_add`

#### Make blast databases for each new transcriptome that's being added:

    for FA in *.fasta ; do makeblastdb -in $FA -dbtype prot -out $FA.DB ; done
    
##### Make a list of all the blast databases for future programs to reference:

	  for fasta in *.fasta ; do echo $fasta'.DB' >> blastdb_list ; done  
    
##### Setup work space:
- Copy the most up-to-date 263 gene dataset--> 263genes_jan17_comb_VJR_IN into my working directory
- Copy the running_blast.pl script to the 263 genes directory and alter the path to where the blastdb_list is.

##### In a screen and inside the 263_genes dir, run the running_blast.pl script:
- This will use the original fasta files and the databases (new transcriptomes) to find genes
- This will create (number of new taxa) x263 blast outputs (GENE*DATABASE.blastout)
	
		for i in *.fasta; do perl running_blast.pl $i /Data/victoria/263_api_dino_2022/transcriptomes_to_add/blastdb_list ; done

###### Parse the BLAST outputs (e-value = 1e-20, query coverage = 50%) and compile all parsed output into one file containing only non redundant sequences: (This new list is called `seq_list_toadd.txt`)
- `batch_parsing_blast_v2.pl` will generate a file of the blastout with _parsed.txt
- `get_list_toadd.pl` will print something like `Gene name: ABCE.fasta` for all the genes and output `seq_list_toadd.txt`

		perl /Data/victoria/263_api_dino_2022/batch_parsing_blast_v2.pl *.blastout
		
		perl /Data/victoria/263_api_dino_2022/get_list_toadd.pl *_parsed.txt

- We now need to trim the identified genes from our added dataset so that they don't have extensions accidentally added in the alignment.
- To do this we blast the sequences against the swissprot data base and trim any parts of the sequences that are not well aligned with the well curated proteins.


###### First, make a master file of added datasets:
	
	cat *.fasta > added_datasets.fa

- Before running the gene finding script, make sure to move the original 263_genes fasta files (I moved them up one directory) to a new folder or else they will be overwritten/confused with new ones. You need to have the added_datasets.fa and seq_list_toadd.txt file in the same directory.
- by now your added_datasets.fa. will be indexed
- this will also create the gene.fasta output with only your sequences from the databases in it

		perl /Data/victoria/263_api_dino_2022/lookup_single_genes.pl added_datasets.fa

##### Create a Blast database using Swissprot proteins to blast cleaned ORFs against.
- First, make a directory:

		mkdir blast_db

###### Inside this directory, download the most recent version of the Swissprot database:

		wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
		gunzip uniprot_sprot.fasta.gz

##### Make a blast database with the Swissprot fasta:

	makeblastdb -in /Data/victoria/263_api_dino_2022/blast_db/uniprot_sprot.fasta -dbtype prot -out swissprot_feb2022

- This script requires the seq_list_toadd.txt file and the added_datasets.fa file - it will find the sequences it needs to get in the seq_list_toadd.txt file and then find them in the added_datasets.fa file. Fasta files will then be exported for each of the 263 genes.		
- Blast the extracted sequences against the swissprot database. Put swissprot fasta in a different directory so it doesn't get caught up in the next steps.

###### In a screen, do a BLAST search:
- this will output a blastout of your gene.fasta file

		for i in *fasta; do blastp -query $i -db /Data/victoria/263_api_dino_2022/blast_db/swissprot_feb2022 -evalue 1e-5 -num_alignments 5 -out $i.sprot.blastout -num_threads 25 ; done

- This may not produce 263 blastout files; that's ok.
- The blast search will identify sequences that should be trimmed.

###### Use the rm_craps.pl script to trim those sequences and output a new set of trimmed fasta files ending with '.nocrap':

	for i in *.fasta ; do perl ../rm_craps.pl $i ; done

###### Move the nocrap files one directory up (or wherever you moved the original fastas to get them out of the way earlier):
```bash
mv *nocrap ../original_fastas
cd ../original_fastas
for i in *.fasta ; do cat $i $i.nocrap > $i.new ; done
```
- this will result in some "no such file or directory" alerts if there were less than 263 blastout files. again, this is fine.)
- This will output new fasta files for each gene called *.new
- **Note: no files for genes: COQ4mito, CTU1, MMAAmito, PLS3**

#### SINGLE GENE TREES CONSTRUCTION
In a screen, make and trim alignments for all genes with new, cleaned sequences added using Linsi and Trimal:

	for i in *.new ; do linsi --thread 24 $i > $i.linsi ; done
	for i in *.linsi ; do trimal -in $i -out $i.trimal -gt 0.8 ; done

##### for EFTUD1 (where alignments were so trimmed, that only 3 aa length sequences remained)
	
	for i in *.linsi ; do trimal -in $i -out $i.trimal -gt 0.7 ; done

##### Remove any headers that lack sequences (may not exist):

	for i in *.trimal ; do awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' $i > $i.clean ; done

```julia
function dir_remove_empty_seq(dir)
	for f in readdir(dir)
		file_remove_empty_seq(f)
	end
end

function file_remove_empty_seq(file)
	lines = split(read(file, String), ">"; keepempty = false)
	clean = file * ".no_empty_seq"
	open(clean, "w") do io
		for line in lines
			split_line = split(line, "\n"; limit = 2)
			header = split_line[1]
			seq = split_line[2]
			if !isempty(replace(replace(seq, "\n" => ""), "-" => ""))
				write(io, ">" * header * "\n" * seq)
			end
		end
	end
	return clean
end

dir_remove_empty_seq("/Data/victoria/263_api_dino_2022/original_263/clean/")

```
```bash
rename 's/.no_empty_seq//' *seq
```

##### Run the trees on compute canada
- Upload *.trimal.clean files to ssh vjackor@graham.computecanada.ca to scratch folder # done on march 1 2022

graham:

	python /home/vjackor/scripts/RAxML_263_April2020_2.py '*.trimal.clean'
	for SUBFILE in *.sh ; do sbatch $SUBFILE ; done
		
cedar:	

	make -f Makefile.SSE3.PTHREADS.gcc
	python /home/vjackor/scripts/single_gene_tree_sub_scripts_RAxML_SSE3_PTHREADS.py '*.trimal.clean'
	for SUBFILE in *.sh ; do sbatch $SUBFILE ; done

##### Check for failed computations, move to re-runs folder
- Make a new output folder called `rerun` composed of all the files that did not generate final RAxML trees.

```julia
function check_for_reruns(dir, re_run_dir)
    prefix_w1 = "RAxML_info."
    prefix_w2 = "RAxML_bootstrap."
    prefix_w3 = "RAxML_bipartitionsBranchLabels."
    prefix_w4 = "RAxML_bestTree."
    prefix_w5 = "sub_"
    suffix_w1 = ".out"
    suffix_w2 = ".reduced"
    prefix_c = "RAxML_bipartitions."
    for f in readdir(dir)
        if isfile(f) && !startswith(f, prefix_w1) && !startswith(f, prefix_w2) && !startswith(f, prefix_w3) && !startswith(f, prefix_c) &&
                !startswith(f, prefix_w4) && !startswith(f, prefix_w5) && !endswith(f, suffix_w1) && !endswith(f, suffix_w2)
            if !any(d -> startswith(d, prefix_c * f), readdir(dir))
                cp(f, re_run_dir * "/" * f; force=false)
            end
        end
    end
end
check_for_reruns("/home/vjackor/scratch/263_genes_Feb2022", "/home/vjackor/scratch/263_genes_Feb2022/rerun")
```
    


##### CLEAN TREES

Make a new directory called coloured_trees
- Copy files ending in .new from 1_pre_tree_pipeline: these files contain
untrimmed protein sequences from the orginial files + the new dino sequences
- Make sure original fasta files (prior to trimming) and coloured trees are in one folder
- Make single line format files from the .new files (will produce .sl files)


`for file in *.new; do python /Data/victoria/scripts/single_line_fasta.py $file; done`
`rename 's/_coloured.tre/.fasta.new.linsi.trimal.clean.treefile_coloured/' *.tre`

- Use the julia script to keep colored sequences from .new files using the corresponding colored trees:

Move to rosetta

	/Data/victoria/software/julia-1.5.3/bin/julia
		include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
		dir_extract_original_and_coloured_taxa("/Data/victoria/263_june2022/coloured/", "ff0000")

**here i downloaded the files and added in the merged sequences and removed a few sequences from previous analyses**
**do do in the future, re-run Gpn1 - might have a reverse sequence to fix & figure out Rps6 there are two Psp sequences added, so i renamed the orthologs as name_ortho ** look into



Move the cleaned fasta files (.cleaned) to a new directory (renamed to 263_genes_11June2022_VJR). Enter this new directory and remove all empty files:

	find . -size 0 -delete

#####  Files in 263_genes_11June2022_VJR  will be the new 263_genes files that you will start this whole giant process with if you want to add subsequent species to protein phylogenies

Then make sure no contaminations or paralogs were missed by counting the number of sequences for the newly added taxa in each fasta (they should all be 1 or 0).

download all the *.cleaned files to make it easier to delete the more than 1 OTU; may count the # times the organism is present in the dataset
- downloaded and re-uploaded modified fastas here

			grep -c 'file' *.cleaned
notes for this round:
grep -c 'GB_who_decH9' *.cleaned
grep -c 'Cur_GB1_H2803' *.cleaned
grep -c 'Cur_GB2_H2806' *.cleaned
grep -c 'Cur_GB3_H2810' *.cleaned
grep -c 'Cur_GB4_H2813' *.cleaned
grep -c 'Cur_GB5_H2815' *.cleaned
grep -c 'GB_gonspo_decH7' *.cleaned #dont add to tree
grep -c 'GB_spo_decH7' *.cleaned #dont add to tree
grep -c 'mysid_parasite_3273_H1' *.cleaned #dont add to tree
grep -c 'mysid_parasite_3537_H2' *.cleaned #dont add to tree
grep -c 'copepod_infesting_alv_3147' *.cleaned #did i even add this one?**
grep -c 'A_ocellatum' *.cleaned 
grep -c 'GZC_10102021' *.cleaned 
grep -c 'GZC_10142021' *.cleaned 
grep -c 'M6' *.cleaned 
grep -c 'Psp_2017_2020' *.cleaned 
grep -c 'C34_2017_2020' *.cleaned 


##### Rename cleaned fasta files
If you want to rename -->
**did not do it for the this round of analyses**
rename: thecate_dino_exH11 to unknown_dino_exH11
			unknown_112120-11 to unknown_dino_112120-11

		for file in *.cleaned ; do sed 's/thecate_dino_exH11/unknown_dino_exH11/g' $file > $file.fixed ; done
		rename 's/.cleaned.fixed/.cleaned/' *.fixed

		 Rename these cleaned_fastas

		grep -c 'unknown_dino_112120-11' *.cleaned

		rename 's/.fasta.new.sl.cleaned/.fasta/' *.cleaned



#### Make your tree
##### tree 1 on soyouz

##### Scafos : Concatentate 263 genes
Soyouz: /opt/scafos/scafos
jezero: /usr/local/bin/scafos (1.25)

##### Aligned all clean alignments with linsi (can do overnight in a screen)

	for fasta in *.fasta; do linsi --thread 24 $fasta > $fasta.linsi; done


##### Trim using trimal (jezero)

    for i in *.linsi ; do trimal -in $i -out $i.trimal -gt 0.7 ; done


##### Move all trimmed files to a new folder named "scafos"

	mv *.trimal /Data/victoria/scafos/

##### Change the names to GENE.fasta

	rename 's/.fasta.linsi.trimal/.fasta/' *.trimal

##### Open scafos jezero 
	
    ssh -Y victoria@soyouz.zoology.ubc.ca
    /usr/local/bin/scafos

##### Troubleshooting (must fix before first step):
- To fix bad characters, remove them and replace them with a space:
- Get rid of original .fastas and remove the ".clean" from the end of the new files:

	for file in *.fasta ; do sed 's/|/ /g' $file > $file.clean ; done
	rename 's/.fasta.clean/.fasta/' *.clean

-Other issues: files with empty sequences from trimming **just ran the script written in julia from above**


##### SPECIES PRESENCE
- Download the species_presence-freq.otu file and delete the taxa that you do not want:
- Generally, remove all low % data species. The exception is a species that's really important for your species (in this case, I kept Oxyrrhis marina even though it was 32%)
- Figure out what your outgroup is (ciliates for me in this case) and get rid of everything more distantly related than that
- Keep the relevant species that have the highest data

#### Stats of coverage across 263 genes

##### FILE SELECTION
- Run 'file selection' in SCAFOS giving it the new species presence list (after removing taxa) as the OTU file (call the output directory 'file_selection').
- Before doing this, remove random "bak" folder from your 'scafos' directory if there is one.
- Files with ONLY the selected species will be in the output directory.

***did not select the less than 40% OTUs for my first trees- i used all 263 gene trees***
##### DATASET ASSEMBLING 
- Run 'dataset assembling' on the directory containing the fasta with the selected species (it also needs the OTU file) and select 'longer sequence' in the selection criteria options
- In the output of dataset assembly download the .stat file. Open it in excel. Flip (or "transpose") the species and gene table so that genes are the rows. Sort the rows by missing data or otus. Highlight the genes above a certain threshold (ex. 40% missing data) and copy them into a file (ex. cat > genes_over_40.txt).
- Genes over 40 actually means values less than 40 (since it is ' % missing OTUs')
- Move the good genes to a new directory (make a directory called 'good_genes' first):

	for i in `less good_genes_may.txt` ; do cp /Data/victoria/file_selection/$i /Data/victoria/good_genes/ ; done

	for i in `less good_genes.txt` ; do cp /Data/victoria/scafos/file_selection/$i /Data/victoria/scafos/good_genes/ ; done

##### SPECIES PRESENCE (Concatenate)
- Re-run 'species presence' on this new directory containing the good genes. When it asks if you want to concatenate, say yes. Name output directory 'concatenation_date'.

##### MAKE TREE

	FastTree cat_april162021.fasta > cat_april162021.fasta.fasttree

	for file in *.fasta ; do iqtree-omp -m TEST -bb 1000 -nt 20 -s $file ; done


**do this for my next tree**
Clean up organism headers before running tree 
**remove dino_**
	
    for file in *.fasta ; do sed 's/dino_//g' $file > $file.1 ; done
	for file in *.1 ; do sed 's/Dinophyta_Dinophyceae-//g' $file > $file.2 ; done
	for file in *.2 ; do sed 's/Pyrrophycophyta_Dinophyceae-//g' $file > $file.3 ; done
	for file in *.3 ; do sed 's/Dinoflagellate-//g' $file > $file.4 ; done
	for file in *.4 ; do sed 's/A_ocellatum/Amyloodinium_ocellatum/g' $file > $file.5 ; done
	for file in *.5 ; do sed 's/Sarai//g' $file > $file.6 ; done

	for file in *.6 ; do sed 's/C34_2017_2020/Colp-34_2017_2020/g' $file > $file.7 ; done
	for file in *.7 ; do sed 's/Colp37_culture/Colp37_RNA/g' $file > $file.8 ; done
	for file in *.8 ; do sed 's/Colp37_singlecell/Colp37_sc/g' $file > $file.9 ; done
	for file in *.9 ; do sed 's/TGD_/TGD/g' $file > $file.10 ; done
	for file in *.10 ; do sed 's/Polykrikos_lebourae_GG/Polykrikos_lebourae/g' $file > $file.11 ; done
	for file in *.11 ; do sed 's/Psp_2017_2020/Psammosa_pacifica_2017_2020/g' $file > $file.12 ; done
	for file in *.12 ; do sed 's/Colp-34_2017_2020/Colp34_2017_2020/g' $file > $file.13 ; done





### other phylogenetic estimations
	
    /home/vjackor/scripts/Phylobayes_script.sh
    sbatch /home/vjackor/scripts/Phylobayes_script.sh  tree4_concat.fasta

**remove >**
	for file in *.3 ; do sed 's/>//g' $file > $file.4 ; done

**place seq beside name -- tried to remove carriage return and spaces but nothing changed**
	for file in *.4 ; do sed 's/ //g' $file > $file.5 ; done








### Combine datasets on Rosetta
#### combine Ina's most current version of the 263 gene tree on google drive:
#### 263genes_nov012021_IN
#### 263_genes_24april2021_VJR

#### check each folder has 263 genes using `ls | wc -l`

#### script:

location on rosetta: `/Data/victoria/263_db/263_genes_comb_VJR_IN`


#### Run combination script that adds both header and sequence **
    Example on how to use the new combine datasets julia script (rosetta)

      /Data/victoria/software/julia-1.5.3/bin/julia
      include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
      combine_datasets("/Data/victoria/263_db/263_genes_comb_VJR_IN/263_combined_out_false/", "/Data/victoria/263_db/263_genes_comb_VJR_IN/263_trees_28april2021_VJR", "/Data/victoria/263_db/263_genes_comb_VJR_IN/263genes_nov012021_IN"; check_id = true

#### notes from previous run:
- This script will first use my dataset to make a list of all the headers and sequences.
      then it will add any header or sequence that does not match the one my
      dataset. Sometimes you will get copies of the same sequence, beacuse of slight changes
      in the header (font or format etc). My method to remove copies is
      - maintain order of your fasta addition (i.e. remove the newest duplicate)
      - you can see where your last added transcriptomes are and remove all
      the unwanted duplicates following yours. Do not change order of sequences.

***under construction**
#### cleaned up version is on jezero at : 

I will collect all of Ina and My new sequences and add them to Liz's 263 curated database
- Liz has removed many old sequences, and fixed up paralogies :) 
- Use scafos to pull out my sequences and add to her dataset

#### to do after this phylo update
Another thing I want to do: re-add transcriptomes using combined datasets (parasites from the same host infection)


