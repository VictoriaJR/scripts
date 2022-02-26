## DOCUMENTATION FOR THE 263 PHYLOGENOMIC DATASET COMBINED APIS AND DINOS - ADDITION OF NEW PARASITES, CORALLICOLIDS, BLASTODINIALES (CURACAO & BC), UNKNOWN COPEPOD INFESTING DINOFLAGELLATE, AMYLOODINIUM (SRA SRR8776921), AND PSAMMOSA PACIFICA + C34

#### Jan 14 2022
last update Feb 24 2022


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




old notes below, edit as you go 
______________


###### First, make a master file of added datasets:
	
	cat *.fasta > added_datasets.fa

- Before running the gene finding script, make sure to move the original 263_genes fasta files (I moved them up one directory) to a new folder or else they will be overwritten/confused with new ones. You need to have the added_datasets.fa and seq_list_toadd.txt file in the same directory.
- by now your added_datasets.fa. will be indexed
- this will also create the gene.fasta output with only your sequences from the databases in it

		perl /Data/victoria/scripts/lookup_single_genes.pl added_datasets.fa


##### Create a Blast database using Swissprot proteins to blast cleaned ORFs against.
- First, make a directory:

		mkdir blast_db

###### Inside this directory, download the most recent version of the Swissprot database:

		wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
		gunzip uniprot_sprot.fasta.gz

##### Make a blast database with the Swissprot fasta:

	makeblastdb -in /Data/victoria/april2021_263_genes/blast_db/uniprot_sprot.fasta -dbtype prot -out swissprot_april2021

- This script requires the seq_list_toadd.txt file and the added_datasets.fa file - it will find the sequences it needs to get in the seq_list_toadd.txt file and then find them in the added_datasets.fa file. Fasta files will then be exported for each of the 263 genes.		
- Blast the extracted sequences against the swissprot database. Put swissprot fasta in a different directory so it doesn't get caught up in the next steps.

###### In a screen, do a BLAST search:
- this will output a blastout of your gene.fasta file

		for i in *fasta; do blastp -query $i -db /Data/victoria/april2021_263_genes/blast_db/swissprot_april2021 -evalue 1e-5 -num_alignments 5 -out $i.sprot.blastout -num_threads 25 ; done

- This may not produce 263 blastout files; that's ok.
- The blast search will identify sequences that should be trimmed.

###### Use the rm_craps.pl script to trim those sequences and output a new set of trimmed fasta files ending with '.nocrap':

	for i in *.fasta ; do perl /Data/victoria/scripts/rm_craps.pl $i ; done

###### Move the nocrap files one directory up (or wherever you moved the original fastas to get them out of the way earlier):
```bash
mv *nocrap ../original_fastas
cd ../original_fastas
for i in *.fasta ; do cat $i $i.nocrap > $i.new ; done
```

- this will result in some "no such file or directory" alerts if there were less than 263 blastout files. again, this is fine.)
- This will output new fasta files for each gene called *.new




*tree construction*
SINGLE GENE TREES CONSTRUCTION

In a screen, make and trim alignments for all genes with new, cleaned sequences added using Linsi and Trimal:

	for i in *.new ; do linsi --thread 24 $i > $i.linsi ; done
###rosetta
	for i in *.linsi ; do /opt/trimAl/source/trimal -in $i -out $i.trimal -gt 0.8 ; done
###soyouz
  for i in *.linsi ; do trimal -in $i -out $i.trimal -gt 0.8 ; done

  for EFTUD1 (where alignments were so trimmed, that only 3 aa length sequences remained)
	for i in *.linsi ; do trimal -in $i -out $i.trimal -gt 0.7 ; done


Remove any headers that lack sequences (may not exist):

	for i in *.trimal ; do awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' $i > $i.clean ; done

	Run the trees: SWTICH TO COMPUTE CANADA now!

	Upload *.trimal.clean files to ssh vjackor@graham.computecanada.ca to scratch folder #DONE ON MARCH 29 2021

		python /home/vjackor/scripts/RAxML_263_April2020_2.py '*.trimal.clean'
		for SUBFILE in *.sh ; do sbatch $SUBFILE ; done


	*** some files were taking a long time to run on graham, so i made a list what
	was remaining in the squeue -u vjackor que, and used that to pull the *.clean
	files from soyouz and put in their own folder and put on cedar
	```julia
	/Data/victoria/software/julia-1.5.3/bin/julia
	fasta_list_file = "/Data/victoria/march2021/263_genes_27march2021/fasta_list_file.txt"
	for line in eachline(fasta_list_file)
	    line_name = line * ".fasta.new.linsi.trimal.clean"
	        run(`cp $line_name /Data/victoria/march2021/263_genes_27march2021/to_download_on_cedar`)
	end
	```
	remaining on cedar:

	#run this first in the folder of the RAxML program
	#compile the type of raxml software I want to use on my local system
	make -f Makefile.SSE3.PTHREADS.gcc

	python /home/vjackor/scripts/single_gene_tree_sub_scripts_RAxML_SSE3_PTHREADS.py '*.trimal.clean'
	for SUBFILE in *.sh ; do sbatch $SUBFILE ; done

### Here I used a newly written julia script to identify trimmed files on cc that
### were not run beacuse of empty sequences: re-run_fastas_for_263.jl
rename 's/.no_empty_seq.no_empty_seq.no_empty_seq.no_empty_seq.no_empty_seq//' *seq

		#### CLEAN TREES

		Make a new directory called coloured_trees
		- Copy files ending in .new from 1_pre_tree_pipeline: these files contain
			untrimmed protein sequences from the orginial files + the new dino sequences
		- Make sure original fasta files (prior to trimming) and coloured trees are in one folder
		- Make single line format files from the .new files (will produce .sl files)


			for file in *.new; do python /Data/victoria/scripts/single_line_fasta.py $file; done
			rename 's/_coloured.tre/.fasta.new.linsi.trimal.clean.treefile_coloured/' *.tre

		- Use the julia script to keep colored sequences from .new files using the corresponding colored trees:
		- Move to rosetta
		/Data/victoria/software/julia-1.5.3/bin/julia
		include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
		dir_extract_original_and_coloured_taxa("/Data/victoria/temp/", "ff0000")

		Move the cleaned fasta files (.cleaned) to a new directory (ex. 263_genes_24april2021_VJR). Enter this new directory and remove all empty files:

			find . -size 0 -delete

		#####  Files in 263_genes_24april2021_VJR  will be the new 263_genes files that you will start this whole giant process with if you want to add subsequent species to protein phylogenies

		Then make sure no contaminations or paralogs were missed by counting the number of sequences for the newly added taxa in each fasta (they should all be 1 or 0).

		download all the *.cleaned files to make it easier to delete the more than 1 OTU

			grep -c 'file' *.cleaned


		may count the # times the organism is present in the dataset

			grep -c 'Name' *.cleaned
			grep -c 'Name' *.out

		#Rename cleaned fasta files

rename: thecate_dino_exH11 to unknown_dino_exH11
			unknown_112120-11 to unknown_dino_112120-11

		for file in *.cleaned ; do sed 's/thecate_dino_exH11/unknown_dino_exH11/g' $file > $file.fixed ; done
		rename 's/.cleaned.fixed/.cleaned/' *.fixed

		for file in *.cleaned ; do sed 's/unknown_112120-11/unknown_dino_112120-11/g' $file > $file.fixed ; done
		rename 's/.cleaned.fixed/.cleaned/' *.fixed

		 Rename these cleaned_fastas

		grep -c 'unknown_dino_112120-11' *.cleaned
		grep -c 'unknown_dino_exH11' *.cleaned



		grep -c 'yellow_egg' *.cleaned
		grep -c 'thecate_dino_20201206-3' *.cleaned
		grep -c 'Unknown_112120-10' *.cleaned
		grep -c 'white_blasto_sporo_exH8' *.cleaned
		grep -c 'oodinium_1031' *.cleaned
		grep -c 'white_blast2_who_exH6' *.cleaned
		grep -c 'white_blasto_exH10' *.cleaned
		grep -c 'white_blast_gon_exH9' *.cleaned
		grep -c 'white_blasto_spo_exH9' *.cleaned

		rename 's/.fasta.new.sl.cleaned/.fasta/' *.cleaned



		##### Scafos : Concatentate 263 genes
		(on Soyouz: /opt/scafos/scafos)


		    screen -r 263

		Aligned all clean alignments with linsi (can do overnight in a screen)

		    for fasta in *.fasta; do linsi --thread 24 $fasta > $fasta.linsi; done


		## do next
		Trim using trimal

		    for linsi in *.linsi; do /opt/trimAl/source/trimal -in $linsi -out $linsi.trimal -gt 0.8; done
		#soyouz
			  for linsi in *.linsi; do trimal -in $linsi -out $linsi.trimal -gt 0.7; done

		Move all trimmed files to a new folder named "scafos"

		    mkdir scafos closest to home directory
		    mv *.trimal /Data/victoria/scafos/

		Change the names to GENE.fasta

		    rename 's/.fasta.linsi.trimal/.fasta/' *.trimal

		##### Open scafos on soyouz: (log in with -Y)

		xs
		LOGIN:

		    ssh -Y victoria@soyouz.zoology.ubc.ca
		    /opt/scafos/scafos

		##### Troubleshooting (must fix before first step):

		- To fix bad characters, remove them and replace them with a space:

		      for file in *.fasta ; do sed 's/|/ /g' $file > $file.clean ; done


		- Get rid of original .fastas and remove the ".clean" from the end of the new files:

		      rename 's/.fasta.clean/.fasta/' *.clean

		- Other issues: files with empty sequences from trimming

		  - NAA15
		  - EFTUD1
		  - BSM1
		  - PRPF8
		  - RPS19 (double abedinium seq: kept liz's)


		##### SPECIES PRESENCE
		Download the species_presence-freq.otu file and delete the taxa that you do not want:
		-Generally, remove all low % data species. The exception is a species that's really important for your species (in this case, I kept Oxyrrhis marina even though it was 32%)
		-Figure out what your outgroup is (ciliates for me in this case) and get rid of everything more distantly related than that
		-Keep the relevant species that have the highest data

		#### Stats of coverage across 263 genes
		oodinium : oodinium (42%)
		GB_spo_exH1 : GB_spo_exH1 (47%)
		GB_gon_exH4 : GB_gon_exH4 (41%)
		GB_spo_exH4 : GB_spo_exH4 (30%)
		GB_mas_exH1 : GB_mas_exH1 (29%)
		WBM_BH3 : WBM_BH3 (26%)
		WB_gon_exH3 : WB_gon_exH3 (15%)
		WB_spo_exH3 : WB_spo_exH3 (14%)
		GB_spo_exH2 : GB_spo_exH2 (6%)

		##### FILE SELECTION

		2)Run 'file selection' in SCAFOS giving it the new species presence list (after removing taxa) as the OTU file (call the output directory 'file_selection').

		Before doing this, remove random "bak" folder from your 'scafos' directory if there is one.
		Files with ONLY the selected species will be in the output directory.


		##### DATASET ASSEMBLING

		3)Run 'dataset assembling' on the directory containing the fasta with the selected species (it also needs the OTU file) and select 'longer sequence' in the selection criteria options

		In the output of dataset assembly download the .stat file. Open it in excel. Flip (or "transpose") the species and gene table so that genes are the rows. Sort the rows by missing data or otus. Highlight the genes above a certain threshold (ex. 40% missing data) and copy them into a file (ex. cat > genes_over_40.txt).
		Genes over 40 actually means values less than 40 (since it is ' % missing OTUs')
		Move the good genes to a new directory (make a directory called 'good_genes' first):

		for i in `less good_genes_may.txt` ; do cp /Data/victoria/file_selection/$i /Data/victoria/good_genes/ ; done

		for i in `less good_genes.txt` ; do cp /Data/victoria/scafos/file_selection/$i /Data/victoria/scafos/good_genes/ ; done
		/Data/victoria/scafos/dataset_assembling

		##### SPECIES PRESENCE (Concatenate)

		Re-run 'species presence' on this new directory containing the good genes. When it asks if you want to concatenate, say yes. Name output directory 'concatenation_date'.

		##### MAKE TREE

		    FastTree cat_april162021.fasta > cat_april162021.fasta.fasttree

		    for file in *.fasta ; do iqtree-omp -m TEST -bb 1000 -nt 5 -s $file ; done

		    ended up running RAxML for trees on compute canada
		    /home/vjackor/scripts/Phylobayes_script.sh

		    sbatch /home/vjackor/scripts/Phylobayes_script.sh  tree4_concat.fasta


		transform data for phylobayes
		#taxa #sites
		name1 seq1.....
		name2 seq2.....

		61
		currently looks like
		>name1
		seq1

		add a ! at the end of each name that starts with a header

		remove dino_
		for file in *.fasta ; do sed 's/dino_//g' $file > $file.fixed ; done
		for file in *.fixed ; do sed 's/Dinophyta_Dinophyceae-//g' $file > $file.2 ; done
		for file in *.2 ; do sed 's/Pyrrophycophyta_Dinophyceae-//g' $file > $file.3 ; done

		remove >
		for file in *.3 ; do sed 's/>//g' $file > $file.4 ; done

		place seq beside name -- tried to remove carriage return and spaces but nothing changed
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


