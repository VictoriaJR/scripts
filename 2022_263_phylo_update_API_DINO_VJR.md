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
##### SCREEN screen -r neq TO RERUN THIS

```julia
  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
  dir = "/Data/victoria/parasites_proj/Dec_2021_parasites/assemblies/"
for d in readdir(dir; join=true)
        d_name = relpath(d, dir)
        transcriptome_assembly_paired_jezero(
                d, d_name,
                "alveolata_odb10", ["Chordata", "Bacteria", "Malassezia"], "NoPrey", ["contamination_removal", "transdecoder"])
end
```

#### Look into older transcriptomes that failed transdecoder:
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

### Collect Psammosa pacifica + Colp-34? transcriptomes (document whats added and what has not yet been added, should I combine these?):

##### 5 psammosa transcriptomes: Psp_2017_2020
- Combine all Psammosa transcriptomes (sequenced using miseq or nextseq)
- these raw data were taken from my scratch folder in soyouz
- I need to rename these files (thats ok because these are not original fastqs)

#### Victoria (2019) 5 individual cells
mv Ppac1_S1_R1_001.fastq.gz Psp_2017_2020_S1_R1_001.fastq.gz
mv Ppac1_S1_R2_001.fastq.gz Psp_2017_2020_S1_R2_001.fastq.gz
mv Undetermined_S0_R1_001.fastq.gz Psp_2017_2020_S2_R1_001.fastq.gz
mv Undetermined_S0_R2_001.fastq.gz Psp_2017_2020_S2_R2_001.fastq.gz

#### Denis (2017) 20 cells
mv PsammosaS-C_S1_L001_R1_001.fastq.gz Psp_2017_2020_S3_R1_001.fastq.gz
mv PsammosaS-C_S1_L001_R2_001.fastq.gz Psp_2017_2020_S3_R2_001.fastq.gz

#### Victoria (2020) 3 transcriptomes
mv Psp2020_S14_R1_001.cutadapt.fastq Psp_2017_2020_S4_R1_001.cutadapt.fastq
mv Psp2020_S14_R2_001.cutadapt.fastq Psp_2017_2020_S4_R2_001.cutadapt.fastq
mv Psp2020_S18_R1_001.cutadapt.fastq Psp_2017_2020_S5_R1_001.cutadapt.fastq
mv Psp2020_S18_R2_001.cutadapt.fastq Psp_2017_2020_S5_R2_001.cutadapt.fastq
mv Psp2020_S19_R1_001.cutadapt.fastq Psp_2017_2020_S6_R1_001.cutadapt.fastq
mv Psp2020_S19_R2_001.cutadapt.fastq Psp_2017_2020_S6_R2_001.cutadapt.fastq


#### Combine two C-34 transcriptomes: C34_2017_2020

#### Denis (2017) whole culture transcriptome
mv Colp-34RNA_S2_L001_R1_001.fastq.gz C34_2017_2020_S1_R1_001.fastq.gz
mv Colp-34RNA_S2_L001_R2_001.fastq.gz C34_2017_2020_S1_R2_001.fastq.gz

#### Victoria (2020) single cell transcriptome
mv C34-3_S2_R1_001.cutadapt.fastq C34_2017_2020_S2_R1_001.cutadapt.fastq
mv C34-3_S2_R2_001.cutadapt.fastq C34_2017_2020_S2_R2_001.cutadapt.fastq


```julia
## Psammosa pacifica ** Do clean up step by hand

  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
  transcriptome_assembly_paired_3_jezero("/Data/victoria/psammosa/Psp_2017_2020/", "Psp_2017_2020", "alveolata_odb10", ["Chordata", "Bacteria"], "Spumella", ["rnaspades", "blastn_megablast", "diamond_blastx", "bowtie2"])


### Colp-34 ** Do clean up step by hand

  include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
transcriptome_assembly_paired_2_jezero("/Data/victoria/psammosa/C34_2017_2020/", "C34_2017_2020", "alveolata_odb10", ["Chordata", "Bacteria"], "Procryptobia", ["rnaspades", "blastn_megablast", "diamond_blastx", "bowtie2"])

## Cleaning step

makeblastdb -in both_Procyrptobia.fa -dbtype nucl -out both_Procyrptobia.DB

blastn -task megablast -query soft_filtered_transcripts_no_Chordata_no_Bacteria.fasta -db /Data/victoria/transcriptomes/Procryptobia/Procryptobia_both.fa.DB -outfmt 6 -num_threads 24 -evalue 1e-25 -max_target_seqs 1 -out c34_clean_vs_Procryptobia.blastnout

#get script from soyouz
cut -f 1 c34_clean_vs_Procryptobia.blastnout > Kinetoplastid_contigs.list
996 Kinetoplastid_contigs.list

perl /Data/victoria/scripts/lookup_reverse.pl soft_filtered_transcripts_no_Chordata_no_Bacteria.fasta Kinetoplastid_contigs.list

mv lookup_out.fasta Trinity_noBac_noPrey.fasta


* do this for psammosa pacifica

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

downloaded parasites:
amyloodinium_ocellatum

Corallicolids:
GZC_10102021
GZC_10142021
M6


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
file_rename_headers_phylo("soft_filtered_transcripts_no_Chordata_no_Bacteria.fasta.transdecoder.pep", "amyloodinium_ocellatum")
exit()
mv soft_filtered_transcripts_no_Chordata_no_Bacteria.fasta.transdecoder.pep.renamed amyloodinium_ocellatum.pep.fasta
cp amyloodinium_ocellatum.pep.fasta ../transcriptomes_to_add

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

```




### Combine datasets on Rosetta
# combine Ina's most current version of the 263 gene tree on google drive:
# 263genes_nov012021_IN
# 263_genes_24april2021_VJR

# check each folder has 263 genes using `ls | wc -l`

## script:

location on rosetta: `/Data/victoria/263_db/263_genes_comb_VJR_IN`


### Run combination script that adds both header and sequence **
    Example on how to use the new combine datasets julia script (rosetta)

      /Data/victoria/software/julia-1.5.3/bin/julia
      include("/Data/victoria/scripts/Transcriptome.jl"); using .Transcriptome;
      combine_datasets("/Data/victoria/263_db/263_genes_comb_VJR_IN/263_combined_out_false/", "/Data/victoria/263_db/263_genes_comb_VJR_IN/263_trees_28april2021_VJR", "/Data/victoria/263_db/263_genes_comb_VJR_IN/263genes_nov012021_IN"; check_id = true

### notes from previous run:
- This script will first use my dataset to make a list of all the headers and sequences.
      then it will add any header or sequence that does not match the one my
      dataset. Sometimes you will get copies of the same sequence, beacuse of slight changes
      in the header (font or format etc). My method to remove copies is
      - maintain order of your fasta addition (i.e. remove the newest duplicate)
      - you can see where your last added transcriptomes are and remove all
      the unwanted duplicates following yours. Do not change order of sequences.

***under construction**
## cleaned up version is on jezero at : 

I will collect all of Ina and My new sequences and add them to Liz's 263 curated database
- Liz has removed many old sequences, and fixed up paralogies :) 
- Use scafos to pull out my sequences and add to her dataset

### to do after this phylo update
Another thing I want to do: re-add transcriptomes using combined datasets (parasites from the same host infection)


