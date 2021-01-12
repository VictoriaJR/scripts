#!/bin/bash

## SCRIPT: rnaspades.sh

## USAGE: For running the assembly program rnaspades for one or more samples on Cedar.
##        Will then move assembly results to a directory called "assemblies".

## INPUT: Forward and reverse trimmed raw reads in .fastq format, and a pre-made directory called "assemblies".

#SBATCH --time=01-20      		# Mandatory! DD-HH:MM:SS  or DD-HH:MM or HH:MM:SS or MM:SS or DD-HH or MM.
#SBATCH --account=def-keeling   # Mandatory!
#SBATCH --mem=180G           	# Memory requested. Can specify in G, M, K (binary prefixes). Default is 256M.
#SBATCH --cpus-per-task=40      # Number of cores requested. Default is probably 1?

# Output and Stderr
#SBATCH --job-name=rnaspades  		# Default is name of the batch script
#SBATCH --output=%x.out         # Default is "slurm-%j.out" where %j is jobid.
#SBATCH --error=%x.error        # Default is written to same file as standard out

# Mail Options
#SBATCH --mail-user=victoriak.jacko-reynolds@botany.ubc.ca      # Input your email
#SBATCH --mail-type=BEGIN,END,FAIL           # Receive an email when your script begins, ends, or fails

module load nixpkgs/16.09
module load gcc/7.3.0                  		#load your modules
module load spades/3.13.1

spades.py --rna --pe1-1 "$i"_R1_001.cutadapt.fastq --pe1-2 "$i"_R2_001.cutadapt.fastq --threads 40 --memory 180 -o "$i"_rnaspades_assembly ; done
