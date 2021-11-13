#!/bin/bash

#SBATCH --job-name=S2_fastq_to_fasta
#SBATCH --output=S2_fastq_to_fasta.output
#SBATCH --error=S2_fastq_to_fasta.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=2
#SBATCH --mem=64000 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#notes on script contained in "S2_Master_V4_Scripts_Notes_17July17.txt"

fastools -fastqtofasta -in /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed.assembled.fastq -out /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed_assembled.fasta