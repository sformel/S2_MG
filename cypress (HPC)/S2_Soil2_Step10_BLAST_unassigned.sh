#!/bin/bash

#SBATCH --job-name=S2_Step10_BLAST_unassigned
#SBATCH --output=S2_Step10_BLAST_unassigned.output
#SBATCH --error=S2_Step10_BLAST_unassigned.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#notes on script contained in "S2_Master_Scripts_Notes_11May17.txt"

module load ncbi-blast/2.5.0+

blastn -db /lustre/project/svanbael/steve/nt_db/nt -query /lustre/project/svanbael/steve/nt_db/rep_set.fna -outfmt '7 qseqid sseqid sskingdoms sscinames sblastnames length pident qstart qend sstart send evalue bitscore sgi sacc staxids sscinames' -num_threads 10 -out /lustre/project/svanbael/steve/unassigned_100_blast_results.txt
