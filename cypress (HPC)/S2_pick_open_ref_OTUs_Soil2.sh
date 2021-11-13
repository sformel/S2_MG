#!/bin/bash

#SBATCH --job-name=S2_S2_Soil2
#SBATCH --output=S2_S2_Soil2.output
#SBATCH --error=S2_S2_Soil2.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=10
#SBATCH --mem=64000 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#notes on script contained in "S2_Master_V4_Scripts_Notes_17July17.txt"

module load qiime/1.9.1

pick_open_reference_otus.py -i /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed_assembled.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2 -a -O 200 -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta --suppress_step4 --suppress_taxonomy_assignment --suppress_align_and_tree



