#!/bin/bash

#SBATCH --job-name=S2_fastqc
#SBATCH --output=S2_fastqc.output
#SBATCH --error=S2_fastqc.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#Running FastQC, a fastq quality control program written in java.  program also has a wrapper function to run it in the command line.  All the lines are commented out here so I can also run this as a job on cypress.

#18 May 2017
#by Steve Formel

#here I'm running it on cypress

#installation instructions for FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
#wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip --no-check-certificate

#unzip file
#unzip fastqc_v0.11.5.zip

#make sure to give the wrapper command, "fastqc" in the FastQC folder, execution permission
#chmod +x fastqc

#go into idev mode to run it, or run as a submitted job using SLURM script

#run batch analyses - if files are in .gz format, use "zcat" instead of "cat"
cat ./S2_fastq/*fastq | FastQC/fastqc stdin --outdir=/lustre/project/svanbael/steve/S2_fastq/FastQC_out



