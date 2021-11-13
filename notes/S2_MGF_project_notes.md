# Notes describing the S2_MG Project.

### Last updated Nov 12, 2021
### by Stephen Formel
----------------------
## Project History

_This notebook contains references to 16S sequences, which ultimately weren't included in the paper: Spatial and temporal comparisons of salt marsh soil fungal communities following the deepwater horizon spill in Wetlands Ecology and Management.  DOI: 10.1007/s11273-021-09848-y_

_Thus the original name of the project was S2_MG, but was changed to S2_MGF when we decided to focus on Fungal analyses._

_Portions of the 16S sequences were used in another paper:_

_Rhizosphere microbial communities reflect genotypic and trait variation in a salt marsh ecosystem engineer
DOI: 10.1002/ajb2.1497_

 


This project was begun in 2010 or 2011 (I think) by Mike Blum and was funded by GOMRI (I believe in the first call for proposals).  I have not seen the original proposal, but to my understanding, it was somewhat different from what was ultimately achieved.  At some point grad student Elizabeth Jarrell began taking soil cores, but eventually the project was passed off to post-doctoral researcher, Dr. Demetra Kandalepas. Demetra finished taking soil cores and some samples of *Spartina alterniflora* as well.  She extracted DNA from a subset of the samples and sent them to a ACGT inc. for library prep and sequencing.  The work was done, but for various reasons was not analyzed until Kim Mighell and myself picked up the project in the spring of 2016.

Kim helped me get started in analyzing the sequences,  but both of us were pretty inexperienced so it took us some time to understand what we were doing.  Furthermore there were very few notes regarding the collection and processing of samples, but eventually I felt satisfied that I understood where and when every sample came from.

After several false starts, and a lot of dead ends, I finally reached a point where I was happy with the technical filtering and analysis of the sequences, around October 2017.  In February of 2018 I began writing the manuscript.

##Sequence processing

In the end I only worked with the soil samples.  From the ACGT report, it wasn't clear if files were already quality trimmed and filtered.  I wanted to make sure I was using raw files and not files that had been trimmed or filtered even though they were demultiplexed.  So I analyzed the sequences we received and compared them to the to the numbers in the report:

ATGC\_Table 8\_Read\_counts\_SF\_check.xls

#### EXAMINE READS FOR LENGTH AND QUALITY
note: Q scores are dependent on platform used.  Most modern sequencers use ASCII base 33 characters, but older illumina platforms used ASCII base 64.  I checked, these are ASCII base 33 because I saw characters like "." and "<" and "7" and ")" in my fastq files.

	https://en.wikipedia.org/wiki/FASTQ_format#Encoding
	http://drive5.com/usearch/manual/quality_score.html

#### CONCATENATE SEQUENCE FILES
	#unzip gz files: gunzip *.gz
	#merge seqs: merge fastq files from each sample into one file
		a. If all your fastq files are in individual folders, use this code.  The command "find" makes a recursive file list if your files are located in subfolders
			
			find |cat *.fastq > allseqs.allsamples.fastq
			
		b. if all your fastq files are in one folder use this code.  Here I've decided to make a merge of R1, R2, and all reads so I can examine the stats by forward and reverse read as well.
			
			cat *R1*.fastq > all_R1.allsamples.fastq  	#merge R1 reads
			cat *R2*.fastq > all_R2.allsamples.fastq	#merge R2 reads
			cat *.fastq > allseqs.allsamples.fastq	#merge all reads
			
#### ANALYZE FASTQ FILES

Analyze all R1 & R2 samples
each merge of all the fastq files (see above) is about 35 gigs

    FastQC/fastqc ./S2_fastq/all_R1.allsamples.fastq --outdir=/lustre/project/svanbael/steve/S2_fastq/FastQC_R1_out
    FastQC/fastqc ./S2_fastq/all_R2.allsamples.fastq --outdir=/lustre/project/svanbael/steve/S2_fastq/FastQC_R2_out
    
	#RUNTIME = 27 min



#### Count sequence length distribution using unix commands
oddly, only usearch seems to do this of all the fastq analysis programs, and my file is too big to work on the free 32-bit version of usearch.

But I found this nice awk command someone else wrote:

    # ~35 Gig file, took about 10 min to run on cypress login node
    awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' all_R1.allsamples.fastq >> all.R1.fastq_length_dist.txt 
    awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' all_R2.allsamples.fastq >> all.R2.fastq_length_dist.txt 
    
#### Graphs of seq distribution made in R with script "seq\_length\_distribution.R"

    #Script for making histogram of fastq sequence length distribution.
    #18 May 2017
    #by Steve Formel
    
    #Run in command line, not in R
      
      #Make distribution with the Unix command:
      #awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' all_R1.allsamples.fastq >> all.R1.fastq_length_dist.txt
    
    #Above took ~ 10 min to run on the cypress login node.
    
    setwd("C:/Users/Stephen and Moppa/Google Drive/Van Bael Google Drive/Boggs443 Data & User folders/Users/Grad Students/Steve/S2_11May2017/Results/FastQC_out/soil_only/")
    
    #read in data
    r1 <- read.table("all.R1.fastq_length_dist.txt", col.names = c("Seq_Length", "Count"), colClasses = c("numeric", "numeric"))
    r2 <- read.table("all.R2.fastq_length_dist.txt", col.names = c("Seq_Length", "Count"))
    
    #sort by seq_length
    r1.ordered <- r1[order(r1$Seq_Length),]
    r2.ordered <- r2[order(r2$Seq_Length),]
    
    #plot, not technically a histogram, but serves as a freq. dist all the same
    library(ggplot2)
    
    # Colorblind friendly palette
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    p <- ggplot(r1.ordered, aes(Seq_Length, Count, fill = Seq_Length))
    p + geom_bar(stat = "identity") + scale_y_log10() +ylab("log10(Count)") + labs(title = "S2 FastQ Read Length Distribution for all R1 reads")
    ggsave("r1_seq_length_dist.png", width = 10, height = 8, device = "png")
    
    p <- ggplot(r2.ordered, aes(Seq_Length, Count, fill = Seq_Length))
    p + geom_bar(stat = "identity") + scale_y_log10() +ylab("log10(Count)") + labs(title = "S2 FastQ Read Length Distribution for all R2 reads")
    ggsave("r2_seq_length_dist.png", width = 10, height = 8, device = "png")
    

##### Rename files
1. changed seq filenames using the unix command "rename" to read S2.1_R1.fastq instead of 1.R1.fastq.  Also changed the several sequences with longer names (named differently by the sequencing company for no apparent reason) to follow this convention.  Ultimately this step wasn't necessary, but I have left it in for documentation.

2. Using Candice's friend's program called FAST (https://github.com/ZeweiSong/FAST), I generated a mapping file for R1 and R2, ultimately this wasn't necessary either.	
	
	```fast.py -generate_mapping -i ./R1 -o read1_map.txt
	fast.py -generate_mapping -i ./R2 -o read2_map.txt```


3. add sample labels to all sequences

	```fast.py -add_labels -m read1_map.txt -i ./R1 -o read1_labeled -t 4
	fast.py -add_labels -m read2_map.txt -i ./R2 -o read2_labeled -t 4```

4. output file is named labeled\_S2.1\_R1.fastq or labeled\_S2.1\_R2.fastq

#### count number of seqs in labeled fastq files (all seqs)

		 awk '{s++}END{print s/4}' S2_R1_merged.fastq
			
			#69,348,965 sequences from 143 fastq files
 
		 awk '{s++}END{print s/4}' S2_R2_merged.fastq

			#69,348,965 sequences from 143 files

#### each individual fastq file

	for file in *.fastq
	do 
		awk '{s++}END{print s/4}' $file
	done

##### Results in file "ATGC\_Table8\_Read\_counts\_SF\_check.xls"

All files match table 8 read count for raw data except for:

58,76, and 81.  They have significantly fewer reads in the files we received than are listed on the table.  

I did not check plant sample read counts.  I checked the count number on the files given to us named "76.R1_val_1.fq" and found that this matched the trimmed reads in the Table 8 xls file delivered to us.  So this means for samples 58,76, and 81 our raw data does not match what the supposed output was.  But it is larger than the trimmed data, so I'm going to assume it is raw data for now.

#### Adapter and quality trimming.  

It was difficult for me to understand what I was trying to do.
Below are some thoughts and links I found useful:

1. Illumina 5' end sequencing cycle starts right from the start (5' end) of the actual sequence. So adapters sequence contamination in the final read would be only on 3' end.
2. Tutorial: http://www.ark-genomics.org/events-online-training-eu-training-course/adapter-and-quality-trimming-illumina-data
3. http://seqanswers.com/forums/showthread.php?t=17939
4. http://onetipperday.sterding.com/2012/08/three-ways-to-trim-adaptorprimer.html
5. explains how to trim "linked" primers with cutadapt:
	the second and fourth primers are the reverse complement of the primers in primers.csv  In theory the seq is only trimmed if it occurs with the first primer.  If left untrimmed, it is discarded. 

	https://github.com/marcelm/cutadapt/issues/237

###### To compare the overlap parameter in cutadapt, I tested 3 and 5 for the sample #100 reads:

	Overlap 5:
		
	=== Summary ===

	Total read pairs processed:            348,289
  	Read 1 with adapter:                 106,831 (30.7%)
  	Read 2 with adapter:                  96,006 (27.6%)
	Pairs that were too short:              10,871 (3.1%)
	Pairs written (passing filters):        37,511 (10.8%)

	Overlap 3:


	=== Summary ===

	Total read pairs processed:            348,289
  	Read 1 with adapter:                 106,831 (30.7%)
  	Read 2 with adapter:                  96,006 (27.6%)
	Pairs that were too short:              10,871 (3.1%)
	Pairs written (passing filters):        37,511 (10.8%)

Didn't seem to make a difference and 5 is considered more conservative.  I didn't have enough space on the cluster to run both for the entire sample set.

##### Install cutadapt
	CUTADAPT - http://cutadapt.readthedocs.io/en/stable/guide.html
					
	To install cutadapt - from Candice Lumibao's lesson in Lorena Torres' genomics class

			conda create -n software --clone=/share/apps/anaconda/2/2.5.0
			source activate software
			conda install -c bioconda cutadapt
			#https://wiki.hpc.tulane.edu/trac/wiki/cypress/AnacondaInstallPackage
			conda remove conda-build
			conda remove conda-env
				
			module load cutadapt

##### Load cutadapt
You may need to add cutadapt variables to your path

	module load anaconda
	export CONDA_ENVS_PATH=/lustre/project/svanbael/steve/software/conda-env
	source activate software

##### Run cutadapt

	This is cutadapt 1.14 with Python 2.7.13
	Command line parameters: -a AYTGGGYDTAAAGNG...GGATTAGATACCCBNGTA -a ACCTGCGGARGGATCA...AACTTTYARCAAYGGATCTC -A TACNVGGGTATCTAATCC...CNCTTTAHRCCCART -A GAGATCCRTTGYTRAAAGTT...TGATCCYTCCGCAGGT --discard-untrimmed --match-read-wildcards -e 0.1 -O 5 -m 50 -o S2_R1_merged.cut.fastq -p S2_R2_merged.cut.fastq S2_R1_merged.fastq S2_R2_merged.fastq
	Trimming 4 adapters with at most 10.0% errors in paired-end mode ...
	Finished in 3700.14 s (53 us/read; 1.12 M reads/minute).

	=== Summary ===

	Total read pairs processed:         69,348,965
  	Read 1 with adapter:              56,060,960 (80.8%)
  	Read 2 with adapter:              51,928,162 (74.9%)
	Pairs that were too short:           1,460,564 (2.1%)
	Pairs written (passing filters):    46,955,072 (67.7%)

	Total basepairs processed: 20,432,280,128 bp
  	Read 1: 10,232,982,504 bp
  	Read 2: 10,199,297,624 bp
	Total written (filtered):  12,433,567,385 bp (60.9%)
  	Read 1: 6,308,312,147 bp
  	Read 2: 6,125,255,238 bp

##### Used FastQC to examine quality profile of trimmed seqs

	./software/FastQC/fastqc ./seqfiles/soil_only/cleaned_seqs/S2_R1_merged.cut.fastq --outdir=/lustre/project/svanbael/steve/seqfiles/soil_only/cleaned_seqs/FastQC_R1_trimmed_out

	./software/FastQC/fastqc ./seqfiles/soil_only/cleaned_seqs/S2_R2_merged.cut.fastq --outdir=/lustre/project/svanbael/steve/seqfiles/soil_only/cleaned_seqs/FastQC_R2_trimmed_out


#### Used PEAR to align seqs

Note: this is so fast that I just ran it in idev, not a script.

	pear-0.9.10-bin-64 -f ./S2_R1_merged.cut.fastq -r ./S2_R2_merged.cut.fastq -o ./S2_trimmed -n 100 -j 20

	#Started at 3:08pm ended at 3:18pm

	PEAR v0.9.10 [May 30, 2016]

	Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: ./S2_R1_merged.cut.fastq
    Reverse reads file.................: ./S2_R2_merged.cut.fastq
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 100
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 20
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.230545
      C: 0.263754
      G: 0.252193
      T: 0.253507
      464961 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 45,708,023 / 46,955,072 (97.344%)
    Discarded reads ...................: 0 / 46,955,072 (0.000%)
    Not assembled reads ...............: 1,247,049 / 46,955,072 (2.656%)
    	Assembled reads file...............: ./S2_trimmed.assembled.fastq
    	Discarded reads file...............: ./S2_trimmed.discarded.fastq
    	Unassembled forward reads file.....: ./S2_trimmed.unassembled.forward.fastq
    	Unassembled reverse reads file.....: ./S2_trimmed.unassembled.reverse.fastq

#### Convert fastq to fasta file for use in QIIME

	Convert to fasta with fastools (add fastools to your path after installing)

		fastools -fastqtofasta -in /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed.assembled.fastq -out /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed_assembled.fasta

#### Split results by Sample for Statistics in QIIME

	module load qiim/1.9.1
	split_sequence_file_on_sample_ids.py -i ./S2_trimmed_assembled.fasta -o ./	--file_type fasta

##### count reads per sample

	for fasta files: grep -c "^>" file.fasta
	for fastq files: awk '{s++}END{print s/4}' file.fastq

#### Compare My results with ACGT

File: "S2\_seq\_analysis_1Aug17.R"

![](https://i.imgur.com/2J78C4V.png)
![](https://i.imgur.com/Ymn7SMq.png)


#### Seq Stats to report in paper

143 Soil Samples


Raw reads pairs

	Total: 69,348,965
	Avg quality per read: Q35
	Avg read length: 150
	Avg reads per sample: 484,947.79 
	Median reads per sample: 400,643
	Range: 41,262 - 2,821,283

After trimming and filtering and merging pairs

	45,708,023
	Avg quality per read: Q37
	Avg read length: 207
	Avg reads per sample: 319,636.5
	Median reads per sample: 222,187
	Range: 19494 - 2,641,849
	
	

## Picking OTUs in QIIME

#### CHECK QIIME MAP FILE

	validate_mapping_file.py -m QIIME_MAP_S2_metagenomics_sample_map_11May17.txt -o S2_map_check

	#only warnings are for the fake barcodes


-------------------------------------------------------------------------------
#### PIPELINE DESCRIPTION

Note that the combined databases is because ITS and 16S were sequenced together and I was unable to satisfactorily separate them with bioinformatics.

1. 	Pick OTUs - using QIIME default pick_open_reference_otus.py (uclust) against a combined GreenGenes and UNITE database
2. 	Assign Taxonomy - using QIIME default assign_taxonomy.py (uclust) against a combined GreenGenes and UNITE database
3. 	Chimera Filtering - using VSEARCH
4. 	Make OTU table - using QIIME make_otu_tables.py
5. 	Split OTU tables into Fungi, Bacteria, Archaea, and Unassigned - using QIIME filter_taxa_from_otu_table.py
6. 	Summarize OTU tables - using biom summarize-table
7. 	Filter out sequences with fewer than 100 Seqs
8. BLAST "Unassigned" taxa against NCBI's nt database.

12. Ecological analysis in QIIME, R(vegan), R(phyloseq)

Note that the QIIME pipeline typically builds trees for phylogenetic analysis, but I couldn't do this because I'm working with Fungi also.

-----------------------------------------------------------------------------
Note: Below are the notes and code for what I used in my final analysis.  I also did a comparison to understand how different decisions affected the outcome of picking OTUs.  The results are in the file:
	
	S2_comparison_of_bioinformatics.xlsx

Ultimately I went with the result named "Soil2".

##### 1.	PICK OTUs

######Make combined reference database

	#combine GreenGenes and UNITE databases

	cat ./GG_db/gg_13_8_otus/rep_set/97_otus.fasta ./UNITE_db/sh_refs_qiime_ver7_97_20.11.2016.fasta >> ./combined_db/GG_UNITE_db_20Nov2016.fasta
			
This next step is run as a SLURM script ("S2_pick_open_ref_OTUs_Soil2.sh") in QIIME on Tulane's HPC "Cypress" The version of QIIME they have loaded is 1.9.1  Here is notes about running QIIME jobs in parallel: [http://qiime.org/tutorials/parallel_qiime.html](http://qiime.org/tutorials/parallel_qiime.html)

Important:  I had to modify part of the qiime_config file

	Use an editor like nano to change the "cluster_jobs_fp" parameter from "start_parallel_jobs.py" to "start_parallel_jobs_slurm.py" in order for parallel functions to actually run in parallel.

Notes about certain parameters in the OTU picking function (AKA flags)

1. --suppress\_step4 is an option to remove a denovo picking step of the algorithm.  Kim and others have found that we can't get enough time on Cypress for it to complete this step.
2. --suppress\_taxonomy\_assignment was necessary because I didn't understand how to pass in the combined taxonomy database at the time. 
3. --suppress\_align\_and\_tree, same as above.

command from script:

	`pick_open_reference_otus.py -i /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed_assembled.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2 -a -O 200 -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta --suppress_step4 --suppress_taxonomy_assignment --suppress_align_and_tree`

	#Runtime: (hours:minutes:seconds) 07:40:37

Note:  The default algorithm requires each OTU to contain at least 2 sequences, any OTUs which failed to meet this criteria would be removed from the final_otu_map.txt to produce: final_otu_map_mc2.txt

The final_otu_map_mc2.txt is used to build the final representative set: rep_set.fna
	
---------------------------------------------------------------------------
##### 2.	ASSIGN TAXONOMY

###### make GG and UNITE taxonomy database 

1. first I added a space (by find and replace in a text editor) in UNITE taxa labels so that it would match Greengenes.  Not sure if this was necessary, but I didn't want it to cause problems.
	
		cat 97_otu_taxonomy.txt sh_taxonomy_qiime_ver7_97_20.11.2016_with_space_after_semicolon.txt >>GG_UNITE_tax_db.txt

2. Originally I tried to run this in parallel but the parallel version of this command kept hanging up, but this only took ~ 3 min on the cypress login node  (probably should have been run in idev)
	
    `assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2/Step2`

------------------------------------------------------------------------------------
##### 3.	CHIMERA FILTERING

After doing some searching, it seems that no one can decide on a good way to check for chimeras.

	#https://github.com/torognes/swarm/issues/43
	#http://drive5.com/usearch/manual/uchime2_algo.html
	#https://groups.google.com/forum/#!topic/vsearch-forum/fTjFTXJq3BU
	#https://github.com/torognes/swarm/wiki/Frequently-Asked-Questions#when-clustering-with-swarm-when-is-an-appropriate-time-to-check-for-chimeras
	#https://groups.google.com/forum/#!searchin/vsearch-forum/chimera%7Csort:relevance/vsearch-forum/Hh-AaYVTQpg/C7HqESYrDQAJ
	#Good ideas: https://groups.google.com/forum/#!topic/qiime-forum/zRiEF3wmtQo
	#http://geoffreyzahn.com/getting-started-with-qiime-for-fungal-its-cleaning-its-reads-and-picking-otus/

But I like the idea of using vsearch, an opensource version of usearch.

	vsearch v2.4.3_linux_x86_64, 63.0GB RAM, 40 cores
	https://github.com/torognes/vsearch  

First had to combine the representative sets for GG and UNITE.  UNITE has a UCHIME chimera database here: [https://unite.ut.ee/repository.php#uchime](https://unite.ut.ee/repository.php#uchime)  

Download it and combine it with the unaligned GG 97_otus.fasta

	cat uchime_reference_dataset_01.12.2016.ITS1.fasta 97_otus.fasta >> uchime_GG_chimera_db.fasta

Don't forget to add vsearch bin to your PATH!

	'mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/Step3

	vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil2/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil2/Step3/Soil2_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil2/Step3/Soil2_non_chimeras.fna --strand plus

Yielded  2284 chimeras

Note about VSEARCH in parallel: 

Many of the VSEARCH commands support multiple threads. It is therefore generally a good idea to run it on a computer with several cores. **It does not support communication between nodes (e.g. MPI) in a cluster though.**  By default VSEARCH will start a number of worker threads equal to the number of cores detected in the computer. To override this, use the "--threads" option with the number of threads as the argument.

RUNTIME: 30 secs

---------------------------------------------------------------------------
###### 4.Make OTU Tables

Run in cypress login node
	
There was a glitch:

An OTU that was assigned from UNITE has an e with a dierisis (umlaut), and QIIME doesn't like it.  So to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/ASCII_find_Soil2.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments_clean.txt
	
After that was fixed:


	mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables

	make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil2/Step3/Soil2_chimeras.fna

---------------------------------------------------------------------------	
##### 5. Break OTU Table into Bacteria, Fungi, Archaea, and "Unassigned"
		
remove mitochondria & chloroplast
	
	filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

break out tables by taxa
	
Unassigned OTUs
	
	filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass.biom -p Unassigned
	
16S (bacteria and archaea)	

	filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S.biom -p k__bacteria,k__archaea  
	
fungi

	filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS.biom -p k__fungi
	
---------------------------------------------------------------------------
##### 6. Make summaries of OTU Tables

Note:  	Num observations = total number of OTUs
 		Total count = total number of seqs

	mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries
	
Unassigned OTUs

	biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_unass.txt

	Results:

	Num samples: 143
	Num observations: 2115
	Total count: 12549002
	Table density (fraction of non-zero values): 0.220
	
16S OTUs

	biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_16S.txt

	Results:

	Num samples: 143
	Num observations: 18992
	Total count: 18372382
	Table density (fraction of non-zero values): 0.237	

Fungal OTUs

	biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_ITS.txt

	Results: 

    Num samples: 143
    Num observations: 535
    Total count: 1867258
    Table density (fraction of non-zero values): 0.099
-------------------------------------------------------------------------
##### 7. FILTER OUT SEQUENCES WITH FEWER THAN 100 SEQS

Unassigned

	filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass_100.biom

	Results:

	Num samples: 143
	Num observations: 1323
	Total count: 12525366
	Table density (fraction of non-zero values): 0.316

16S

	filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S_100.biom

	Results:

	Num samples: 143
	Num observations: 8117
	Total count: 18134630
	Table density (fraction of non-zero values): 0.472

Fungi

	filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS_100.biom
		
	Results:
	
	Num samples: 143
	Num observations: 194
	Total count: 1860060
	Table density (fraction of non-zero values): 0.233

-------------------------------------------------------------------------
##### 8. BLAST Unassigned OTUs against NCBI nt database

######Make local ncbi blast database to explore unassigned otus

download files.  Need 34G of space

		wget ftp://ftp.ncbi.nih.gov/blast/db/nt.*tar.gz
		wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz #This is the taxid for all entries		
extract files, tar can't batch extract, you need to do it in a loop.  

		for i in *.tar.gz ; do tar -zxvf $i ; done  

As long as the files are all in the same directory, the database should already be "made", you don't need to run the command "makeblastdb"

		#Also add the blast db to your path
		nano ~/.bash_profile

		export BLASTDB=$BLASTDB:/media/Data/path/to/your/database/files/

filter representative set fasta using biom tables to only include representative seqs from Unassigned OTUs - run in login node
	
		filter_fasta.py -f ./Soil2/rep_set.fna -b ./Soil2/biom_tables/Soil2_unass_100.biom -o ./unassigned_rep_set.fna

BLAST - run as script "S2\_Step10\_BLAST\_unassigned.sh"

*Note:  At the time, this was step 10 in the process, not step 8.* 

For description of output format options: [https://www.ncbi.nlm.nih.gov/books/NBK279675/](https://www.ncbi.nlm.nih.gov/books/NBK279675/)

	module load ncbi-blast/2.5.0+

	#it looks like you have to be explicit in your filepaths

	blastn -db /lustre/project/svanbael/steve/nt_db/nt -query /lustre/project/svanbael/steve/nt_db/rep_set.fna -outfmt '7 qseqid sseqid sskingdoms sscinames sblastnames length pident qstart qend sstart send evalue bitscore sgi sacc staxids sscinames' -num_threads 10 -out /lustre/project/svanbael/steve/unassigned_100_blast_results.txt

###### filter blast output by column

filter output  $X is column number X

	awk -F'\t' '{ if ( $1 ~ /^#/ || ($6>=150 && $7 >=90)) {print}}' ./unassigned_100_blast_results.txt > ./filtered_blast_align_100_ident_90.txt

filter for top hits only (much easier to look at)
	awk '!x[$1]++' ./filtered_blast_align_150_ident_90.txt > top_blast_hit_only.txt
			 
#RERUN THE ABOVE WITH NEW BLAST SCRIPT
-------------------------------------------------------------------------------------------			 
#Steve pick up here
####Generate Rarefaction Curves

Analysis in R using phyloseq and custom functions borrowed from other people.  This is a minimal version of the code, the file below has more comments and specific filepaths.

Filename: "make_rarefaction curve.R"

    #set working directory & load libraries------ 
    
    	library(phyloseq)
    	library(ggplot2)
		library(tidyr)
	    library(Rmisc)
	    library(cowplot)
    
    #import and clean data-------------
    
    	source("S2_import_clean_data_used_in_final_analysis.R")
    
    #Make curves-----
    
    	#taken from: https://github.com/joey711/phyloseq/issues/143
    
    #rename phyloseq object for rarefaction function
	    psdata <- fung.2season
	    
	    set.seed(42)
	    
	    calculate_rarefaction_curves <- function(psdata, measures, depths) {
	      require('plyr') # ldply
	      require('reshape2') # melt
	      
	      estimate_rarified_richness <- function(psdata, measures, depth) {
	    if(max(sample_sums(psdata)) < depth) return()
	    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
	    
	    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
	    
	    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
	    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
	    
	    molten_alpha_diversity
	      }
	      
	      names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
	      rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
      
      # convert Depth from factor to numeric
	      rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
	      
	      rarefaction_curve_data
	    }
    
    #make list of reads per sample up to 100,000 reads - This was added by Steve, change the number in the line generating "sample_list_100000" if you want a lower or higher number.
	    sample_list <- as.list(sort(sample_sums(psdata)))
	    sample_list_100000 <- unlist(sample_list[sort(sample_sums(psdata))<=100000])
    
    #make curve using actual read numbers of samples 
	    rarefaction_curve_data <- calculate_rarefaction_curves(psdata, c('Observed'), rep(c(1, 10, 100, min(sort(sample_sums(psdata))),sample_list_100000), each = 10))
	    
    #make summary
	    rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
	    
    #add sample data
	    rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')
	    
    #plot results - messy version meant for exploration
    
	    ggplot(
	      data = rarefaction_curve_data_summary_verbose,
	      mapping = aes(
	    x = Depth,
	    y = Alpha_diversity_mean,
	    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
	    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
	    colour = season,
	    group = X.SampleID
	      )
	    ) + geom_line(
	    ) + facet_wrap(~site)
	    
      
	#Change site and season to whatever variables you're interested in exploring
	    df2 <- summarySE(data = rarefaction_curve_data_summary_verbose, measurevar = "Alpha_diversity_mean", groupvars = c('site', 'season', 'Depth'))
	 
	#plot rarefaction curves   
	    p1 <- ggplot(data = df2, aes(x = Depth, y = Alpha_diversity_mean, color = season)) +
	      geom_point() +
	      geom_line() +
	      theme_bw() +
	      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	      facet_wrap(~ site)
	 #plot number of samples included as a function of rarefaction depth
   
	    p2 <- ggplot(data = df2, aes(x = Depth, y = N, color = season)) +
	      geom_point() +
	      geom_line() +
	      theme_bw() +
	      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	      scale_y_continuous(breaks = c(10:25)) +
	      facet_wrap(~ site)
	
	#being the plotting voodoo with the package cowplot	    
	    p1.2 <- p1 + theme(legend.position = "none")
	    p2.2 <- p2 + theme(legend.position = "none")
	    legend <- get_legend(p1)
	    
	    plots1 <- align_plots(p1.2,p2.2, align = 'v', axis = '1')
	    
	    top_row <- plot_grid(plots1[[1]],plots1[[2]], labels = c("A","B"))
	    
	    p <- plot_grid(top_row, legend, ncol = 2, rel_widths = c(4,0.5))
	    
	# now add the title
	    title <- ggdraw() + draw_label("Rarefaction Curves (A) and Number of Samples included (B) for Fungi ", fontface='bold')
	    plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
	    
	#don't forget to change file name according to input phyloseq object!
	    
	    ggsave("./images/manuscript/rare_curves_fung_2season.png", dpi = 300, width = 10, height = 6, units = "in", device = "png")

#### Results of Rarefaction Curves

####### Note the below is only for the winter and summer 2013 dates that were used in the paper. 


![](https://i.imgur.com/H8Tlpwn.png)
![](https://i.imgur.com/oIZnHu7.png)
------------------------------------------------------------------------------------------
S13. Make alpha rarefaction curves in QIIME


		#After looking at OTU table summary, 20,000 sequences would be a good cutoff for rarefaction.  I would throw out 7 Fourchon and 2 Bay Jimmy samples.
			
			filter_samples_from_otu_table.py -i SOIL_ONLY_open_ref_no4_no_unassigned_nochimera_otu_table.biom -o SOIL_ONLY_open_ref_no4_no_unassigned_nochimera_otu_table_20000.biom -n 20000

		#run as script in cypress - "S2_Step13_alpha)rare.sh" - for this to work you must add a config file because cypress has no graphics generator
		
					#add config file: https://groups.google.com/forum/#!searchin/qiime-forum/no$20$24DISPLAY$20/qiime-forum/6srHDScyh_0/yVJQfjo1V-EJ
					
						nano ~/.config/matplotlib/matplotlibrc
						
							1. type (without quotes) "backend : agg"
							2. save and exit

			alpha_rarefaction.py -i ../biom_tables/soil_only/SOIL_ONLY_open_ref_no4_no_unassigned_nochimera_otu_table_20000.biom -m QIIME_MAP_S2_metagenomics_sample_map_11May17.txt -o ../alpha_rare_20000 -a -O 20 -n 20 -e 20000
			
				#params (saved as separate file "alpha_params.txt"
					
					make_rarefaction_plots:d 300 #output hi-res plot
					alpha_diversity:m chao1,observed_otus,simpson,shannon,simpson_e,simpson_reciprocal #diversity metrics to calculate, important not to have spaces after commas
			
			RUNTIME = 14 min
			
		#Do the same for the bacteria/archaea biom table
			
			filter_samples_from_otu_table.py -i SOIL_ONLY_open_ref_no4_bact_archaea_nochimera_otu_table.biom -o SOIL_ONLY_open_ref_no4_bact_archaea_nochimera_otu_table_22000.biom -n 22000
			
			alpha_rarefaction.py -i ../biom_tables/soil_only/SOIL_ONLY_open_ref_no4_bact_archaea_nochimera_otu_table_22000.biom -m ../QIIME_MAP_S2_metagenomics_sample_map_11May17.txt -o ../alpha_rare_bac_22000 -a -O 20 -e 22000 -p alpha_params.txt
			
			RUNTIME = 12 min
			
		#do the same for the fungal biom table
 			
			filter_samples_from_otu_table.py -i SOIL_ONLY_open_ref_no4_fungi_nochimera_otu_table.biom -o SOIL_ONLY_open_ref_no4_fungi_nochimera_otu_table_10000.biom -n 10000

			alpha_rarefaction.py -i ../biom_tables/soil_only/SOIL_ONLY_open_ref_no4_fungi_nochimera_otu_table_10000.biom -m ../QIIME_MAP_S2_metagenomics_sample_map_11May17.txt -o ../alpha_rare_fun_10000 -a -O 20 -e 10000 -p alpha_params.txt

			RUNTIME = 10 min
			
#### R analysis

The remaining analysis was conducted in R and songbird.  See the R folder of scripts for additional information.  Songbird notes are below.


### Songbird - Multinomial Regression

See: 

Morton, J.T., Marotz, C., Washburne, A., Silverman, J., Zaramela, L.S., Edlund, A., Zengler, K., and Knight, R. (2019). Establishing microbial composition measurement standards with reference frames. Nature Communications 10, 2719.
 
##### Install on Cypress

On Cypress had to make the condarc file point to a new folder with more space for packages on lustre
https://docs.anaconda.com/anaconda/user-guide/tasks/shared-pkg-cache/

	conda create -n songbird_env python=3.7.3
	unset PYTHONPATH

	source activate songbird_env

	conda install songbird -c conda-forge

#had to fix this bug

https://github.com/matplotlib/matplotlib/issues/12439/


Run model in folder called songbird in the S2 project 

	songbird multinomial \
	--input-biom S2_MG/songbird/Soil2_16S.biom \
	--metadata-file S2_MG/songbird/S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
	--formula "site+season+site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--training-column Testing \
	--summary-interval 1 \
	--summary-dir results

Not working test with package data downloaded from github

	songbird multinomial \
	--input-biom songbird-master/data/redsea/redsea.biom \
	--metadata-file songbird-master/data/redsea/redsea_metadata.txt \
	--formula "Depth+Temperature+Salinity+Oxygen+Fluorescence+Nitrate" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir results


####See this

https://github.com/biocore/songbird/issues/128

Downgrade pandas

	conda install pandas==0.25.3


# New try On Mac laptop


	conda create -n songbird_env python=3.7.3
	unset PYTHONPATH

	source activate songbird_env

	conda install songbird -c conda-forge

	conda install pandas==0.25.3

#load env
	
	source activate songbird_env
	
	cd /Users/stephenformel/Google Drive/VB_lab/VBL_data/Collections, Projects/Collections/S2_MG collection/songbird


Test worked well

	songbird multinomial \
		--input-biom songbird-master/data/redsea/redsea.biom \
		--metadata-file songbird-master/data/redsea/redsea_metadata.txt \
		--formula "Depth+Temperature+Salinity+Oxygen+Fluorescence+Nitrate" \
		--epochs 10000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--summary-dir result

Run S2 simple model

	songbird multinomial \
		--input-biom Soil2_16S_100.biom \
		--metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
		--formula "site+season" \
		--epochs 10000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--summary-dir S2_results
		
In S2_16S_result folder:

	tensorboard --logdir .
	
Results look ok, but loss (graph 2) is very chunky, unstabler and not close to 0.  So I'm going to increase epochs by an order of magnitude and include an interaction term.  


	songbird multinomial \
		--input-biom Soil2_16S_100.biom \
		--metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
		--formula "site+season + site*season" \
		--epochs 10000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--summary-dir S2_16s_results
		
Version with increase in magnitude

	songbird multinomial \
		--input-biom Soil2_16S_100.biom \
		--metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
		--formula "site+season + site*season" \
		--epochs 100000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--summary-dir S2_16s_results
		
It turns out that Tensor can compare different runs as long as they are all in the same root folder.  So I'm going to start with my null model.  Anyway, doing 1 million iterations didn't seem to solve it. Note that it's important to make a subfolder into which your results run.  Aparently tensoflow:

"This is a known issue, TensorBoard doesn't like it when you write multiple event files from separate runs in the same directory. It will be fixed if you use a new subdirectory for every run (new hyperparameters = new subdirectory)."



	songbird multinomial \
		--input-biom Soil2_16S_100.biom \
		--metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
		--formula "1" \
		--epochs 10000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--summary-dir logdir/S2_NULL

site + season

	songbird multinomial \
		--input-biom Soil2_16S_100.biom \
		--metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
		--formula "site+season" \
		--epochs 10000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--summary-dir logdir/S2_site_season
		
site*season

	songbird multinomial \
    --input-biom Soil2_16S_100.biom \
    --metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \
    --formula "site+season + site*season" \
    --epochs 10000 \
    --differential-prior 0.5 \
    --summary-interval 1 \
    --summary-dir logdir/S2_sitexseason

Enter log directory
    
	cd logdir
	tensorboard --logdir .
	
Produced really noisy loss graphs.  Based off this conversation, I should reduce the learning rate:

https://forum.qiime2.org/t/songbird-loss-function-plot-looks-like-a-bad-ekg/14329/11

	songbird multinomial \
    --input-biom Soil2_16S_100.biom \
    --metadata-file S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt \--formula "site+season" \
    --epochs 10000 \
    --differential-prior 0.5 \
    --summary-interval 1 \
    --summary-dir logdir/S2_site_season_lowlearn \
    --learning-rate 0.00001
	
Didn't imrpove the loss beyond the Null, although it was a better model than the null.  But that took a long time to reach so I'm going to run it again with more epochs.  Then I realized my data included many timepoints that I was technically throwing in different seasons (so convoluting the data). So I subset the data to the 4 time points I was interested in and added a Test column.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_NULL
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_site_season
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_site_x_season
	
Definitely improved things.  Not sure if the interaction makes the model better, but it would be good info to know.  Try increasing the epochs by 5 times.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 150000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_site_x_season_epochx5

The epochs didn't really improve the model.  The loss is still a little high, and then I realized I had probably not specified enough samlpes for testing.  I had sepcified 8, and the rule of thumb is 10-20%.  With 8 I was a little shy of 10%.  So I specified 4 more for a total of 12 and then Reran the above models.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_NULL_batch12

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_site_season_batch12

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_site_x_season_batch12
	
It doesn't look like the interaction is improving the model much, but it's an interesting part of the hypothesis.  Since it's not making it worse, I'm going to include it.  I'm also going to try more epochs and a lower learning rate to see if that improves results at all.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 50000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir logdir/S2_site_x_season_batch12_epochx5

I originally ran the low learning model with 10000 epochs and it wasn't enough to approach a plateau, so I ran it again with 50,000 epochs
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 50000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--learning-rate 0.00001 \
	--summary-dir logdir/S2_site_x_season_batch12_lowlearn
	
So this didn't improve things and then I noticed that the "random test samples" category (in the HPARAMS part of Tensorboard) was only 5 for all models.  I thought I had specified it in the testing column!  Turns out that I'm a dummy.  The testing column is if you want to restrict which samples are used for training and which for testing.  I'm not sure if I do because it might be best if it was random.  But I've also set aside a nice balanced set for testing.  Anyway, if I want to change the test batch size I do that with:

	--num-random-test-examples
	
And then there is

	--batch-size
	
which is the number of samples to be evaluated per training iteration.  Both default to 5.  Which means all my "batch 12" stuff is meaningless.  The reason it didn't improve the model much was because all I was doing was specifying a balanced training set. 

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--summary-dir logdir/S2_site_x_season_batch12_test12

This one below seemed to run much faster than the others
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_batch12_realbatch12
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_batch12_test_and_batch12
	
Although the above model is sufficiently good, I want to see what happens if I adjust the prior.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season" \
	--epochs 10000 \
	--differential-prior 0.9 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_batch12_test_and_batch12_highprior
	
The prior actually didn't change the model much, although the loss was slightly worse.  Enough that I'm not going to use it, but good to know. Now I want to see what happens if I add in oil.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site + season + site*season + Total_relevant_PAHs" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_batch12_test_and_batch12_PAH
	
Wow, loss got down to almost zero.  But the curve was rather jagged.  

I did some reading and convinced myself that scaling is probably a good idea for continuous variables also.  

https://github.com/biocore/songbird/issues/135

It also looks like I could have used site\*season, just like in R (i.e. it expands to site + season + site*season).  Which might help, because my differential files haven't been including three levels for site x season.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site*season + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_batch12_test_and_batch12_PAH_newform
	
Hmm, this ended up making the model crummy again (probably scaling the oil, and reducing it's extreme values as well as introducing some NA values.)  But it also didn't fix the number of columns the model is returning from the interaction term.  So I'm going to try a test with the readsea data

	songbird multinomial \
	--input-biom songbird-master/data/redsea/redsea.biom \
	--metadata-file songbird-master/data/redsea/redsea_metadata_SFmod.txt \
	--formula "Dummy1*Dummy2" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir int_test_results
	
Ok, that test showed that it's something wrong with my metadata that is keeping it from being read as nominal.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "site*season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_batch12_test_and_batch12_newform


Figured it out.  site\*season doesn't return all combinations of variables (k-1 levels).  Turns out Patsy (the backend for the formulas is much more complicated than that (but for good and interesting reasons!)  This is the way I came up with the specifically get the combinations I need.

	songbird multinomial \
	--input-biom songbird-master/data/redsea/redsea.biom \
	--metadata-file songbird-master/data/redsea/redsea_metadata_SFmod.txt \
	--formula "Dummy1 + \
	Dummy2 + \
	C(Dummy1, Treatment('Top')):Dummy2 + \
	C(Dummy1, Treatment('Bottom')):Dummy2" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--summary-dir int_test_results
	
So this should work, but I'm going to verify its legitimacy on the songbird github

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + \
	C(season, Treatment('WINTER')) +\
	C(site, Treatment('BJ')):season + \
	C(site, Treatment('F')):season" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_newform

With Oil

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) +\
	C(site, Treatment('BJ')):season + \
	C(site, Treatment('F')):season + \
	Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_newform_oil
	
NULL

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/S2_site_x_season_newform_NULL
	
	
So I asked the author of songbird about it and he gave me a nice explanation:  https://github.com/biocore/songbird/issues/136

I also realized I need to filter out the two outliers (S2.36, S2.23) I didn't include in the R analysis, so I did that in QIIME on Cypress.

	filter_samples_from_otu_table.py -i Soil2_16S_100_songbird.biom -o Soil2_16S_100_no_outliers_songbird.biom -m S2_metagenomics_sample_map_songbird.txt --output_mapping_fp S2_metagenomics_sample_map__no_outliers_songbird.txt -s 'sampleID:*,!S2.36, !S2.23'
	
	filter_samples_from_otu_table.py -i Soil2_ITS_100_songbird.biom -o Soil2_ITS_100_no_outliers_songbird.biom -m S2_metagenomics_sample_map_songbird.txt --output_mapping_fp S2_metagenomics_sample_map__no_outliers_songbird.txt -s 'sampleID:*,!S2.36, !S2.23'
	
Here is the final model selection (null, with and without oil) using BJ and Winter as the reference treatments.

	songbird multinomial \
	--input-biom Soil2_16S_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/final_form/S2_MG_bac_NULL
	
	songbird multinomial \
	--input-biom Soil2_16S_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ'))*C(season, Treatment('WINTER'))" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/final_form/S2_MG_bac_site_x_season
	
	songbird multinomial \
	--input-biom Soil2_16S_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ'))*C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/final_form/S2_MG_bac_site_x_season_with_oil

	songbird multinomial \
	--input-biom Soil2_16S_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ'))*C(season, Treatment('WINTER'))" \
	--epochs 50000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--learning-rate 0.00001 \
	--summary-dir logdir/final_form/S2_MG_bac_site_x_season_lowlearn

Low learning rate didn't improve the model.  I'm going to go with site x season, as oil had slightly better loss values, but slightly worse error values, and so it doesn't seem like adding a term to the model improves the model overall.

Same site*season model for fungi

	songbird multinomial \
	--input-biom Soil2_ITS_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ'))*C(season, Treatment('WINTER'))" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/final_form/S2_MG_fung_site_x_season
	
	songbird multinomial \
	--input-biom Soil2_ITS_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/final_form/S2_MG_fung_NULL
	
	songbird multinomial \
	--input-biom Soil2_ITS_100_no_outliers_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_no_outliers_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ'))*C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/final_form/S2_MG_fung_site_x_season_oil
	
Like bacteria, oil did not improve the model

##### Sep 22, 2020

It looks like I lost my notes on the final models when the google drive re-synced.  Sigh.  I'll recreate them and run them here and check the results.

Also, I realized that log transforming the PAHs might create be a more linear relationship.
	
##### Bacteria Standalone
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_2/S2_MG_bac_null
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER'))" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_2/bac/S2_MG_bac_site_season
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_2/S2_MG_bac_site_season_oil

This last one took a really long time to fit.  It may have been related to google drive trying to back things up.
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + std_log_PAHs" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_2/S2_MG_bac_site_season_std_log_oil
	
The results made this model look overfit (CV increased after plateau).  Since the prior was already low, I decided to not use the log of PAHs even though the loss was lower.  But I suppose I should get Q2 values from QIIME2

##### Fungi

	songbird multinomial \
	--input-biom Soil2_ITS_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_2/S2_MG_fungi_null
	
	songbird multinomial \
	--input-biom Soil2_ITS_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_2/S2_MG_fungi_site_season_oil


##### Check results in qiime2 (v2020.8, on my laptop)

	qiime tools import \
		--input-path Soil2_ITS_100_songbird.biom \
		--output-path Soil2_ITS_100_songbird.biom.qza \
		--type FeatureTable[Frequency]
		
	qiime songbird multinomial \
	--i-table Soil2_ITS_100_songbird.biom.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/fung_site_season_oil_differentials.qza \
	--o-regression-stats q2/fung_site_season_oil_regression-stats.qza \
	--o-regression-biplot q2/fung_site_season_oil_regression-biplot.qza \
	
	qiime songbird multinomial \
	--i-table Soil2_ITS_100_songbird.biom.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER'))" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/fung_site_season_differentials.qza \
	--o-regression-stats q2/fung_site_season_regression-stats.qza \
	--o-regression-biplot q2/fung_site_season_regression-biplot.qza \
	
	qiime songbird multinomial \
	--i-table Soil2_ITS_100_songbird.biom.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "1" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/fung_null_differentials.qza \
	--o-regression-stats q2/fung_null_regression-stats.qza \
	--o-regression-biplot q2/fung_null_regression-biplot.qza \
	
##### Visualize the first model's regression stats *and* the null model's regression stats

site+season against null (as baseline)
		
	qiime songbird summarize-paired \
		--i-regression-stats q2/fung_site_season_regression-stats.qza \
		--i-baseline-stats q2/fung_null_regression-stats.qza \
		--o-visualization q2/fungi_site_season_paired-summary.qzv


	qiime tools view q2/fungi_site_season_paired-summary.qzv
	
site+season+oil against site+season (as baseline)

	qiime songbird summarize-paired \
		--i-regression-stats q2/fung_site_season_oil_regression-stats.qza \
		--i-baseline-stats q2/fung_site_season_regression-stats.qza \
		--o-visualization q2/fungi_oil_paired-summary.qzv

	qiime tools view q2/fungi_oil_paired-summary.qzv
	
site+season+oil against null (as baseline)

	qiime songbird summarize-paired \
		--i-regression-stats q2/fung_site_season_oil_regression-stats.qza \
		--i-baseline-stats q2/fung_null_regression-stats.qza \
		--o-visualization q2/fungi_oil_null_paired-summary.qzv

	qiime tools view q2/fungi_oil_null_paired-summary.qzv
	
### Bacteria
		
	qiime tools import \
		--input-path Soil2_16S_100_songbird.biom \
		--output-path Soil2_16S_100_songbird.biom.qza \
		--type FeatureTable[Frequency]
		
	qiime songbird multinomial \
	--i-table Soil2_16S_100_songbird.biom.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/bac_site_season_oil_differentials.qza \
	--o-regression-stats q2/bac_site_season_oil_regression-stats.qza \
	--o-regression-biplot q2/bac_site_season_oil_regression-biplot.qza \
	
	qiime songbird multinomial \
	--i-table Soil2_16S_100_songbird.biom.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER'))" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/bac_site_season_differentials.qza \
	--o-regression-stats q2/bac_site_season_regression-stats.qza \
	--o-regression-biplot q2/bac_site_season_regression-biplot.qza \
	
	qiime songbird multinomial \
	--i-table Soil2_16S_100_songbird.biom.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "1" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/bac_null_differentials.qza \
	--o-regression-stats q2/bac_null_regression-stats.qza \
	--o-regression-biplot q2/bac_null_regression-biplot.qza \
	
##### Visualize the first model's regression stats *and* the null model's regression stats

site+season against null (as baseline)

	qiime songbird summarize-paired \
	--i-regression-stats q2/bac_site_season_regression-stats.qza \
	--i-baseline-stats q2/bac_null_regression-stats.qza \
	--o-visualization q2/bac_site_season_paired-summary.qzv

	qiime tools view q2/bac_site_season_paired-summary.qzv

site+season+oil against site+season (as baseline)
	
	qiime songbird summarize-paired \
		--i-regression-stats q2/bac_site_season_oil_regression-stats.qza \
		--i-baseline-stats q2/bac_site_season_regression-stats.qza \
		--o-visualization q2/bac_oil_paired-summary.qzv

	qiime tools view q2/bac_oil_paired-summary.qzv

site+season+oil against null (as baseline)

	qiime songbird summarize-paired \
	--i-regression-stats q2/bac_site_season_oil_regression-stats.qza \
	--i-baseline-stats q2/bac_null_regression-stats.qza \
	--o-visualization q2/bac_oil_null_paired-summary.qzv
		
	qiime tools view q2/bac_oil_null_paired-summary.qzv


I notices that my results don't exactly line up with my past results, and I think it might have to do with the random seed.  I think I had set it to one in the notes I lost, but most recently I used the default which is zero, so I'm going to try and run it with 1 to see what happens.

##### Bacteria Standalone
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_null \
	--random-seed 1
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER'))" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season \
	--random-seed 1
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_oil \
	--random-seed 1

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + std_log_PAHs" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_std_log_oil \
	--random-seed 1
	
The results made this model look overfit (CV increased after plateau).  Since the prior was already low, I decided to not use the log of PAHs even though the loss was lower.  But I suppose I should get Q2 values from QIIME2

##### Fungi

	songbird multinomial \
	--input-biom Soil2_ITS_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "1" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/fungi/S2_MG_fungi_null \
	--random-seed 1
	
	songbird multinomial \
	--input-biom Soil2_ITS_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER'))" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/fungi/S2_MG_fungi_site_season \
	--random-seed 1
	
	songbird multinomial \
	--input-biom Soil2_ITS_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/fungi/S2_MG_fungi_site_season_oil \
	--random-seed 1


I figured it out.  The seed doesn't make a difference in results (which is good, that means the model isn't arbitrary).  But my values for scaled PAHs are different in my metadata files.  I must have used the data set without the two outliers for the "new_final" models.  I don't know if I did this on purpose, but I definitely desire to use the outliers now.  I think this because the scaled numbers are different when you don't use the outliers.  Also, I double checked log-transforming and then scaling the PAHs.  The numbers aren't exactly the same if I do it in excel or R, but they're pretty close.

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + std_naph" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_naph \
	--random-seed 1
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + std_phen" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_phen \
	--random-seed 1
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + std_diben" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_diben \
	--random-seed 1
	
Using individual PAH classes doens't improve the models more than Total PAHs.


Try using log ratios of 3-ring/Chrysene PAHs

Naphthalenes

	songbird multinomial \
		--input-biom Soil2_16S_100_songbird.biom \
		--metadata-file S2_metagenomics_sample_map_songbird.txt \
		--training-column Testing \
		--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + log_naph_chry" \
		--epochs 10000 \
		--differential-prior 0.5 \
		--summary-interval 1 \
		--num-random-test-examples 12 \
		--batch-size 12 \
		--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_log_naph_chry \
		--random-seed 1

Dibenzothiophenes
	
	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + log_dibenz_chry" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_log_dibenz_chry \
	--random-seed 1

Phenanthrenes

	songbird multinomial \
	--input-biom Soil2_16S_100_songbird.biom \
	--metadata-file S2_metagenomics_sample_map_songbird.txt \
	--training-column Testing \
	--formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + log_phen_chry" \
	--epochs 10000 \
	--differential-prior 0.5 \
	--summary-interval 1 \
	--num-random-test-examples 12 \
	--batch-size 12 \
	--summary-dir logdir/new_final_3/bac/S2_MG_bac_site_season_log_phen_chry \
	--random-seed 1
	
Visualize

		tensorboard --logdir . --host=127.0.0.1
		
Ratios don't improve the model fit.  In general oil doesn't improve fit over site and season.  How about if I use just dominant taxa?

Filter in qiime

	qiime feature-table filter-features \
	--i-table Soil2_16S_100_songbird.biom.qza \
	--m-metadata-file bac_dom_featureID.txt \
	--o-filtered-table S2_bac-filtered-table.qza
	
	qiime feature-table summarize --i-table S2_bac-filtered-table.qza --o-visualization S2_bac-filtered-table.qzv
  
Correctly kept 719 OTUs

	qiime songbird multinomial \
	--i-table S2_bac-filtered-table.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "1" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/filtered/bac_null_differentials.qza \
	--o-regression-stats q2/filtered/bac_null_regression-stats.qza \
	--o-regression-biplot q2/filtered/bac_null_regression-biplot.qza \
	--verbose \
	
	qiime songbird multinomial \
	--i-table S2_bac-filtered-table.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER')) + Total_relevant_PAHs_scaled" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/filtered/bac_site_season_oil_differentials.qza \
	--o-regression-stats q2/filtered/bac_site_season_oil_regression-stats.qza \
	--o-regression-biplot q2/filtered/bac_site_season_oil_regression-biplot.qza \
	--verbose \
	
	qiime songbird multinomial \
	--i-table S2_bac-filtered-table.qza \
	--m-metadata-file S2_metagenomics_sample_map_songbird.txt \
	--p-formula "C(site, Treatment('BJ')) + C(season, Treatment('WINTER'))" \
	--p-epochs 10000 \
	--p-differential-prior 0.5 \
	--p-training-column Testing \
	--p-summary-interval 1 \
	--o-differentials q2/filtered/bac_site_season_differentials.qza \
	--o-regression-stats q2/filtered/bac_site_season_regression-stats.qza \
	--o-regression-biplot q2/filtered/bac_site_season_regression-biplot.qza \
	--verbose \
	
##### Visualize

site+season+oil against null (as baseline)

	qiime songbird summarize-paired \
	--i-regression-stats q2/filtered/bac_site_season_oil_regression-stats.qza \
	--i-baseline-stats q2/filtered/bac_null_regression-stats.qza \
	--o-visualization q2/filtered/bac_oil_null_paired-summary.qzv
		
	qiime tools view q2/filtered/bac_oil_null_paired-summary.qzv
	
site+season against null (as baseline)

	qiime songbird summarize-paired \
	--i-regression-stats q2/filtered/bac_site_season_regression-stats.qza \
	--i-baseline-stats q2/filtered/bac_null_regression-stats.qza \
	--o-visualization q2/filtered/bac_SS_null_paired-summary.qzv
		
	qiime tools view q2/filtered/bac_SS_null_paired-summary.qzv