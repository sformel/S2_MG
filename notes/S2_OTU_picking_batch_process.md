#August 1, 2016
#Steve Formel

#Here I picked OTUs 3 different times to get a sense of how robust my data was with the OTU picking algorithm.  One set, Soil2, was picked against the 97_s version of UNITE.  The other set, Soil3, was picked against the dynamic_s version of UNITE.  I'm not sure that I understand the implications of picking against the dynamic version, so I think I will continue to use the version that is clustered at 97%.

#Ecept for OTU picking, all of it run in the login node of cypress

#first analyze the final sequences used in OTU picking
./software/FastQC/fastqc ./S2_July17/sequences/S2_trimmed.assembled.fastq --outdir=/lustre/project/svanbael/steve/S2_July17/sequences/FastQC_analysis
#sequence length distribution
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' S2_trimmed.assembled.fastq >> S2_trimmed.assembled.fastq_length_dist.txt 

#to compare to soil of ATGC seqs
#NR is the awk comman for number of records.  This tells it to count the length of every two lines starting at the second line.
awk 'NR%2 == 0 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' S2_ATGC_all_oiled_soil_samples.fasta >> S2_ATGC_all_oiled_soil_samples.fasta_length_dist.txt 

#something isn't right in the file that has all the smaples combined.  I'm going to start by trying to remove blank lines.
grep -c '^$' S2_ATGC_all_oiled_soil_samples.fasta #no blank lines
#remove white space but not new lines
cat *.fasta | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//' >> S2_ATGC_all_oiled_soil_samples.fasta  #didn't have enough space so I have to remake the file when I do it.
#above didn't make a difference in the weird things showing up in the ilne distirbution.
for i in *.fasta ; do awk 'NR%2 == 0 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' $i >> $i.txt; done 

#I've hunted down in file 98.merge.fasta there is a strange line that comes up as being "31501706" in length.
grep -nr '^3' filename
awk 'length>1000' 98.merge.fasta
grep -nr '^ATTGGGTATAAAGCGCACGTAGGCGGCGAGGCAAGTGTCGGGTGAAATCCCACAGCTCACCCGGGGAACTGCCCGGCAAACTGCA' 98.merge.fasta

#of course its the binary garbage at the end of the file
#to clean:
tr -cd '\11\12\15\40-\176' < 98.merge.fasta > 98.cleaned.merge.fasta
#did not work.
#delete the last line and add the sequence back in.
sed -i '$d' 98.merge.fasta >> 98.deleted_last_line.merge.fasta
#above wasn't outputting the file for some reason, though it appears to be the best solution.  An alternative:
head -n -1 98.merge.fasta >> 98.deleted_last_line.merge.fasta

#check that it's gone
awk 'length>1000' 98.deleted_last_line.merge.fasta
#add sequence back in
cat 98.deleted_last_line.merge.fasta 98.merge.fasta_lastline.txt >> 98.cleaned.merge.fasta
#one last check
awk 'length>1000' 98.cleaned.merge.fasta

#merge all files
cat

#rerun length distribution
awk 'NR%2 == 0 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' S2_ATGC_all_oiled_soil_samples.fasta >> S2_ATGC_all_oiled_soil_samples.fasta_length_dist.txt
#still geting funny results, going to try counting lines that don't start with >
awk '!/>/{print length}' S2_ATGC_all_oiled_soil_samples.fasta | sort | uniq -c >> S2_ATGC_all_oiled_soil_samples_length_dist.txt
#worked!
#to compare apples to apples I'm oing to do the other fasta file the same way
awk '!/>/{print length}' S2_trimmed_assembled.fasta | sort | uniq -c >> S2_trimmed_assembled.fasta_length_dist.txt

-------------------------------------------
S1. picked OTUs against combined Nov 2016 UNITE 97_s (Soil2) or dynamic_s (Soil3) combined with GreenGenes 13_8

#script looked something like this:

	pick_open_reference_otus.py -i /lustre/project/svanbael/steve/S2_July17/sequences/S2_trimmed_assembled.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2 -a -O 200 -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta --suppress_step4 --suppress_taxonomy_assignment --suppress_align_and_tree



S2.	ASSIGN TAXONOMY

#with the UNITE developed taxonomy which is untrimmed - https://groups.google.com/forum/#!topic/qiime-forum/cvuJzDZF9AI
		#Note: the dynamic UNITE clusters at 98.5, not 97% (or maybe multiple levels?)

	#Soil 2: not the developer version
	assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2/Step2

	assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step2

	assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step2

	#Soil 3: the developer version

	assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil3/Step2

	assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step2

	assign_taxonomy.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/rep_set.fna -t /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_tax_db.txt -r /lustre/project/svanbael/steve/S2_July17/databases/combined_db/GG_UNITE_db_20Nov2016.fasta -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step2
		
------------------------------------------------------------------------------------
S3.	CHIMERA FILTERING - Note this variation is due to variation in the OTU picking process. VSEARCH seems to come up with the same results every time for chimera detection.

#Soil2

mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/Step3

vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil2/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil2/Step3/Soil2_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil2/Step3/Soil2_non_chimeras.fna --strand plus

	#2284 chimeras
	
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step3

vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil2_1/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step3/Soil2_1_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step3/Soil2_1_non_chimeras.fna --strand plus

	#2240 chimeras
	
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step3

vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil2_2/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step3/Soil2_2_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step3/Soil2_2_non_chimeras.fna --strand plus

	#2228 chimeras
	
#Soil3

mkdir /lustre/project/svanbael/steve/S2_July17/Soil3/Step3

vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil3/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil3/Step3/Soil3_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil3/Step3/Soil3_non_chimeras.fna --strand plus

	# 2281 chimeras
	
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step3

vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil3_1/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step3/Soil3_1_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step3/Soil3_1_non_chimeras.fna --strand plus

	# 2232 chimeras
	
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step3

vsearch --uchime_ref /lustre/project/svanbael/steve/S2_July17/Soil3_2/rep_set.fna --db /lustre/project/svanbael/steve/S2_July17/databases/chimera_db/uchime_GG_chimera_db.fasta --chimeras /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step3/Soil3_2_chimeras.fna --nonchimeras /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step3/Soil3_2_non_chimeras.fna --strand plus

	# 2265 chimeras
	
---------------------------------------------------------------------------
S4.	Make OTU Tables

#Soil2	

mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables

#An OTU has an e with a dierisis (umlaut) so to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/ASCII_find_Soil2.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments_clean.txt
	
make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil2/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil2/Step3/Soil2_chimeras.fna

mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables

#An OTU has an e with a dierisis (umlaut) so to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step2/ASCII_find_Soil2_1.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step2/rep_set_tax_assignments_clean.txt
	
make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil2_1/Step3/Soil2_1_chimeras.fna

mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables

#An OTU has an e with a dierisis (umlaut) so to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step2/ASCII_find_Soil2_2.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step2/rep_set_tax_assignments_clean.txt
	
make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil2_2/Step3/Soil2_2_chimeras.fna

#Soil3

mkdir /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables

#An OTU has an e with a dierisis (umlaut) so to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil3/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil3/Step2/ASCII_find_Soil3.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil3/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil3/Step2/rep_set_tax_assignments_clean.txt
	
make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil3/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil3/Step3/Soil3_chimeras.fna

mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables

#An OTU has an e with a dierisis (umlaut) so to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step2/ASCII_find_Soil3_1.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step2/rep_set_tax_assignments_clean.txt
	
make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil3_1/Step3/Soil3_1_chimeras.fna

mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables

#An OTU has an e with a dierisis (umlaut) so to find and replace it I used PERL and these commands in cypress:
	
	module load pcre/8.39
	grep --color='auto' -P -n '[^\x00-\x7F]' /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step2/rep_set_tax_assignments.txt >> /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step2/ASCII_find_Soil3_2.txt
	sed 's/ë/e/g' /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step2/rep_set_tax_assignments.txt > /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step2/rep_set_tax_assignments_clean.txt
	
make_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/final_otu_map_mc2.txt -t /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step2/rep_set_tax_assignments_clean.txt -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_orig.biom -e /lustre/project/svanbael/steve/S2_July17/Soil3_2/Step3/Soil3_2_chimeras.fna

---------------------------------------------------------------------------	
S5. #Break OTU Table in tables of: Bacteria, Fungi, Archaea, and "Unassigned"

#Soil2	
#remove mitochondria & chloroplast
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

#break out tables by taxa
	
#Unassigned OTUs
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass.biom -p Unassigned
	
#16S (bacteria and archaea)	
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S.biom -p k__bacteria,k__archaea  
	
#fungi
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS.biom -p k__fungi	

#Soil2_1	

#remove mitochondria & chloroplast
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

#break out tables by taxa

#Unassigned OTUs
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_unass.biom -p Unassigned
	
#16S (bacteria and archaea)	
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_16S.biom -p k__bacteria,k__archaea  
	
#fungi
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_ITS.biom -p k__fungi	

#Soil2_2	

#remove mitochondria & chloroplast
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

#break out tables by taxa

#Unassigned OTUs
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_unass.biom -p Unassigned
	
#16S (bacteria and archaea)	
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_16S.biom -p k__bacteria,k__archaea  
	
#fungi
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_ITS.biom -p k__fungi


#Soil3	
#remove mitochondria & chloroplast
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

#break out tables by taxa
	
#Unassigned OTUs
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_unass.biom -p Unassigned
	
#16S (bacteria and archaea)	
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_16S.biom -p k__bacteria,k__archaea  
	
#fungi
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_ITS.biom -p k__fungi	

#Soil3_1	

#remove mitochondria & chloroplast
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

#break out tables by taxa

#Unassigned OTUs
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_unass.biom -p Unassigned
	
#16S (bacteria and archaea)	
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_16S.biom -p k__bacteria,k__archaea  
	
#fungi
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_ITS.biom -p k__fungi	

#Soil3_2	

#remove mitochondria & chloroplast
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_orig.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_no_chloro_no_mito.biom -n f__mitochondria,c__Chloroplast

#break out tables by taxa

#Unassigned OTUs
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_unass.biom -p Unassigned
	
#16S (bacteria and archaea)	
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_16S.biom -p k__bacteria,k__archaea  
	
#fungi
filter_taxa_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_no_chloro_no_mito.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_ITS.biom -p k__fungi
---------------------------------------------------------------------------
S6. #Make summaries of OTU Tables

#Soil2
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_unass.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_16S.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_ITS.txt

#Soil2_1
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries/Soil2_1_unass.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries/Soil2_1_16S.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries/Soil2_1_ITS.txt

#Soil2_2
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries/Soil2_2_unass.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries/Soil2_2_16S.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries/Soil2_2_ITS.txt

#Soil3
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries/Soil3_unass.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries/Soil3_16S.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries/Soil3_ITS.txt

#Soil3_1
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries/Soil3_1_unass.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries/Soil3_1_16S.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries/Soil3_1_ITS.txt

#Soil3_2
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_unass.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries/Soil3_2_unass.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_16S.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries/Soil3_2_16S.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_ITS.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries/Soil3_2_ITS.txt


#FILTER OUT SEQUENCES WITH FEWER THAN 100 SEQS
#Soil2
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS_100.biom

#Soil2_1
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_unass_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_16S_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_ITS_100.biom

#Soil2_2
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_unass_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_16S_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_ITS_100.biom

#Soil3
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_unass_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_16S_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_ITS_100.biom

#Soil3_1
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_unass_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_16S_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_ITS_100.biom

#Soil3_2
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_unass.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_unass_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_16S.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_16S_100.biom
filter_otus_from_otu_table.py -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_ITS.biom -n 101 -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_ITS_100.biom

#Make summaries of filtered OTU Tables

#Soil2
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_unass_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_unass_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_16S_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_16S_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2/biom_tables/Soil2_ITS_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2/biom_summaries/Soil2_ITS_100.txt

#Soil2_1
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_unass_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries/Soil2_1_unass_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_16S_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries/Soil2_1_16S_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_tables/Soil2_1_ITS_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_1/biom_summaries/Soil2_1_ITS_100.txt

#Soil2_2
mkdir /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_unass_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries/Soil2_2_unass_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_16S_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries/Soil2_2_16S_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_tables/Soil2_2_ITS_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil2_2/biom_summaries/Soil2_2_ITS_100.txt

#Soil3
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_unass_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries/Soil3_unass_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_16S_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries/Soil3_16S_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3/biom_tables/Soil3_ITS_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3/biom_summaries/Soil3_ITS_100.txt

#Soil3_1
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_unass_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries/Soil3_1_unass_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_16S_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries/Soil3_1_16S_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_tables/Soil3_1_ITS_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_1/biom_summaries/Soil3_1_ITS_100.txt

#Soil3_2
mkdir /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_unass_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries/Soil3_2_unass_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_16S_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries/Soil3_2_16S_100.txt
biom summarize-table -i /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_tables/Soil3_2_ITS_100.biom -o /lustre/project/svanbael/steve/S2_July17/Soil3_2/biom_summaries/Soil3_2_ITS_100.txt

#not done as of 11pm on 29 July 17
-------------------------------------------------------------------------
S7. Align Sequences - using QIIME align_seqs.py (GreenGenes template) - run as script "S2_Step7_Alignment.sh"  RUNTIME = < 10 minutes

	align_seqs.py -i ./S2_Step3/S2_both_db_nonchimeras.fna -o ./S2_Step7_Alignment
------------------------------------------------------------------------
S8.	Filter alignment - using QIIME filter_alignment.py
	
		filter_alignment.py -i ./S2_Step7_Alignment/S2_both_db_nonchimeras_aligned.fasta -o ./S2_Step8/filtered_alignment/
-----------------------------------------------------------------------------
S9. Build tree for OTUs that passed alignment - QIIME (fasttree) make_phylogeny.py

	#Run as script "S2_Step9.sh" RUNTIME = < 10 minutes
	
		make_phylogeny.py -i ./S2_Step8/filtered_alignment/S2_both_db_nonchimeras_aligned_pfiltered.fasta -r midpoint -o ./S2_Step9/S2_tree.tre
		
-------------------------------------------------------------------------
S10. BLAST Unassigned OTUs against NCBI nt database

		MAKE NCBI BLAST DATABASE TO EXPLORE UNASSIGNED OTUS

		#download files.  Need 34G of spaces
		wget ftp://ftp.ncbi.nih.gov/blast/db/nt.*tar.gz
		wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz #This is the taxid for all entries		

		#nt.*tar.gz | nucleotide sequence database, with entries from all traditional divisions of GenBank, EMBL, and DDBJ excluding bulk divisions (gss, sts, pat, est, and htg divisions. wgs entries are also excluded. Not non-redundant.  Downloaded 22 May 2017.

		for i in *.tar.gz ; do tar -zxvf $i ; done  #extract files, tar can't batch extract, you need to do it in a loop.  As long as the files are in the same directory, the databse should already be "made", you don't need to run the ocmmand "makeblastdb"

		#Also add the blast db to your path
		nano ~/.bash_profile

		export BLASTDB=$BLASTDB:/media/Data/path/to/your/database/files/

#pick representative sets to blast (i.e. a sequence for each OTU)

	#filter fasta using biom tables - run in login node
	
		filter_fasta.py -f ./ASAS_open_otus_both_db/rep_set.fna -b ./biom_tables/soil_only/SOIL_ONLY_open_ref_no4_unassigned_nochimera_otu_table.biom -o ./S2_Step10/unassigned_OTUs_rep_set.fna

		#BLAST - run as script "S2_Step10_BLAST_unassigned.sh" RUNTIME = 

			#it looks like you have to be explicit in your filepaths
			#describes output format options: https://www.ncbi.nlm.nih.gov/books/NBK279675/

			blastn -db /lustre/project/svanbael/steve/S2/db/ncbi_nt/nt -query /lustre/project/svanbael/steve/S2/S2_Step10/unassigned_rep_set.fna -outfmt '7 qseqid sseqid sskingdoms length pident qstart qend sstart send evalue bitscore sgi sacc staxids sscinames' -num_threads 10 -out /lustre/project/svanbael/steve/S2/S2_Step10/unassigned_blast_results.txt
			
			#For some reason it did not include taxonomy although the script is the same as when it previously included taxonomy.
			
		#filter blast output by column
			#filter output  $X is column number X

					awk '$4>100' unassigned_blast_results.txt >> filtered_blast_align_100.txt
					awk '$5>90' filtered_blast_align_100.txt >> filtered_blast_align_100_ident_90.txt
					