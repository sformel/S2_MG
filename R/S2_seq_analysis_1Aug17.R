#Script for making histogram of fastq sequence length distribution.
#1 Aug 2017
#by Steve Formel

#Run in command line, not in R
  
  #Make distribution with the Unix command:
      #awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' all_R1.allsamples.fastq >> all.R1.fastq_length_dist.txt

# or for fasta
#awk '!/>/{print length}' S2_trimmed_assembled.fasta | sort | uniq -c >> S2_trimmed_assembled.fasta_length_dist.txt

#Above took ~ 10 min to run on the cypress login node.

setwd("C:/Users/Stephen and Moppa/Google Drive/Van Bael Google Drive/Boggs443 Data & User folders/Users/Grad Students/Steve/S2_21July2017/Seq_analysis_1Aug17/")

#read in data
ATGC_dist <- read.table("S2_ATGC_all_oiled_soil_samples_length_dist.txt", col.names = c("Count","Seq_Length"), colClasses = c("numeric", "numeric"))
SF_dist <- read.table("S2_trimmed.assembled.fastq_length_dist.txt", col.names = c("Seq_Length", "Count"))

#sort by seq_length
ATGC.ordered <- ATGC_dist[order(ATGC_dist$Seq_Length),]
SF.ordered <- SF_dist[order(SF_dist$Seq_Length),]

#plot, not technically a histogram, but serves as a freq. dist all the same
library(ggplot2)

# Colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(ATGC.ordered, aes(Seq_Length, Count, fill = Seq_Length))
p + geom_bar(stat = "identity") + scale_y_log10() +ylab("log10(Count)") + labs(title = "S2 Fasta Read Length Distribution for ATGC final reads")
ggsave("ATGC_seq_length_dist.png", width = 10, height = 8, device = "png")

p <- ggplot(SF.ordered, aes(Seq_Length, Count, fill = Seq_Length))
p + geom_bar(stat = "identity") + scale_y_log10() +ylab("log10(Count)") + labs(title = "S2 FastQ Read Length Distribution for SF final reads")
ggsave("SF_seq_length_dist.png", width = 10, height = 8, device = "png")

