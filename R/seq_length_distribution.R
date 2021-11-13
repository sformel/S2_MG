#Script for making histogram of fastq sequence length distribution.
#Last updated Nov 8 2021
#by Steve Formel

#Run in command line, not in R
  
#Make distribution with the Unix command:

#awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' all_R1.allsamples.fastq >> all.R1.fastq_length_dist.txt

#Above took ~ 10 min to run on the cypress login node.

#read in data
r1 <- read.table("../data/Seq_analysis_output/all.R1.fastq_length_dist.txt", 
                 col.names = c("Seq_Length", "Count"))

r2 <- read.table("../data/Seq_analysis_output/all.R2.fastq_length_dist.txt", 
                 col.names = c("Seq_Length", "Count"))

#sort by seq_length
r1.ordered <- r1[order(r1$Seq_Length),]
r2.ordered <- r2[order(r2$Seq_Length),]

#plot, not technically a histogram, but serves as a freq. dist all the same
library(ggplot2)

# Colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(r1.ordered, aes(Seq_Length, Count, fill = Seq_Length))
p + geom_bar(stat = "identity") + scale_y_log10() +ylab("log10(Count)") + labs(title = "S2 FastQ Read Length Distribution for all R1 reads")

p <- ggplot(r2.ordered, aes(Seq_Length, Count, fill = Seq_Length))
p + geom_bar(stat = "identity") + scale_y_log10() +ylab("log10(Count)") + labs(title = "S2 FastQ Read Length Distribution for all R2 reads")

