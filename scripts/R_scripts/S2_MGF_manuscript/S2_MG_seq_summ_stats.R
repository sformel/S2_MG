#S2_MG Sequencing Summary Stats

#Last updated Sep 15, 2020

#by Steve Formel

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)

#Summarize Bacteria

sample_sums(bac.2season_with_outliers) %>%
  sum() 

sample_sums(bac.2season_with_outliers) %>%
  mean()

sample_sums(bac.2season_with_outliers) %>%
  median()

sample_sums(fung.2season_with_outliers) %>%
  sum() 

sample_sums(fung.2season_with_outliers) %>%
  mean()

sample_sums(fung.2season_with_outliers) %>%
  median()

#OTU Histogram-------

#How many OTUs are in 1,2,3..89 samples?

library(genefilter)

#A = The value of seqs you want to exceed.
#k = The number of samples that have to exceed A

voi <- bac.2season_with_outliers
  
bac.hist.df <- data.frame("n_samples" = c(1:89), "n_otus" = NA)

for (j in 1:length(sample_sums(voi))){
    flist_a <- filterfun(kOverA(A = 0, k = j))
    bac.hist.df[j,2] <- sum(filter_taxa(voi, flist_a))
    } #end j

ggplot(data = bac.hist.df,
       aes(x = n_samples,
           y = n_otus)) +
  geom_col()

voi <- fung.2season_with_outliers

fung.hist.df <- data.frame("n_samples" = c(1:89), "n_otus" = NA)

for (j in 1:length(sample_sums(voi))){
  flist_a <- filterfun(kOverA(A = 0, k = j))
  fung.hist.df[j,2] <- sum(filter_taxa(voi, flist_a))
} #end j

ggplot(data = fung.hist.df,
       aes(x = n_samples,
           y = n_otus)) +
  geom_col()

#What are the bacteria that are common to all samples?

flist_a <- filterfun(kOverA(A = 0, k = 89))
a <- filter_taxa(bac.2season_with_outliers, flist_a)

bac.89 <- prune_taxa(a, bac.2season_with_outliers)

taxa_sums(bac.89)
tax_table(bac.89)

tax.means <- apply(X = otu_table(bac.89), MARGIN = 1, FUN = mean)
tax.sd <- apply(X = otu_table(bac.89), MARGIN = 1, FUN = sd)

hist(unlist(tax.means))
sort(tax.means)

hist(unlist(tax.sd))
sort(tax.sd)

plot(tax.means, tax.sd)

tax_table(bac.89)[rownames(tax_table(bac.89)) %in% c("592957", "589224", "1111247", "790834", "320758"),]

otu_table(bac.89)[rownames(otu_table(bac.89)) %in% c("592957", "589224", "1111247", "790834", "320758"),]

#The above doesn't say much, nothing is that even across all samples, some of it is probably cross contamination from sequencing and other library prep processes.

