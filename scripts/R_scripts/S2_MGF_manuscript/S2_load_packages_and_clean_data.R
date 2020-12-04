#S2_MG Project

#Main script to load packages, make phyloseq objects, and subset objects by site and season

#by Steve Formel
#Last updated July 15, 2020

#This is analysis for the S2_2season_metagenomics project mainly done in the package "phyloseq" with the help of many other packages.
#Biom tables were created in QIIME (1.9.1), script is "S2_Master_V3_Scripts_Notes_25May17.txt"

#Note: updated on Dec 5, 2018 by SF to remove OTUs with < 100 total seqs after subsetting to the 2 season data set.

# Load Packages and Notes -------------------------------------------------

#load some packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(data.table)
library(tidyr)
library(genefilter)

#import biom and variable map--------------------------------------

#Import Bioms

bac.biom <- import_biom("data/QIIME_output/Soil2/biom_tables/Soil2_16S_100.biom", parseFunction=parse_taxonomy_greengenes)
fung.biom <- import_biom("data/QIIME_output/Soil2/biom_tables/Soil2_ITS_100.biom", parseFunction=parse_taxonomy_greengenes)

#Gives warning: 

# Warning message:
# In strsplit(conditionMessage(e), "\n") :
#   input string 1 is invalid in this locale

#As far as I can tell, this isn't critical. Names seem to match.  Oddly the parse_taxonomy_greengenes function adds several additional empty ranks, but I delete those below.

#See: https://github.com/joey711/phyloseq/issues/986


#reorder sample map to match order of OTU table - only need to do the first time, it is now fixed, code is left in as a record

  #A <- read.csv("S2_metagenomics_sample_map_30Oct17.csv", row.names = 1)
  #B <- A[(sample_names(bac.biom)),]
  #write.csv(B, "S2_metagenomics_sample_map_ordered_to_match_OTU_table_31Oct17.csv")

#import map file
S2_MAP <-import_qiime_sample_data("data/Environmental_data/S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt")

#Make interaction term for site and season
S2_MAP$site_season <- interaction(S2_MAP$site, S2_MAP$season)


#Merge bac biom, and map into one phyloseq object
bac.all <- merge_phyloseq(bac.biom,S2_MAP)

#Merge fung biom and map into one phyloseq object
fung.all <- merge_phyloseq(fung.biom,S2_MAP)

#Fix pecadillos---------------------------------

#Make phyloseq objects a list so I can use lapply for cleaning, this keeps the code from being excessive.

psq.list <- c(bac.all, fung.all)

#make genotype a factor

psq.list <- lapply(psq.list, function(x) {sample_data(x)$genotype <- factor(sample_data(x)$genotype) ; x})

#arrange date in order

psq.list <- lapply(psq.list, function(x) {sample_data(x)$date_collected <- factor(sample_data(x)$date_collected, levels = c( "8/9/2012", "1/18/2013", "1/19/2013", "6/26/2013", "7/10/2013", "11/15/2013", "11/20/2013", "12/3/2013")) ; x})

#arrange seasons

psq.list <- lapply(psq.list, function(x) {sample_data(x)$season <- factor(sample_data(x)$season, levels = c("WINTER", "SUMMER", "FALL")) ; x})

#the parse_taxonomy function in biom import adds ranks to tax_table, but oddly only for 16S remove them here.

colnames(tax_table(bac.all))
any(!is.na(tax_table(bac.all)[,8:13]))

colnames(tax_table(fung.all)) #hmm that's fine

#All cells ar NA, remove empty columns
tax_table(bac.all) <- tax_table(bac.all)[,1:7]

#any OTUs less than 100 total or empty OTUs?  There shouldn't be.
any(taxa_sums(bac.all) <= 100)  
any(taxa_sums(fung.all) <= 100)  

#Subset to phyloseq objects of all dates by SITE--------

bac.all.BJ <- bac.all %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.all.F <- bac.all %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.BJ <- fung.all %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.F <- fung.all %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to phyloseq objects of all dates by WINTER-------

bac.all.BJ.winter <- bac.all %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.all.F.winter <- bac.all %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.BJ.winter <- fung.all %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.F.winter <- fung.all %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to phyloseq objects of all dates by SUMMER-------

bac.all.BJ.summer <- bac.all %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.all.F.summer <- bac.all %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.BJ.summer <- fung.all %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.F.summer <- fung.all %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to Jan 2013 and Summer (June July) 2013 sampling points only----------------

#includes outliers

bac.2season_with_outliers <- subset_samples(bac.all, date_collected%in%c("1/18/2013", "1/19/2013", "6/26/2013", "7/10/2013"))
fung.2season_with_outliers <- subset_samples(fung.all, date_collected%in%c("1/18/2013", "1/19/2013", "6/26/2013", "7/10/2013"))

#probably a good idea to avoid these samples.  36 is suspiciously depauperate, 23 is not, but I know that oil is magnitudes higher in each sample.

bac.2season <- prune_samples(!sample_data(bac.2season_with_outliers)$SampleID %in% c("S2.36", "S2.23"), bac.2season_with_outliers)

#So I can test the effects of filtering
bac.2season.not_filtered <- bac.2season

fung.2season <- prune_samples(!sample_data(fung.2season_with_outliers)$SampleID %in% c("S2.36", "S2.23"), fung.2season_with_outliers)

#So I can test the effects of filtering
fung.2season.not_filtered <- fung.2season

#Filter taxa that don't have at least 100 sequences in 5 samples-----------

psq.list.2season <- c("bac.2season"= bac.2season, "fung.2season" = fung.2season)

#A = The value you want to exceed.
#k = The number of elements that have to exceed A

for (i in 1:2){
  voi <- psq.list.2season[[i]]
  flist_a <- filterfun(kOverA(A = 0.001*mean(sample_sums(voi)), k = 0.05*(length(sample_sums(voi)))))
  a <- filter_taxa(voi, flist_a)
  psq.list.2season[[i]] <- prune_taxa(a, voi)
}

bac.2season <- psq.list.2season[[1]]
fung.2season <- psq.list.2season[[2]]

# The below was commented out to show that rare species strongly affect patterns in the relationship between diversity and oil

#now filter with outliers
psq.list.2season_with_outliers <- c("bac.2season_with_outliers"= bac.2season_with_outliers, "fung.2season_with_outliers" = fung.2season_with_outliers)

for (i in 1:2){
  voi <- psq.list.2season_with_outliers[[i]]
  flist_a <- filterfun(kOverA(A = 0.001*mean(sample_sums(voi)), k = 0.05*(length(sample_sums(voi)))))
  a <- filter_taxa(voi, flist_a)
  psq.list.2season_with_outliers[[i]] <- prune_taxa(a, voi)
}

bac.2season_with_outliers_filtered <- psq.list.2season_with_outliers[[1]]
fung.2season_with_outliers_filtered <- psq.list.2season_with_outliers[[2]]

# Subset "2season" by Site ------------

bac.2season.BJ <- bac.2season %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F <- bac.2season %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ <- fung.2season %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F <- fung.2season %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

bac.2season.BJ.winter <- bac.2season %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.winter <- bac.2season %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.winter <- fung.2season %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter <- fung.2season %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to "2season" by SUMMER-------

bac.2season.BJ.summer <- bac.2season %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.summer <- bac.2season %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.summer <- fung.2season %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer <- fung.2season %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset with outliers------

# Subset "2season" by Site ------------

bac.2season.BJ_with_outliers <- bac.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F_with_outliers <- bac.2season_with_outliers %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

bac.2season.BJ.winter_with_outliers <- bac.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.winter_with_outliers <- bac.2season_with_outliers %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.winter_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to "2season" by SUMMER-------

bac.2season.BJ.summer_with_outliers <- bac.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.summer_with_outliers <- bac.2season_with_outliers %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.summer_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#FILTERED SUBSET------

#Same but with filtered species------

# Subset "2season" by Site using filtered taxa------------

bac.2season.BJ_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

bac.2season.BJ.winter_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.winter_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.winter_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to "2season" by SUMMER-------

bac.2season.BJ.summer_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.summer_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.summer_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#rarefy and subset phyloseq objects for the 2season----------------

set.seed(28132)

bac.2season_even = rarefy_even_depth(bac.2season, sample.size = 17149)
bac.2season_even_with_outliers = rarefy_even_depth(bac.2season_with_outliers, sample.size = 17149)
bac.2season_even_not_filtered = rarefy_even_depth(bac.2season.not_filtered, sample.size = 17149)

fung.2season_even = rarefy_even_depth(fung.2season, sample.size = 1041)
fung.2season_even_with_outliers = rarefy_even_depth(fung.2season_with_outliers, sample.size = 1041)
fung.2season_even_not_filtered = rarefy_even_depth(fung.2season.not_filtered, sample.size = 1041)

