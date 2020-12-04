#S2_MGF Project

#Main script to load packages, make phyloseq objects, and subset objects by site and season

#by Steve Formel
#Last updated Dec 3, 2020

#Biom tables were created in QIIME (1.9.1), script is "S2_Master_V3_Scripts_Notes_25May17.txt"

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

fung.biom <- import_biom("~/Google Drive/VB_lab/VBL_data/Collections/S2_MG collection/data/QIIME_output/Soil2/biom_tables/Soil2_ITS_100.biom", parseFunction=parse_taxonomy_greengenes)

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
S2_MAP <-import_qiime_sample_data("~/Google Drive/VB_lab/VBL_data/Collections/S2_MG collection/data/Environmental_data/S2_metagenomics_sample_map_ordered_to_match_OTU_table_QIIME_31Oct17.txt")

#Make interaction term for site and season
S2_MAP$site_season <- interaction(S2_MAP$site, S2_MAP$season)

#Merge fung biom and map into one phyloseq object
fung.all <- merge_phyloseq(fung.biom,S2_MAP)

#Fix pecadillos---------------------------------

#Make phyloseq objects a list so I can use lapply for cleaning, this keeps the code from being excessive.

#make genotype a factor

sample_data(fung.all)$genotype <- factor(sample_data(fung.all)$genotype)

#arrange date in order

sample_data(fung.all)$date_collected <- factor(sample_data(fung.all)$date_collected, 
                                               levels = c( "8/9/2012", 
                                                           "1/18/2013", 
                                                           "1/19/2013", 
                                                           "6/26/2013", 
                                                           "7/10/2013", 
                                                           "11/15/2013", 
                                                           "11/20/2013", 
                                                           "12/3/2013"))

#arrange seasons

sample_data(fung.all)$season <- factor(sample_data(fung.all)$season, 
                                       levels = c("WINTER", 
                                                  "SUMMER", 
                                                  "FALL"))

#any OTUs less than 100 total or empty OTUs?  There shouldn't be.
any(taxa_sums(fung.all) <= 100)  

#Subset to phyloseq objects of all dates by SITE--------

fung.all.BJ <- fung.all %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.F <- fung.all %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to phyloseq objects of all dates by WINTER-------

fung.all.BJ.winter <- fung.all %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.F.winter <- fung.all %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to phyloseq objects of all dates by SUMMER-------

fung.all.BJ.summer <- fung.all %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.all.F.summer <- fung.all %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to Jan 2013 and Summer (June July) 2013 sampling points only----------------

#includes outliers
fung.2season_with_outliers <- subset_samples(fung.all, date_collected%in%c("1/18/2013", "1/19/2013", "6/26/2013", "7/10/2013"))

#probably a good idea to avoid these samples.  36 is suspiciously depauperate, 23 is not, but I know that oil is magnitudes higher in each sample.

#So I can test the effects of filtering
fung.2season <- prune_samples(!sample_data(fung.2season_with_outliers)$SampleID %in% c("S2.36", "S2.23"), fung.2season_with_outliers)

fung.2season.not_filtered <- fung.2season

#Filter taxa that don't have at least 100 sequences in 5 samples-----------

#A = The value you want to exceed.
#k = The number of elements that have to exceed A

#Commented out because our conclusiosn were the same whether rare species were filtered out or not.

# voi <- fung.2season
# flist_a <- filterfun(kOverA(A = 0.001*mean(sample_sums(voi)), k = 0.05*(length(sample_sums(voi)))))
# a <- filter_taxa(voi, flist_a)
# fung.2season <- prune_taxa(a, voi)

#Filter with outliers-----

voi <- fung.2season_with_outliers
flist_a <- filterfun(kOverA(A = 0.001*mean(sample_sums(voi)), k = 0.05*(length(sample_sums(voi)))))
a <- filter_taxa(voi, flist_a)
fung.2season_with_outliers_filtered <- prune_taxa(a, voi)

# Subset "2season" by Site ------------

fung.2season.BJ <- fung.2season %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F <- fung.2season %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

fung.2season.BJ.winter <- fung.2season %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter <- fung.2season %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by SUMMER-------

fung.2season.BJ.summer <- fung.2season %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer <- fung.2season %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

# Subset "2season" by Site ------------

fung.2season.BJ_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

fung.2season.BJ.winter_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by SUMMER-------

fung.2season.BJ.summer_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer_with_outliers <- fung.2season_with_outliers %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#FILTERED SUBSET------

#Same but with filtered species

fung.2season.BJ_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

fung.2season.BJ.winter_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to "2season" by SUMMER-------

fung.2season.BJ.summer_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#rarefy and subset phyloseq objects for the 2season----------------

set.seed(28132)

fung.2season_even = rarefy_even_depth(fung.2season, sample.size = 1041)
fung.2season_even_with_outliers = rarefy_even_depth(fung.2season_with_outliers, sample.size = 1041)
fung.2season_even_not_filtered = rarefy_even_depth(fung.2season.not_filtered, sample.size = 1041)

#Clear unwanted Objects

rm(a)
rm(flist_a)
rm(voi)
