#Ranked Differential plots

#by Steve Formel
#Last updated 11 Oct, 2020
#Make ranked abundance plot for songbird analysis

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(ggforce)

#Bacteria-----

#load Songbird output-----

sb <- read.delim("songbird/logdir/new_final_3/bac/S2_MG_bac_site_season/differentials.tsv")

#rename columns
names(sb) <- c("featureid", "Intercept", "Fourchon_rel_BJ", "Summer_rel_Wint")

#Using not filtered because the biom used in songbird wasn't filtered this way 
bac.2season.sub <- subset_taxa(bac.2season_with_outliers, rownames(tax_table(bac.2season_with_outliers)) %in% sb$featureid)
 
sb.filtered <- sb[sb$featureid %in% rownames(tax_table(bac.2season_with_outliers_filtered)), ]

#bind taxonomy with songbird results
sb <- cbind(sb, tax_table(bac.2season.sub))

#Factor for whether something is dominant or not

sb$rare_dom <- factor(ifelse(sb$featureid %in% sb.filtered$featureid,
                      "dominant",
                      "rare"))

#Replace NA with "unidentified"
sb[,6:12] <- lapply(sb[,6:12], replace_na, replace = "unidentified")

#I'm basically imitating the QURRO plots: https://www.biorxiv.org/content/10.1101/2019.12.17.880047v1.full

#rename in case I need to come back
sb.mod <- sb

#Add rank columns for x axis-----
sb.mod$rank_site <- rank(sb.mod$Fourchon_rel_BJ)
sb.mod$rank_season <- rank(sb.mod$Summer_rel_Wint)

#Add percentile columns 
sb.mod$percentile_site <- cut(sb.mod$rank_site , 
                         breaks = quantile(sb.mod$rank_site, 
                                           seq(from = 0, to = 1, by = 0.1)), 
                         labels=1:10, 
                         include.lowest=TRUE)

sb.mod$percentile_season <- cut(sb.mod$rank_season , 
                         breaks = quantile(sb.mod$rank_season, 
                                           seq(from = 0, to = 1, by = 0.1)), 
                         labels=1:10, 
                         include.lowest=TRUE)


#plot taxonomic ranks as pdf for browsing-----

#Site-----

#How many orders are there?
length(unique(sb.mod$Order)) #222
#If I divide by 6, I get 37 pages

pdf("images/manuscript/fig7new/S2_MG_explore_bac_order_rankings_by_site.pdf")
for (i in 1:37){
  print(ggplot(data = sb.mod,
               aes(x = rank_site,
                   y = Fourchon_rel_BJ,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#site by family

#How many families are there?
length(unique(sb.mod$Family)) #233
#If I divide by 6, I get 39 page, but then I need to somehow account for sharing orders.  I'm going to guess at 60

pdf("images/manuscript/fig7new/S2_MG_explore_bac_family_rankings_by_site.pdf")
for (i in 1:60){
  print(ggplot(data = sb.mod,
               aes(x = rank_site,
                   y = Fourchon_rel_BJ,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels

pdf("images/manuscript/fig7new/S2_MG_explore_bac_phylum_to_genus_rankings_by_site.pdf")
for (i in 1:100){
  print(ggplot(data = sb.mod,
               aes(x = rank_site,
                   y = Fourchon_rel_BJ,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Season-----

#How many orders are there?
length(unique(sb.mod$Order)) #222
#If I divide by 6, I get 37 pages

pdf("images/manuscript/fig7new/S2_MG_explore_bac_order_rankings_by_season.pdf")
for (i in 1:37){
  print(ggplot(data = sb.mod,
               aes(x = rank_season,
                   y = Summer_rel_Wint,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Season by family

#How many families are there?
length(unique(sb.mod$Family)) #233
#If I divide by 6, I get 39 page, but then I need to somehow account for sharing orders.  I'm going to guess at 60

pdf("images/manuscript/fig7new/S2_MG_explore_bac_family_rankings_by_season.pdf")
for (i in 1:60){
  print(ggplot(data = sb.mod,
               aes(x = rank_season,
                   y = Summer_rel_Wint,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels

pdf("images/manuscript/fig7new/S2_MG_explore_bac_phylum_to_genus_rankings_by_season.pdf")
for (i in 1:100){
  print(ggplot(data = sb.mod,
               aes(x = rank_season,
                   y = Summer_rel_Wint,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()



#Venn Diagram to help me visualize my hypotheses-----
library(venn)
library(ggplot2)
library(ggpolypath)

sb.mod.dom <- sb.mod %>%
  filter(rare_dom=="dominant")

#Site and Season extremes
venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season==1],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season==10],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10])

venn(venn.list,  zcolor = "style")

lapply(venn.list, length)

#Site, Season agnostic
venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(5,6)],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(5,6)],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10])

venn(venn.list,  zcolor = "style")

#Season, Site agnostic

venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site%in% c(5,6)],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(1)],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(10)],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site %in% c(5,6)])

venn(venn.list,  zcolor = "style")








#Stacked Bar plot equivalent-----

#site
ggplot(data = sb.mod,
       aes(x = fct_reorder(Phylum, rank_site, .fun='median'), 
           y = rank_site,
           fill = Phylum,
           color = Phylum)) +
  geom_jitter(alpha = 0.5,
             size = 3) +
  geom_boxplot(alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "none") +
  labs(y = "Rank",
       x = "Phylum") +
  geom_hline(aes(yintercept = min(rank_site[percentile_site==10])),
             size = 0.5,
             linetype=2) +
  geom_hline(aes(yintercept = max(rank_site[percentile_site==1])),
             size = 0.5,
             linetype=2) +
  annotate("rect", 
           ymin=min(sb.mod$rank_site[sb.mod$percentile_site==1]), 
           ymax=max(sb.mod$rank_site[sb.mod$percentile_site==1]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  annotate("rect", 
           ymin=min(sb.mod$rank_site[sb.mod$percentile_site==10]), 
           ymax=max(sb.mod$rank_site[sb.mod$percentile_site==10]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  geom_text(aes(x= 3, 
                y= 500, 
                label= "Heavily Oiled", 
                size = 14),
            color = "black",
            hjust = "left") +
  geom_text(aes(x= 3, 
                y= 7900, 
                label= "Lightly Oiled", 
                size = 14),
            color = "black",
            hjust = "left") +
  facet_wrap(~ rare_dom)

ggsave("images/manuscript/S2_MG_boxplot_rank_phylum_site.png", width = 10, height = 7, units = "in")

#season
ggplot(data = sb.mod,
       aes(x = fct_reorder(Phylum, rank_season, .fun='median'), 
           y = rank_season,
           fill = Phylum,
           color = Phylum)) +
  geom_jitter(alpha = 0.5,
              size = 3) +
  geom_boxplot(alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "none") +
  labs(y = "Rank",
       x = "Phylum") +
  geom_hline(aes(yintercept = min(rank_season[percentile_season==10])),
             size = 0.5,
             linetype=2) +
  geom_hline(aes(yintercept = max(rank_season[percentile_season==1])),
             size = 0.5,
             linetype=2) +
  annotate("rect", 
           ymin=min(sb.mod$rank_season[sb.mod$percentile_season==1]), 
           ymax=max(sb.mod$rank_season[sb.mod$percentile_season==1]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  annotate("rect", 
           ymin=min(sb.mod$rank_season[sb.mod$percentile_season==10]), 
           ymax=max(sb.mod$rank_season[sb.mod$percentile_season==10]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  geom_text(aes(x= 3, 
                y= 500, 
                label= "Winter", 
                size = 14),
            color = "black",
            hjust = "left") +
  geom_text(aes(x= 3, 
                y= 7900, 
                label= "Summer", 
                size = 14),
            color = "black",
            hjust = "left") +
  facet_wrap(~ rare_dom,
             scales = "free")


ggsave("images/manuscript/S2_MG_boxplot_rank_phylum_season.png", width = 10, height = 7, units = "in")

#Overview of Taxa by % OTUs----
sb.mod %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#percent by reads
fid <- prune_taxa(sb.mod$featureid, bac.2season_with_outliers) 
sb.mod$tax.sums <- fid %>%
  taxa_sums() 

sb.mod %>%
  group_by(Phylum) %>%
  summarise (n = sum(tax.sums)) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  group_by(Phylum, rare_dom) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(Phylum %in% c("Proteobacteria",
                       "Planctomycetes",
                       "Chloroflexi",
                       "Acidobacteria",
                       "Bacteroidetes"))

#Taxa Overlap----

#Heavily Oiled-----
#Taxa that are positively associated with Heavily Oiled
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

#Taxa that are positively associated with Heavily Oiled and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Heavily Oiled and Winter
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))
#Taxa that are positively associated with Heavily Oiled and no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Lightly Oiled------
#Taxa that are positively associated with Lightly Oiled
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

#Taxa that are positively associated with Lightly Oiled and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Lightly Oiled and Winter
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Lightly Oiled and no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Winter and no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Summer and no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#12 HC Orders-------

sb.mod %>%
           filter(Order %in% c("Alteromonadales", 
                                 "Chromatiales", 
                                 "Desulfobacterales",
                                 "Desulfovibrionales", 
                                 "Desulfuromonadales", 
                                 "Enterobacteriales", 
                                 "Methylococcales", 
                                 "Oceanospirillales", 
                                 "Pseudomonadales", 
                                 "Rhizobiales",
                                 "Rhodobacterales",
                                 "Thiotrichales")) %>%
           ggplot(aes(x = rank_site,
                      y = Fourchon_rel_BJ,
                      fill = rare_dom,
                      color = rare_dom)) +
           geom_col(width = 1) +
           theme_bw() +
           labs(y = "Differential with Respect to Site",
                x = "Rank") +
           facet_wrap(~ Order) +
           scale_color_manual(values = c("red", "darkgray")) +
           scale_fill_manual(values = c("red", "darkgray")) +
           theme(legend.position = "bottom")

ggsave("images/manuscript/S2_MG_Fig8_v2.png", width = 10, height = 8, units = "in")   

sb.mod %>%
           filter(Order %in% c("Alteromonadales", 
                                 "Chromatiales", 
                                 "Desulfobacterales",
                                 "Desulfovibrionales", 
                                 "Desulfuromonadales", 
                                 "Enterobacteriales", 
                                 "Methylococcales", 
                                 "Oceanospirillales", 
                                 "Pseudomonadales", 
                                 "Rhizobiales",
                                 "Rhodobacterales",
                                 "Thiotrichales")) %>%
           ggplot(aes(x = rank_season,
                      y = Summer_rel_Wint,
                      fill = rare_dom,
                      color = rare_dom)) +
           geom_col(width = 1) +
           theme_bw() +
           labs(y = "Differential with Respect to Season",
                x = "Rank") +
           facet_wrap(~ Order) +
           scale_color_manual(values = c("red", "darkgray")) +
           scale_fill_manual(values = c("red", "darkgray")) +
           theme(legend.position = "bottom")

ggsave("images/manuscript/S2_MG_Fig8_v2_season.png", width = 10, height = 8, units = "in")   

#Plots of Taxa that are positively associated with site regardless of season------
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(1) &
           percentile_season %in% c(5,6)) %>%
  ggplot(aes(x = fct_reorder(Order, rank_site, .fun='median'), 
           y = rank_site,
           fill = Order,
           color = Order)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "none") +
  labs(y = "Rank",
       x = "Order") +
  facet_wrap(~percentile_site,
             scales = "free")

ggsave("images/manuscript/S2_MG_phylum_site_no_season.png", width = 10, height = 7, units = "in")


#Taxa that are positively associated with season regardless of site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1,10)) %>%
           ggplot(aes(x = fct_reorder(Order, rank_season, .fun='median'), 
                      y = rank_season,
                      fill = Order,
                      color = Order)) +
           geom_point(size = 3) +
           theme_bw() +
           theme(axis.text.x = element_text(hjust = 1, angle = 45),
                 legend.position = "none") +
           labs(y = "Rank",
                x = "Order") +
           facet_wrap(~percentile_season,
                      scales = "free")
         
ggsave("images/manuscript/S2_MG_phylum_season_no_site.png", width = 10, height = 7, units = "in")

#Taxa that are positively associated with no season or site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(5,6)) %>%
  ggplot(aes(x = fct_reorder(Order, rank_season, .fun='median'), 
             y = rank_season,
             fill = Order,
             color = Order)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "none") +
  labs(y = "Rank by Season",
       x = "Order")

ggsave("images/manuscript/S2_MG_phylum_no_season_no_site.png", width = 10, height = 7, units = "in")




#Fungi-----

#load Songbird output-----

sb <- read.delim("songbird/logdir/new_final_3/fungi/S2_MG_fungi_site_season/differentials.tsv")

#rename columns
names(sb) <- c("featureid", "Intercept", "Fourchon_rel_BJ", "Summer_rel_Wint")

#Using not filtered because the biom used in songbird wasn't filtered this way 
fung.2season.sub <- subset_taxa(fung.2season_with_outliers, rownames(tax_table(fung.2season_with_outliers)) %in% sb$featureid)
 
sb.filtered <- sb[sb$featureid %in% rownames(tax_table(fung.2season_with_outliers_filtered)), ]

#bind taxonomy with songbird results
sb <- cbind(sb, tax_table(fung.2season.sub))

#Factor for whether something is dominant or not

sb$rare_dom <- factor(ifelse(sb$featureid %in% sb.filtered$featureid,
                      "dominant",
                      "rare"))

#Replace NA with "unidentified"
sb[,6:12] <- lapply(sb[,6:12], replace_na, replace = "unidentified")

#I'm basically imitating the QURRO plots: https://www.biorxiv.org/content/10.1101/2019.12.17.880047v1.full

#rename in case I need to come back
sb.mod <- sb

#Add rank columns for x axis-----
sb.mod$rank_site <- rank(sb.mod$Fourchon_rel_BJ)
sb.mod$rank_season <- rank(sb.mod$Summer_rel_Wint)

#Add percentile columns 
sb.mod$percentile_site <- cut(sb.mod$rank_site , 
                         breaks = quantile(sb.mod$rank_site, 
                                           seq(from = 0, to = 1, by = 0.1)), 
                         labels=1:10, 
                         include.lowest=TRUE)

sb.mod$percentile_season <- cut(sb.mod$rank_season , 
                         breaks = quantile(sb.mod$rank_season, 
                                           seq(from = 0, to = 1, by = 0.1)), 
                         labels=1:10, 
                         include.lowest=TRUE)


#plot taxonomic ranks as pdf for browsing-----

#Site-----

#How many orders are there?
length(unique(sb.mod$Order)) #20
#If I divide by 6, I get 4 pages

pdf("images/manuscript/fig7new/S2_MG_explore_fung_order_rankings_by_site.pdf")
for (i in 1:4){
  print(ggplot(data = sb.mod,
               aes(x = rank_site,
                   y = Fourchon_rel_BJ,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#site by family

#How many families are there?
length(unique(sb.mod$Family)) #31
#If I divide by 6, I get 6 pages

pdf("images/manuscript/fig7new/S2_MG_explore_fung_family_rankings_by_site.pdf")
for (i in 1:6){
  print(ggplot(data = sb.mod,
               aes(x = rank_site,
                   y = Fourchon_rel_BJ,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels

pdf("images/manuscript/fig7new/S2_MG_explore_fung_phylum_to_genus_rankings_by_site.pdf")
for (i in 1:10){
  print(ggplot(data = sb.mod,
               aes(x = rank_site,
                   y = Fourchon_rel_BJ,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Season-----

#How many orders are there?
length(unique(sb.mod$Order)) #20
#If I divide by 6, I get 4

pdf("images/manuscript/fig7new/S2_MG_explore_fung_order_rankings_by_season.pdf")
for (i in 1:4){
  print(ggplot(data = sb.mod,
               aes(x = rank_season,
                   y = Summer_rel_Wint,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Season by family

#How many families are there?
length(unique(sb.mod$Family)) #31
#If I divide by 6, I get 6

pdf("images/manuscript/fig7new/S2_MG_explore_fung_family_rankings_by_season.pdf")
for (i in 1:6){
  print(ggplot(data = sb.mod,
               aes(x = rank_season,
                   y = Summer_rel_Wint,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels

pdf("images/manuscript/fig7new/S2_MG_explore_fung_phylum_to_genus_rankings_by_season.pdf")
for (i in 1:10){
  print(ggplot(data = sb.mod,
               aes(x = rank_season,
                   y = Summer_rel_Wint,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()



#Venn Diagram to help me visualize my hypotheses-----
library(venn)
library(ggplot2)
library(ggpolypath)

sb.mod.dom <- sb.mod %>%
  filter(rare_dom=="dominant")

#Site and Season extremes
venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season==1],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season==10],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10])

venn(venn.list,  zcolor = "style")

lapply(venn.list, length)

#Site, Season agnostic
venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(5,6)],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(5,6)],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10])

venn(venn.list,  zcolor = "style")

#Season, Site agnostic

venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site%in% c(5,6)],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(1)],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(10)],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site %in% c(5,6)])

venn(venn.list,  zcolor = "style")








#Stacked Bar plot equivalent-----

#site
ggplot(data = sb.mod,
       aes(x = fct_reorder(Phylum, rank_site, .fun='median'), 
           y = rank_site,
           fill = Phylum,
           color = Phylum)) +
  geom_jitter(alpha = 0.5,
             size = 3) +
  geom_boxplot(alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "none") +
  labs(y = "Rank",
       x = "Phylum") +
  geom_hline(aes(yintercept = min(rank_site[percentile_site==10])),
             size = 0.5,
             linetype=2) +
  geom_hline(aes(yintercept = max(rank_site[percentile_site==1])),
             size = 0.5,
             linetype=2) +
  annotate("rect", 
           ymin=min(sb.mod$rank_site[sb.mod$percentile_site==1]), 
           ymax=max(sb.mod$rank_site[sb.mod$percentile_site==1]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  annotate("rect", 
           ymin=min(sb.mod$rank_site[sb.mod$percentile_site==10]), 
           ymax=max(sb.mod$rank_site[sb.mod$percentile_site==10]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  geom_text(aes(x= 3, 
                y= 500, 
                label= "Heavily Oiled", 
                size = 14),
            color = "black",
            hjust = "left") +
  geom_text(aes(x= 3, 
                y= 7900, 
                label= "Lightly Oiled", 
                size = 14),
            color = "black",
            hjust = "left") +
  facet_wrap(~ rare_dom)

ggsave("images/manuscript/S2_MG_boxplot_rank_phylum_site.png", width = 10, height = 7, units = "in")

#season
ggplot(data = sb.mod,
       aes(x = fct_reorder(Phylum, rank_season, .fun='median'), 
           y = rank_season,
           fill = Phylum,
           color = Phylum)) +
  geom_jitter(alpha = 0.5,
              size = 3) +
  geom_boxplot(alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "none") +
  labs(y = "Rank",
       x = "Phylum") +
  geom_hline(aes(yintercept = min(rank_season[percentile_season==10])),
             size = 0.5,
             linetype=2) +
  geom_hline(aes(yintercept = max(rank_season[percentile_season==1])),
             size = 0.5,
             linetype=2) +
  annotate("rect", 
           ymin=min(sb.mod$rank_season[sb.mod$percentile_season==1]), 
           ymax=max(sb.mod$rank_season[sb.mod$percentile_season==1]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  annotate("rect", 
           ymin=min(sb.mod$rank_season[sb.mod$percentile_season==10]), 
           ymax=max(sb.mod$rank_season[sb.mod$percentile_season==10]), 
           xmin= -Inf, 
           xmax=Inf, 
           alpha=0.5, 
           fill="gray") +
  geom_text(aes(x= 3, 
                y= 500, 
                label= "Winter", 
                size = 14),
            color = "black",
            hjust = "left") +
  geom_text(aes(x= 3, 
                y= 7900, 
                label= "Summer", 
                size = 14),
            color = "black",
            hjust = "left") +
  facet_wrap(~ rare_dom,
             scales = "free")


ggsave("images/manuscript/S2_MG_boxplot_rank_phylum_season.png", width = 10, height = 7, units = "in")

#Overview of Taxa by % OTUs----
sb.mod %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#percent by reads
fid <- prune_taxa(sb.mod$featureid, fung.2season_with_outliers) 
sb.mod$tax.sums <- fid %>%
  taxa_sums() 

sb.mod %>%
  group_by(Phylum) %>%
  summarise (n = sum(tax.sums)) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa Overlap----

#Heavily Oiled-----
#Taxa that are positively associated with Heavily Oiled
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1) %>%
  group_by(Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

#Taxa that are positively associated with Heavily Oiled and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==10) %>%
  group_by(Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Heavily Oiled and Winter
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season==1) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Heavily Oiled and no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Lightly Oiled------
#Taxa that are positively associated with Lightly Oiled
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10) %>%
  group_by(Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  print.data.frame()

#Taxa that are positively associated with Lightly Oiled and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==10) %>%
  group_by(Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Lightly Oiled and Winter
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season==1) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Lightly Oiled and no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum,Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Winter and no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa that are positively associated with Summer and no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Phylum, Class) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Family) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10)) %>%
  group_by(Class,Family, Genus) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Figure 7 Bacteria-----

#HC Degrading positive
marinobacter <- ggplot(data = sb.mod %>%
                         filter(Genus=="Marinobacter"),
       aes(x = rank_PAHs,
           y = scaled_PAHs,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "darkgray")) +
  scale_fill_manual(values = c("red", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_PAHs[sb.mod$percentile_PAHs==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_PAHs[sb.mod$percentile_PAHs==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_PAHs)), expand = c(0, 0)) +
  labs(y = "Differential by Oil Abundance",
       x = "Rank",
       subtitle = "Genus: Marinobacter",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  ylim(c(-2,2))

#HC Degrading positive with noise
Alteromonadales <- ggplot(data = sb.mod %>%
                         filter(Order=="Alteromonadales"),
                       aes(x = rank_PAHs,
                           y = scaled_PAHs,
                           fill = rare_dom,
                           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_PAHs[sb.mod$percentile_PAHs==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_PAHs[sb.mod$percentile_PAHs==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_PAHs)), expand = c(0, 0)) +
  labs(y = "Differential by Oil Abundance",
       x = "Rank",
       subtitle = "Order: Alteromonadales") +
  ylim(c(-2,2))

#HC Degrading negative
oceanimonas <- ggplot(data = sb.mod %>%
                         filter(Genus=="Oceanimonas"),
                       aes(x = rank_PAHs,
                           y = scaled_PAHs,
                           fill = rare_dom,
                           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_PAHs[sb.mod$percentile_PAHs==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_PAHs[sb.mod$percentile_PAHs==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_PAHs)), expand = c(0, 0)) +
  labs(y = "Differential by Oil Abundance",
       x = "Rank",
       subtitle = "Genus: Oceanimonas") +
  ylim(c(-2,2))

#Site/Season specific

site.plot <- ggplot(data = sb.mod,
                      aes(x = rank_site,
                          y = Fourchon_rel_BJ,
                          fill = rare_dom,
                          color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "darkgray")) +
  geom_vline(aes(xintercept = min(rank_site[percentile_site==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(rank_site[percentile_site==1])),
             size = 0.5,
             linetype=2) +
  labs(y = "Differential by Site",
       x = "Rank",
       subtitle = "Site") +
  annotate("rect", 
           xmin=min(sb.mod$rank_site[sb.mod$percentile_site==5]), 
           xmax=max(sb.mod$rank_site[sb.mod$percentile_site==6]), 
           ymin= -Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="black") +
  annotate("rect", 
           xmin = 0, 
           xmax=max(sb.mod$rank_site[sb.mod$percentile_site==1]), 
           ymin= -Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="gray") +
  annotate("rect", 
           xmin = min(sb.mod$rank_site[sb.mod$percentile_site==10]), 
           xmax=max(sb.mod$rank_site[sb.mod$percentile_site==10]), 
           ymin= -Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="gray")


season.plot <- ggplot(data = sb.mod,
                    aes(x = rank_season,
                        y = Summer_rel_Wint,
                        fill = rare_dom,
                        color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "darkgray")) +
  geom_vline(aes(xintercept = min(rank_season[percentile_season==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(rank_season[percentile_season==1])),
             size = 0.5,
             linetype=2) +
  labs(y = "Differential by Season",
       x = "Rank",
       subtitle = "Season") +
  annotate("rect", 
           xmin=min(sb.mod$rank_season[sb.mod$percentile_season==5]), 
           xmax=max(sb.mod$rank_season[sb.mod$percentile_season==6]), 
           ymin= -Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="gray") +
  annotate("rect", 
           xmin = 0, 
           xmax=max(sb.mod$rank_season[sb.mod$percentile_season==1]), 
           ymin= -Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="black") +
  annotate("rect", 
           xmin = min(sb.mod$rank_season[sb.mod$percentile_season==10]), 
           xmax=max(sb.mod$rank_season[sb.mod$percentile_season==10]), 
           ymin= -Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="black")


#plot-----

library(cowplot)

Top <- plot_grid(marinobacter, Alteromonadales, nrow = 1, labels = c("A", "B"))
D <- plot_grid(site.plot, season.plot, ncol = 1)
Bottom <- plot_grid(oceanimonas, D, labels = c("C", "D"))

legend <- get_legend(marinobacter + theme(legend.position = "bottom"))

quasi.shebang <- plot_grid(Top, Bottom, nrow = 2)
whole.shebang <- plot_grid(quasi.shebang, legend, ncol = 1, rel_heights = c(1,0.1))

whole.shebang

ggsave("images/manuscript/S2_MG_Fig7.png", height = 8, width = 10, units = "in")

#How do the ranks compare?-----

sb.mod <- sb

#Add rank column for x axis
sb.mod$rank_season <- rank(sb.mod$Summer_rel_Wint)
sb.mod$rank_site <- rank(sb.mod$Fourchon_rel_BJ)
sb.mod$rank_oil <- rank(sb.mod$scaled_PAHs)

plot(sb.mod$rank_oil, sb.mod$rank_season)
plot(sb.mod$rank_oil, sb.mod$rank_site)
plot(sb.mod$rank_site, sb.mod$rank_season)
      
#plot of oil-associated log ratios and oil abundance------

oil.ass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           Genus!="unidentified" &
           percentile_PAHs==10) %>%
  rownames()

oil.noass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs %in% c(4:7) &
           percentile_site %in% c(4:7) &
           percentile_season %in% c(4:7)) %>%
  rownames()

oil.ass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.ass.list,]  

oil.noass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.noass.list,]  

oil.ass.sums <- colSums(oil.ass.OT)
oil.noass.sums <- colSums(oil.noass.OT)

log.oil.ass <- log(oil.ass.sums/oil.noass.sums)

plot(log.oil.ass, log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))

mod <- lm(log.oil.ass ~ sample_data(bac.2season_with_outliers)$Total.relevant.PAHs)

summary(mod)
hist(resid(mod))
shapiro.test(resid(mod))
qqnorm(resid(mod))
abline(0,1)

#Weird that deciles 4,5 are normal, but 5,6 are not.  To investigate

#Are any of the oil taxa real weirdos?

#Where do deciles fall around inflection point?
ggplot(data = sb.mod,
                    aes(x = rank_PAHs,
                        y = scaled_PAHs,
                        fill = rare_dom,
                        color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "darkgray")) +
  geom_vline(aes(xintercept = min(rank_PAHs[percentile_PAHs==4])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(rank_PAHs[percentile_PAHs==5])),
             size = 0.5,
             linetype=2) +
  labs(y = "Differential by Oil",
       x = "Rank",
       subtitle = "Oil")

#negative oil taxa is not very convincing

#features: 239671, 825763, 510987, 818915  

oil.ass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1) %>%
  rownames()

oil.noass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10) %>%
  rownames()

#oil.ass.list <- c("239671", "825763", "510987", "818915") #top 4 negative associates

oil.ass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.ass.list,]  

oil.noass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.noass.list,]  

oil.ass.sums <- colSums(oil.ass.OT)
oil.noass.sums <- colSums(oil.noass.OT)

log.oil.ass <- log(oil.ass.sums/oil.noass.sums)

plot(log.oil.ass, log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))

mod <- lm(log.oil.ass ~ sample_data(bac.2season_with_outliers)$Total.relevant.PAHs)

summary(mod)
hist(resid(mod))
shapiro.test(resid(mod))
qqnorm(resid(mod))
abline(0,1)

#What does middle ground look like?

oil.ass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs %in% c(5)) %>%
  rownames()

oil.noass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs %in% c(6)) %>%
  rownames()

oil.ass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.ass.list,]  

oil.noass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.noass.list,]  

oil.ass.sums <- colSums(oil.ass.OT)
oil.noass.sums <- colSums(oil.noass.OT)

log.oil.ass <- log(oil.ass.sums/oil.noass.sums)

plot(log.oil.ass, log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))

mod <- lm(log.oil.ass ~ log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))

summary(mod)
hist(resid(mod))
qqnorm(resid(mod))
abline(0,1)

#none of the above is super convincing.


