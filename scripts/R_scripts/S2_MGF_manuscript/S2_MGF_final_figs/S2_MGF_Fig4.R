#Ranked Differential plots

#by Steve Formel
#Last updated 11 Nov, 2020
#Make ranked abundance plot for songbird analysis

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(ggforce)
library(cowplot)
library(ggrepel)
library(grid)
library(venn)
library(ggplot2)
library(ggpolypath)

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


# Season----

#How many orders are there?
length(unique(sb.mod$Class)) #11
#If I divide by 4, I get 3 pages

text_high <- textGrob("Summer ", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("Winter", gp=gpar(fontsize=12, fontface="bold"), just = "center")
    
pdf("images/manuscript/S2_MGF_final/S2_MG_explore_fung_class_rankings_by_season.pdf")
for (i in 1:3){
  print(
    
    
    ggplot(data = sb.mod,
       aes(x = rank_season,
           y = Summer_rel_Wint,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgrey")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_season[sb.mod$percentile_season==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_season[sb.mod$percentile_season==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_season)), expand = c(0, 0)) +
  labs(y = "Differential by Season",
       x = "Rank",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(percentile_season %in% c(1,10)),
                  aes(x = rank_season,
                      y = Summer_rel_Wint,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5) +
         facet_wrap_paginate(~ Phylum*Class, ncol = 2, nrow = 2, page = i))
  }  #end loop

dev.off()


#Site ----

text_high <- textGrob("Lightly Oiled", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("Heavily Oiled", gp=gpar(fontsize=12, fontface="bold"), just = "center")
    
pdf("images/manuscript/S2_MGF_final/S2_MG_explore_fung_class_rankings_by_site.pdf")
for (i in 1:3){
  print(
    ggplot(data = sb.mod,
       aes(x = rank_site,
           y = Fourchon_rel_BJ,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgrey")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_site[sb.mod$percentile_site==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_site[sb.mod$percentile_site==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_site)), expand = c(0, 0)) +
  labs(y = "Differential by Site",
       x = "Rank",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(percentile_site %in% c(1,10)),
                  aes(x = rank_site,
                      y = Fourchon_rel_BJ,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5) +
         facet_wrap_paginate(~ Phylum*Class, ncol = 2, nrow = 2, page = i))
  }  #end loop

dev.off()


#From the above I think I should plot:

#1. Dothideomycetes by site and season
#2. Sordariomycetes by site
#3. Agaricomycetes by season
#4. Hydrocarbon degraders by site

#Fig 4------

#Dothideomycetes taxa by season

text_high <- textGrob("Summer \u2192", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("\u2190 Winter", gp=gpar(fontsize=12, fontface="bold"), just = "center")

doth.season <- ggplot(data = sb.mod %>%
         filter(Class=="Dothideomycetes"),
       aes(x = rank_season,
           y = Summer_rel_Wint,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_season[sb.mod$percentile_season==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_season[sb.mod$percentile_season==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_season)), expand = c(0, 0)) +
  labs(y = "Differential by Season",
       x = "Rank",
       title = "Dothideomycetes by Season",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Class=="Dothideomycetes" & percentile_season %in% c(1,10)),
                  aes(x = rank_season,
                      y = Summer_rel_Wint,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)

#Dothideoomycetes taxa by site

text_high <- textGrob("Lightly Oiled \u2192", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("\u2190 Heavily Oiled", gp=gpar(fontsize=12, fontface="bold"), just = "center")

doth.site <- ggplot(data = sb.mod %>%
         filter(Class=="Dothideomycetes"),
       aes(x = rank_site,
           y = Fourchon_rel_BJ,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_site[sb.mod$percentile_site==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_site[sb.mod$percentile_site==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_site)), expand = c(0, 0)) +
  labs(y = "Differential by Site",
       x = "Rank",
       title = "Dothideomycetes by Site",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Class=="Dothideomycetes" & percentile_site %in% c(1,10)),
                  aes(x = rank_site,
                      y = Fourchon_rel_BJ,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)

doth.all <- plot_grid(doth.site, doth.season)

#Agaricomycetes taxa by season

text_high <- textGrob("Summer \u2192", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("\u2190 Winter", gp=gpar(fontsize=12, fontface="bold"), just = "center")

agarico <- ggplot(data = sb.mod %>%
         filter(Class=="Agaricomycetes"),
       aes(x = rank_season,
           y = Summer_rel_Wint,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_season[sb.mod$percentile_season==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_season[sb.mod$percentile_season==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_season)), expand = c(0, 0)) +
  labs(y = "Differential by Season",
       x = "Rank",
       title = "Agaricomycetes",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Class=="Agaricomycetes" & percentile_season %in% c(1,10)),
                  aes(x = rank_season,
                      y = Summer_rel_Wint,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)

#Sordariomycetes taxa by site

text_high <- textGrob("Lightly Oiled \u2192", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("\u2190 Heavily Oiled", gp=gpar(fontsize=12, fontface="bold"), just = "center")

sordario <- ggplot(data = sb.mod %>%
         filter(Class=="Sordariomycetes"),
       aes(x = rank_site,
           y = Fourchon_rel_BJ,
           fill = rare_dom,
           color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgray")) +
  geom_vline(aes(xintercept = min(sb.mod$rank_site[sb.mod$percentile_site==10])),
             size = 0.5,
             linetype=2) +
  geom_vline(aes(xintercept = max(sb.mod$rank_site[sb.mod$percentile_site==1])),
             size = 0.5,
             linetype=2) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_site)), expand = c(0, 0)) +
  labs(y = "Differential by Site",
       x = "Rank",
       title = "Sordariomycetes",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Class=="Sordariomycetes" & percentile_site %in% c(1,10)),
                  aes(x = rank_site,
                      y = Fourchon_rel_BJ,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)

#plot figure 4----

plot_grid(doth.site + labs(title = NULL,subtitle = "Dothideomycetes"), doth.season + labs(title = NULL,subtitle = "Dothideomycetes"), sordario + labs(title = NULL,subtitle = "Sordarioomycetes"), agarico + labs(title = NULL,subtitle = "Agaricomycetes"), labels = "AUTO")

ggsave("images/manuscript/S2_MGF_final/S2_MGF_Fig4.png", width = 10, height = 8, units = "in")

#Genera found in Site and Seasons-----

#Heavily Oiled
BJ.genus <- sb.mod %>%
  filter(percentile_site %in% c(1)) %>%
  select(Genus)

BJ.genus$Genus[BJ.genus$Genus %in% all.HC.list]

#Lightly Oiled
F.genus <- sb.mod %>%
  filter(percentile_site %in% c(10)) %>%
  select(Genus)

F.genus$Genus[F.genus$Genus %in% all.HC.list]

#Winter but no site
sb.mod %>%
  filter(percentile_season %in% c(1) & percentile_site %in% c(5,6)) %>%
  select(Family)

#Summer but no site
sb.mod %>%
  filter(percentile_season %in% c(10) & percentile_site %in% c(5,6)) %>%
  select(Genus)

sb.mod %>%
  filter(Genus %in% all.HC.list & percentile_site %in% c(1,10))

sb.mod %>%
  filter(Genus %in% all.HC.list & percentile_season %in% c(1,10))


#Venn Diagram to help me visualize my hypotheses-----

sb.mod.dom <- sb.mod %>%
  filter(rare_dom=="dominant")

#Site and Season extremes
venn.list <- list('Heavily Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season==1],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season==10],
                  'Lightly Oiled' = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10])

venn(venn.list, box = FALSE, ggplot = TRUE, ellipse = TRUE, )

ggsave("images/manuscript/S2_MG_fungi_venn.png", width = 5, height = 5, units = "in")

lapply(venn.list, length)

#Make Table of Extreme Taxa-----

x.tax <- plyr:::ldply(venn.list, data.frame)

x.taxon <- sb.mod[match(x.tax$X..i.., sb.mod$featureid),] %>%
  select(Phylum, Class, Order, Family, Genus, Species)

#They match, but a .1 was added to duplicates
x.all$X..i..[which(x.all$X..i..!=rownames(x.taxon))]
rownames(x.taxon)[which(x.all$X..i..!=rownames(x.taxon))]

x.all <- cbind(x.tax, x.taxon)

names(x.all) <- c("Associate", "OTU id", "Phylum", "Class", "Order", "Family", "Genus", "Species")

write.csv(x.all, "images/manuscript/S2_MGF_final/S2_MGF_extreme_ranks.table.csv", row.names = FALSE)

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