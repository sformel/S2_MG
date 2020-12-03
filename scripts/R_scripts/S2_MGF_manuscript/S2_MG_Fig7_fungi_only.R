#Ranked Differential plots

#by Steve Formel
#Last updated 26 Oct, 2020
#Make ranked abundance plot for songbird analysis

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(ggforce)
library(cowplot)
library(Cairo)
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

library(ggrepel)
library(grid)

#How many orders are there?
length(unique(sb.mod$Class)) #11
#If I divide by 4, I get 3 pages

text_high <- textGrob("Summer ", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("Winter", gp=gpar(fontsize=12, fontface="bold"), just = "center")
    
pdf("images/manuscript/S2_MG_explore_fung_class_rankings_by_season.pdf")
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
  scale_color_manual(values = c("red", "darkgray")) +
  scale_fill_manual(values = c("red", "darkgrey")) +
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
    
pdf("images/manuscript/S2_MG_explore_fung_class_rankings_by_site.pdf")
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
  scale_color_manual(values = c("red", "darkgray")) +
  scale_fill_manual(values = c("red", "darkgrey")) +
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

#Fig 5------

#Dothideomycetes taxa by season
library(ggrepel)
library(grid)

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
library(ggrepel)
library(grid)

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
library(ggrepel)
library(grid)

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
library(ggrepel)
library(grid)

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

#Hydrocarbon Degraders

#From Prenafeta-Oldu et al.
     
Prenafeta.HC <- c("Acremonium", "Alternaria", "Aspergillus", "Aureobasidium", "Botrytis",
"Candida", "Cephalosporium", "Chaetomium", "Chrysosporium",
"Cladosporium", "Curvularia", "Drechslera", "Epicoccum", "Geomyces", "Geotrichum",
"Gliomastix", "Fusarium", "Hansenula", "Helminthosporium", "Humicola", "Mucor", "Paecilomyces", "Penicillium", "Pestalotiopsis","Phialophora", "Phoma", "Phomopsis", "Pseudallescheria", "Scedosporium","Rhinocladiella", "Rhizopus", "Rhodotorula", "Saccharomyces", "Sordaria", "Stemphylium","Thielavia", "Trichoderma", "Trichosporon", "Trichothecium", "Tritirachium",
"Ulocladium")

Prenafeta.PAH.list <- c("Phanerochaete", "Bjerkandera","Trametes","Agrocybe", "Clitocybe", "Collybia", "Marasmius", "Mycena","Stropharia","Aspergillus", "Penicillium","Fusarium", "Beauveria","Verticillium","Paecilomyces","Penicillium")

#Muncnerova et al
Muncnerova <- c("Phanerochaete", "Aspergillus", "Cunninghamella", "Bjerkandera", "Trametes", "Saccharomyces")

#Furuno et al
Furuno <- c("Pythium")

#Blasi et al
Blasi <- c("Exophiala", "Cladophialophora")

daSilva <- ("Cyclothyrium")
#Verkley et al.

Verkley <- c("Paraconiothyrium") 
             
#Hashem et al

Hashem <- c("Yamadazyma", "Rhodotorula", "Pichia", "Candida", "Meyerozyma")

#Cerniglia 1997
Cerniglia <- c("Cunninghamella", "Crinipellis", "Phanerochaete", "Bjerkandera","Trametes", "Pleurotus", "Aspergillus", "Penicillium", "Lentinus", "Naematoloma")

#Atlas 1981
Atlas <- c("Candida", "Rhodotorula", "Sporobolomyces", "Penicillium", "Cunninghamella", "Penicillium", "Verticillium", "Beauveria", "Mortieriella", "Phoma", "Scolecobasidium", "Tolypocladium", "Hansenula", "Aureobasidium", "Cladosporium", "Penicillium", "Aspergillus", "Rhodosporidium", "Saccharomyces", "Trichosporon", "Graphium", "Fusarium", "Paecilomyces", "Acremonium", "Gliocladium", "Trichoderma", "Sphaeropsidales")

#Kirk and Gordon 1988 - alkanes
Kirk <- c("Corollospora", "Dendryphiella", "Lulworthia", "Varicosporina")

#Prince et al (Handbook of Hycrocarbon and Lipid...2010)

Prince <- read.csv("Prince_et_al_fungi_HC_degrader_list.csv")
Prince.genera <- Prince %>%
  select(Genus) 

Prince.genera <- Prince.genera$Genus
Prince.PAH <- Prince %>%
  filter(Typical_Substrate %in% c("Chrysene", "Phenanthrene", "Naphthalene", "Dibenzothiophene")) %>%
  select(Genus)

Prince.PAH <- Prince.PAH$Genus

all.HC.list <- unique(c(Prenafeta.HC, Prenafeta.PAH.list, Muncnerova, Blasi, daSilva, Verkley, Hashem, Cerniglia, Atlas, Kirk, Prince.genera))

#How many of these show up in our data?
in.our.data <- sb.mod$Genus[sb.mod$Genus %in% all.HC.list] %>%
  unique() %>%
  sort()

paste(noquote(sort(unique(in.our.data))),collapse=', ')

#full classification
sb.mod[sb.mod$Genus %in% all.HC.list,] %>%
  select(Phylum, Order, Family, Genus, Species)

#Fungi that have been shown to interact with the PAHs we measured
S2_PAH.list <- unique(c("Aspergillus", "Phanerochaete", "Cunninghamella", "Pleurotus",  "Naematoloma", "Syncephalastrum", "Trametes", "Penicillium", "Saccharomyces", "Candida", "Pichia", "Paraconiothyrium", "Cyclothyrium", Prince.PAH))

PAH.in.our.data <- sb.mod$Genus[sb.mod$Genus %in% S2_PAH.list] %>%
  unique() %>%
  sort()

paste(noquote(sort(unique(PAH.in.our.data))),collapse=', ')
#print list for manuscript

paste(noquote(sort(unique(all.HC.list))),collapse=', ')
paste(noquote(sort(unique(S2_PAH.list))),collapse=', ')

#Orders (Prenafeta)
HC.orders <- c("Dothideales","Chaetothyriales","Hypocreales","Ophiostomales")

text_high <- textGrob("Lightly Oiled \u2192", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("\u2190 Heavily Oiled", gp=gpar(fontsize=12, fontface="bold"), just = "center")

#Genus

HC.site <- ggplot(data = sb.mod %>%
         filter(Genus %in% all.HC.list),
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
       title = "Hydrocarbon Degraders",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.2, ymax = -9.2) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.2, ymax = -9.2) +
  coord_cartesian(ylim=c(-8,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Genus %in% all.HC.list & percentile_site %in% c(1,10)),
                  aes(x = rank_site,
                      y = Fourchon_rel_BJ,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)


#HC by Season

text_high <- textGrob("Summer \u2192", gp=gpar(fontsize=12, fontface="bold"), just = "center")
text_low <- textGrob("\u2190 Winter", gp=gpar(fontsize=12, fontface="bold"), just = "center")

HC.season <- ggplot(data = sb.mod %>%
         filter(Genus %in% all.HC.list),
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
  labs(y = "Differential by season",
       x = "Rank",
       title = "Hydrocarbon Degraders",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -3.4, ymax = -3.4) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -3.4, ymax = -3.4) +
  coord_cartesian(ylim=c(-2.75,3.5), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Genus %in% all.HC.list & percentile_season %in% c(1,10)),
                  aes(x = rank_season,
                      y = Summer_rel_Wint,
                      label = Genus),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)

#HC by Order
ggplot(data = sb.mod %>%
         filter(Order %in% HC.orders),
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
       title = "Hydrocarbon Degraders",
       color = "Relative Abundance of Taxa",
       fill = "Relative Abundance of Taxa") +
  annotation_custom(text_low, xmin = 20, xmax = 20, ymin = -9.5, ymax = -9.5) +
  annotation_custom(text_high, xmin = 90, xmax = 90, ymin = -9.5, ymax = -9.5) +
  coord_cartesian(ylim=c(-8,3), clip="off") +
  geom_label_repel(inherit.aes = FALSE,
                  data = sb.mod %>%
                    filter(Order %in% HC.orders & percentile_site %in% c(1,10)),
                  aes(x = rank_site,
                      y = Fourchon_rel_BJ,
                      label = Order),
                  nudge_x = 0.5,
                  direction = "y",
                  angle = 0,
                  vjust = 0,
                  hjust = 0,
                  segment.size = 0.5)

#Site/Season specific

site.plot <- ggplot(data = sb.mod,
                      aes(x = rank_site,
                          y = Fourchon_rel_BJ,
                          fill = rare_dom,
                          color = rare_dom)) +
  geom_col(width = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "darkgray")) +
  scale_x_continuous(limits = c(0, max(sb.mod$rank_season)), expand = c(0, 0)) +
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

#plot Figures 4 and 5

#1. Dothideomycetes by site and season
#2. Sordariomycetes by site
#3. Agaricomycetes by season
#4. Hydrocarbon degraders by site

Fig5 <- plot_grid(HC.site + labs(title = NULL), HC.season + labs(title = NULL), labels = "AUTO")

ggsave("images/manuscript/S2_MG_fungi_Fig5.png", width = 8, height = 5, units = "in")
  
Fig4 <- plot_grid(doth.site + labs(title = NULL,subtitle = "Dothideomycetes"), doth.season + labs(title = NULL,subtitle = "Dothideomycetes"), sordario + labs(title = NULL,subtitle = "Sordarioomycetes"), agarico + labs(title = NULL,subtitle = "Agaricomycetes"), labels = "AUTO")

ggsave("images/manuscript/S2_MG_fungi_Fig4.png", width = 10, height = 8, units = "in")


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

write.csv(x.all, "images/manuscript/extreme_ranks.table.csv", row.names = FALSE)

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
  scale_color_manual(values = c("black", "darkgray")) +
  scale_fill_manual(values = c("black", "darkgray")) +
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
  scale_color_manual(values = c("black", "darkgray")) +
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
  scale_color_manual(values = c("black", "darkgray")) +
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
  scale_color_manual(values = c("black", "darkgray")) +
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
  scale_color_manual(values = c("black", "darkgray")) +
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
  scale_color_manual(values = c("black", "darkgray")) +
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


