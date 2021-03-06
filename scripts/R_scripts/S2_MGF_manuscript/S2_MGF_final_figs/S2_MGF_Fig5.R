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


#Figure 5-----

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


plot_grid(HC.site + labs(title = NULL), HC.season + labs(title = NULL), labels = "AUTO")

ggsave("images/manuscript/S2_MGF_final/S2_MGF_Fig5.png", width = 8, height = 5, units = "in")
