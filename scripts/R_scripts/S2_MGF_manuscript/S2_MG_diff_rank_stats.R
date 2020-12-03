#Ranked Differential plots

#by Steve Formel
#Last updated 24 Sep, 2020
#Make ranked abundance plot for songbird analysis

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(ggforce)

#Bacteria-----

#load Songbird output-----

sb <- read.delim("songbird/logdir/new_final_3/bac/S2_MG_bac_site_season_oil/differentials.tsv")

#sb <- read.delim("songbird/logdir/new_final_2/bac/S2_MG_bac_site_season_oil/differentials.tsv")

#Old results that don't align with new ones.  I can't figure out what was different, but don't use.
#sb <- read.delim("songbird/logdir/new_final/S2_MG_bac_site_season_oil/differentials.tsv")

#rename columns
names(sb) <- c("featureid", "Intercept", "Fourchon_rel_BJ", "Summer_rel_Wint", "scaled_PAHs")

#Using not filtered because the biom used in songbird wasn't filtered this way 
bac.2season.sub <- subset_taxa(bac.2season_with_outliers, rownames(tax_table(bac.2season_with_outliers)) %in% sb$featureid)
 
#For fun subset to filtered version
#bac.2season.sub.filtered <- subset_taxa(bac.2season.sub, rownames(tax_table(bac.2season.sub)) %in% rownames(tax_table(bac.2season)))
 
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
sb.mod$rank_PAHs <- rank(sb.mod$scaled_PAHs)
sb.mod$rank_site <- rank(sb.mod$Fourchon_rel_BJ)
sb.mod$rank_season <- rank(sb.mod$Summer_rel_Wint)

#Add percentile columns 
sb.mod$percentile_PAHs <- cut(sb.mod$rank_PAHs , 
                         breaks = quantile(sb.mod$rank_PAHs, 
                                           seq(from = 0, to = 1, by = 0.1)), 
                         labels=1:10, 
                         include.lowest=TRUE)

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
  group_by(rare_dom) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#BJ no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  nrow()

#BJ no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==1 &
           percentile_season %in% c(5,6)) %>%
  group_by(Phylum, Family, Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n))

#F no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  nrow()

#F no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site==10 &
           percentile_season %in% c(5,6)) %>%
  group_by(Family, Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n))

#Winter no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_season==1 &
           percentile_site %in% c(5,6)) %>%
  group_by(Phylum, Family, Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n))

#Summer no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_season==10 &
           percentile_site %in% c(5,6)) %>%
  group_by(Phylum, Family, Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n))


#No site or season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(5,6))

#Oil-----

#Taxa that are positively associated with oil and BJ regardless of season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10) %>%
  nrow()

#Oil
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10) %>%
  group_by(Phylum) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

#Oil
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10) %>%
  group_by(Class) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10) %>%
  group_by(Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

#overlap site & season-----

#Taxa that are positively associated with oil and BJ regardless of season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==1 &
           percentile_season %in% c(5,6))

#Oil and Fourchon no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==10 &
           percentile_season %in% c(5,6))

#Taxa that are positively associated with oil and Winter regardless of site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site %in% c(5,6) &
           percentile_season==1)

#Taxa that are positively associated with oil and Summer regardless of site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site %in% c(5,6) &
           percentile_season==10)

#Taxa Overlap----

#Taxa that are positively associated with oil
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10) %>%
  nrow()

#Taxa that are positively associated with oil and BJ
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==1) %>%
  nrow()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==1) %>%
  group_by(Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

#Taxa that are positively associated with oil and F
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==10) %>%
  nrow()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==10) %>%
  group_by(Class, Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==10) %>%
  nrow()
#Taxa that are negatively associated with oil
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1) %>%
  nrow()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1) %>%
  group_by(Genus) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

2#Taxa that are negatively associated with oil and positive with BJ
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           percentile_site==1) %>%
  nrow()

#Taxa that are negatively associated with oil and positive with F
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           percentile_site==10) %>%
  nrow()

#What is the overlap of the oil and the two sites together?  What is common to the two sites?

#Can't answer this because oil may be the reason for the differences between the two sites.  And It's impossible for the sites to share anything, except in the sense of the middle 50th percentile.  That is, the sites had differential oiling, so we can't assume the oiling was a controlled event like winter or summer.

#Taxa that are positively associated with oil and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_season==10) %>%
  nrow()

sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_season==10) %>%
  group_by(Class) %>%
  summarise (n = n()) %>%
  arrange(desc(n)) %>%
  print.data.frame()

#Taxa that are positively associated with oil and Winter
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_season==1)

#Taxa that are negatively associated with oil and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           percentile_season==10)

#Taxa that are negatively associated with oil and Winter
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           percentile_season==1)


#Taxa that are positively associated with site regardless of season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(1,10) &
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
           ggplot(aes(x = rank_PAHs,
                      y = scaled_PAHs,
                      fill = rare_dom,
                      color = rare_dom)) +
           geom_col(width = 1) +
           theme_bw() +
           labs(y = "Differential with Respect to Oil",
                x = "Rank") +
           facet_wrap(~ Order) +
           scale_color_manual(values = c("red", "darkgray")) +
           scale_fill_manual(values = c("red", "darkgray")) +
           theme(legend.position = "bottom")

ggsave("images/manuscript/S2_MG_Fig8_v1.png", width = 10, height = 8, units = "in")         



#14 oil genera OTUs-----
         
         oil.ass.list <- sb.mod %>%
           filter(rare_dom=="dominant" &
                    percentile_PAHs==10 &
                    Genus!="unidentified") %>%
  rownames()
         
         oil.noass.list <- sb.mod %>%
           filter(percentile_PAHs %in% c(5,6)) %>%
           rownames()
         
         oil.ass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.ass.list,]  
         
         oil.noass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.noass.list,]  
         
         oil.ass.sums <- colSums(oil.ass.OT)
         oil.noass.sums <- colSums(oil.noass.OT)
         
         log.oil.ass <- log(oil.ass.sums/oil.noass.sums)
         log.oil.ass[is.infinite(log.oil.ass)] <- NA
         
         plot(log.oil.ass, log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))
         
         mod <- lm(log.oil.ass ~ log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))
         
         summary(mod)
         hist(resid(mod))
         shapiro.test(resid(mod))
         qqnorm(resid(mod))
         abline(0,1)
         
         
#not quite normal

library(brms)

df <- data.frame(log.oil.ass, "oil" = sample_data(bac.2season_with_outliers)$Total.relevant.PAHs)
                 
b0 <- brm(log.oil.ass ~ log(oil),
          data = df,
          family = "skew_normal",
          chains = 4, cores = 4,
          seed = 1)

M <- b0
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
posterior_summary(M)

#R2
bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975))


#Negatively associated oil genera OTUs-----
oil.ass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           Genus!="unidentified") %>%
  rownames()

oil.noass.list <- sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs %in% c(5,6)) %>%
  rownames()

oil.ass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.ass.list,]  

oil.noass.OT <- otu_table(bac.2season_with_outliers)[rownames(otu_table(bac.2season_with_outliers)) %in% oil.noass.list,]  

oil.ass.sums <- colSums(oil.ass.OT)
oil.noass.sums <- colSums(oil.noass.OT)

log.oil.ass <- log(oil.ass.sums/oil.noass.sums)
log.oil.ass[is.infinite(log.oil.ass)] <- NA

plot(log.oil.ass, log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))

mod <- lm(log.oil.ass ~ log(sample_data(bac.2season_with_outliers)$Total.relevant.PAHs))

summary(mod)
hist(resid(mod))
shapiro.test(resid(mod))
qqnorm(resid(mod))
abline(0,1)

#Fungi-----
         
         
#load Songbird output-----
         
sb <- read.delim("songbird/logdir/new_final/S2_MG_fung_site_season_oil/differentials.tsv")
         
         