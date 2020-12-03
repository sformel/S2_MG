#Ranked Differential plots

#by Steve Formel
#Last updated 9 Oct, 2020
#Make ranked abundance plot for songbird analysis

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(ggforce)

#Bacteria-----

#load Songbird output-----

#sb <- read.delim("songbird/logdir/new_final_3/bac/S2_MG_bac_site_season_oil/differentials.tsv")

sb <- read.delim("songbird/logdir/new_final_3/bac/S2_MG_bac_site_season_phen/differentials.tsv")

#rename columns
names(sb) <- c("featureid", "Intercept", "Fourchon_rel_BJ", "Summer_rel_Wint", "scaled_PAHs")

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


#plot taxonomic ranks as pdf for browsing-----

#How many orders are there?
length(unique(sb.mod$Order)) #222
#If I divide by 6, I get 37 pages

pdf("images/manuscript/S2_MG_explore_bac_order_rankings_by_oil.pdf")
for (i in 1:37){
  print(ggplot(data = sb.mod,
       aes(x = rank_PAHs,
           y = scaled_PAHs,
           fill = rare_dom,
           color = rare_dom)) +
         geom_col(width = 1) +
         theme_bw() +
         theme(legend.position = "none",
               axis.text = element_blank()) +
         labs(x = NULL) +
         scale_color_manual(values = c("red", "darkgray")) +
         geom_vline(aes(xintercept = min(rank[percentile==10])),
                    size = 0.5,
                    linetype=2) +
         geom_vline(aes(xintercept = max(rank[percentile==1])),
                    size = 0.5,
                    linetype=2) +
         facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
         
       }

dev.off()

#oil by family----

#How many families are there?
length(unique(sb.mod$Family)) #233
#If I divide by 6, I get 39 page, but then I need to somehow account for sharing orders.  I'm going to guess at 60

pdf("images/manuscript/S2_MG_explore_bac_family_rankings_by_oil.pdf")
for (i in 1:60){
  print(ggplot(data = sb.mod,
               aes(x = rank_PAHs,
                   y = scaled_PAHs,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
    facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels----

pdf("images/manuscript/S2_MG_explore_bac_phylum_to_genus_rankings_by_oil.pdf")
for (i in 1:120){
  print(ggplot(data = sb.mod,
               aes(x = rank_PAHs,
                   y = scaled_PAHs,
                   fill = rare_dom,
                   color = rare_dom)) +
          geom_col(width = 1) +
          theme_bw() +
          theme(legend.position = "none",
                axis.text = element_blank()) +
          labs(x = NULL) +
          scale_color_manual(values = c("red", "darkgray")) +
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()
# 
# #Make table of OTUs that are in the outer 10 percentiles and are "dominant" according to the filtering criteria in "S2_load_packages_and_clean_data.R"
# 
# bac.oil.pos <- sb.mod %>%
#   filter(percentile==10 & filtered=="dominant")
# 
# bac.oil.neg <- sb.mod %>%
#   filter(percentile==1 & filtered=="dominant")

#Site-----

#How many orders are there?
length(unique(sb.mod$Order)) #222
#If I divide by 6, I get 37 pages

pdf("images/manuscript/S2_MG_explore_bac_order_rankings_by_site.pdf")
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
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#site by family

#How many families are there?
length(unique(sb.mod$Family)) #233
#If I divide by 6, I get 39 page, but then I need to somehow account for sharing orders.  I'm going to guess at 60

pdf("images/manuscript/S2_MG_explore_bac_family_rankings_by_site.pdf")
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
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels

pdf("images/manuscript/S2_MG_explore_bac_phylum_to_genus_rankings_by_site.pdf")
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
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()
# 
# #Site Taxa in 10th percentiles and dominant-----
# 
# bac.site.Fourchon <- sb.mod %>%
#   filter(percentile==10 & filtered=="dominant")
# 
# bac.site.BJ <- sb.mod %>%
#   filter(percentile==1 & filtered=="dominant")
# 
# bac.site.middle.pct <- sb.mod %>%
#   filter(percentile %in% c(5,6) & filtered=="dominant")

#Season-----

#How many orders are there?
length(unique(sb.mod$Order)) #222
#If I divide by 6, I get 37 pages

pdf("images/manuscript/S2_MG_explore_bac_order_rankings_by_season.pdf")
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
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Season by family

#How many families are there?
length(unique(sb.mod$Family)) #233
#If I divide by 6, I get 39 page, but then I need to somehow account for sharing orders.  I'm going to guess at 60

pdf("images/manuscript/S2_MG_explore_bac_family_rankings_by_season.pdf")
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
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family, ncol = 2, nrow = 3, page = i))
  
}

dev.off()

#Genus with all other levels

pdf("images/manuscript/S2_MG_explore_bac_phylum_to_genus_rankings_by_season.pdf")
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
          geom_vline(aes(xintercept = min(rank[percentile==10])),
                     size = 0.5,
                     linetype=2) +
          geom_vline(aes(xintercept = max(rank[percentile==1])),
                     size = 0.5,
                     linetype=2) +
          facet_wrap_paginate(~ Phylum*Order*Family*Genus, ncol = 2, nrow = 3, page = i))
  
}

dev.off()


#Oil overlap site & season-----

#Taxa that are positively associated with oil and BJ regardless of season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==1 &
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


#Venn Diagram to help me visualize my hypotheses-----
library(venn)
library(ggplot2)
library(ggpolypath)

sb.mod.dom <- sb.mod %>%
  filter(rare_dom=="dominant")

venn.list <- list(oil.pos = sb.mod.dom$featureid[sb.mod.dom$percentile_PAHs==10],
                  oil.neg = sb.mod.dom$featureid[sb.mod.dom$percentile_PAHs==1],
                  BJ = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Fourchon = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10])

venn(venn.list,  zcolor = "style")

venn.list <- list(oil.pos = sb.mod.dom$featureid[sb.mod.dom$percentile_PAHs==10],
                  oil.neg = sb.mod.dom$featureid[sb.mod.dom$percentile_PAHs==1],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season==1],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season==10])

venn(venn.list,  zcolor = "style")

venn.list <- list(BJ = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Fourchon = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10],
                  Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season==1],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season==10])

venn(venn.list,  zcolor = "style")

venn.list <- list(BJ = sb.mod.dom$featureid[sb.mod.dom$percentile_site==1],
                  Fourchon = sb.mod.dom$featureid[sb.mod.dom$percentile_site==10],
                  Season_NULL = sb.mod.dom$featureid[sb.mod.dom$percentile_season %in% c(5,6)])

venn(venn.list,  zcolor = "style")

venn.list <- list(Winter = sb.mod.dom$featureid[sb.mod.dom$percentile_season==1],
                  Summer = sb.mod.dom$featureid[sb.mod.dom$percentile_season==10],
                  Site_NULL = sb.mod.dom$featureid[sb.mod.dom$percentile_site %in% c(5,6)])

venn(venn.list,  zcolor = "style")


#Oil taxa with site and season agnostic
venn.list <- list(Oil.Pos = sb.mod.dom %>%
                    filter(percentile_PAHs==10) %>%
                    select(featureid) %>%
                    unlist(),
                  Winter_site_NULL = sb.mod.dom %>%
                    filter(percentile_season==1, 
                           percentile_site %in% c(5,6)) %>%
                    select(featureid) %>%
                    unlist(),
                  Summer = sb.mod.dom %>%
                    filter(percentile_season==10, 
                           percentile_site %in% c(5,6)) %>%
                    select(featureid)%>%
                    unlist(),
                  BJ = sb.mod.dom %>%
                    filter(percentile_site %in% c(1), 
                           percentile_season %in% c(5,6)) %>%
                    select(featureid)%>%
                    unlist(),
                  Fourchon = sb.mod.dom %>%
                    filter(percentile_site %in% c(10), 
                           percentile_season %in% c(5,6)) %>%
                    select(featureid)%>%
                    unlist())


venn(venn.list,  zcolor = "style")


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



#Chi-square test on Poisson------

#https://rpubs.com/vidhya36/292819

#make counts of each order by percentile
Order.counts <- rename(count(sb.oil, Order, percentile), Freq = n)
Order.counts$percentile <- as.numeric(Order.counts$percentile)

#Find median rank by Order
median.rank <- sb.oil %>%
  group_by(Order) %>%
  summarise(Median = median(rank))

median.rank$max_rank <- max(sb.oil$rank)

median.rank

#What are the top ranking taxa?
sb.oil[sb.oil$rank>700,]

#Run Poisson Test----
Order.counts <- split(Order.counts, f = Order.counts$Order)

Ctest.results <- list()

#Ctest.results <- lapply(Order.counts, function(x){

i <- 8
  x <- Order.counts[[i]]
  N <- sum(x$Freq)
  RF <- x$Freq/N
  MEAN <- sum(RF * x$percentile)
  VAR <- (sum(x$percentile^2*x$Freq) - N*MEAN^2)/(N-1) # else use (MEAN+MEAN^2/R)
  DISP <- VAR/MEAN # over dispersion
  THETA <- 1/DISP
  R <- MEAN*THETA/(1-THETA) # MEAN = R(1-THETA)/THETA
  cbind(MEAN,VAR,DISP,THETA,R)
  
  E_poi = round(N * dpois(x$percentile, lambda=MEAN),5)
  chisq.test(rbind(x$Freq,E_poi))
  names(Order.counts[i])

  #Overdispered = Evenly dispersed
  #Underdispered = Clumped
  
#Calculate log ratios of Fourchon relative to BJ for Orders of interest

#make list of hydrocarbon degrading orders per Tatariw et al.
HC.orders <- c("Alteromonadales","Chromatiales","Desulfobacterales","Desulfuromonadales","Methylococcales","Oceanospirillales","Rhizobiales","Rhodobacterales","Thiotrichales","Desulfovibrionales", "Pseudomonadales")

# subset to match OTUs thrown out by Songbird

df <- subset_taxa(bac.2season.not_filtered, taxa_names(bac.2season.not_filtered) %in% sb$featureid)

# melt phyloseq data to data frame
mdf <- psmelt(df)

# compute taxa sum according to site
sampletype_abund <- dplyr::group_by(mdf, OTU, site) %>% 
  dplyr::summarise(abundance = sum(Abundance))

site.lograt <-  tidyr::spread(sampletype_abund, site, abundance)

#calculate log ratios
site.lograt$logratio <- log(site.lograt$F/site.lograt$BJ)

#order based on rankings of differentials
featureid.sorted <- sb$featureid[order(sb$Fourchon_rel_BJ)]

site.lograt.reordered <- site.lograt[match(featureid.sorted, site.lograt$OTU),]

#bind taxonomy with results
tax.df.reordered <- as.data.frame(tax_table(df)[match(featureid.sorted, rownames(tax_table(df))),])
site.lograt.reordered <- cbind(site.lograt.reordered, tax.df.reordered)

#make factor so it's not rearranged during plotting
site.lograt.reordered$OTU <- factor(site.lograt.reordered$OTU, levels = site.lograt.reordered$OTU)

#Get rid of all Order names except Orders of interest
Order.others <- unique(site.lograt.reordered$Order)[!(unique(site.lograt.reordered$Order) %in% c(HC.orders))]
site.lograt.reordered$Order <- fct_collapse(site.lograt.reordered$Order, 
                             Other = Order.others)

site.lograt.reordered$Order <- fct_explicit_na(site.lograt.reordered$Order, na_level = "Other")

#plot
ggplot(data = site.lograt.reordered,
       aes(x = OTU,
           y = logratio,
           color = Order)) +
  geom_col(width = 1) +
  theme(legend.position = "bottom",
        axis.text = element_blank()) +
  labs(x = NULL) +
  facet_wrap(~ Order)


#boxplot
library(grDevices)

ggplot(data = site.lograt.reordered,
       aes(x = Order,
           y = logratio)) +
  geom_boxplot(aes(color = Order)) +
  labs(y = expression(paste("Natural Log-Ratio ", frac("Lightly Oiled", "Heavily Oiled"))))
  


#Examples for Songbird Github-----

#HC Degrading positive
ggplot(data = sb.mod %>%
                         filter(Phylum=="Proteobacteria", 
                                rare_dom=="dominant"),
                       aes(x = rank_PAHs,
                           y = scaled_PAHs,
                           fill = Class,
                           color = Class)) +
  geom_col(width = 1) +
  theme_bw() +
  labs(y = "Chemical Differential",
       x = "Rank")

ggsave("sb_example_no_facet.png")

ggplot(data = sb.mod %>%
                        filter(Phylum=="Proteobacteria", 
                               rare_dom=="dominant"),
                      aes(x = rank_PAHs,
                          y = scaled_PAHs,
                          fill = Class,
                          color = Class)) +
  geom_col(width = 1) +
  theme_bw() +
  labs(y = "Chemical Differential",
       x = "Rank") +
  facet_wrap(~ Class)

ggsave("sb_example_with_facet.png")

ggplot(data = sb.mod %>%
                      filter(Phylum=="Proteobacteria", 
                             rare_dom=="dominant"),
                    aes(x = Class,
                        y = rank_PAHs,
                        fill = Class,
                        color = Class)) +
  geom_jitter() +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  labs(y = "Rank",
       x = "Class")

ggsave("sb_example_stripchart.png")

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

#oil
ggplot(data = sb.mod %>%
         filter(rare_dom=="dominant"),
       aes(x = fct_reorder(Phylum, rank_PAHs, .fun='median'), 
           y = rank_PAHs,
           fill = Phylum,
           color = Phylum)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_boxplot(alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  labs(y = "Rank",
       x = "Phylum") +
  geom_hline(aes(yintercept = min(rank_PAHs[percentile_PAHs==10])),
             size = 0.5,
             linetype=2) +
  geom_hline(aes(yintercept = max(rank_PAHs[percentile_PAHs==1])),
             size = 0.5,
             linetype=2)

#ggsave("sb_example_stripchart.png")

#Taxa Overlap----

#Taxa that are positively associated with oil and BJ
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==1)

#Taxa that are positively associated with oil and F
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_site==10)

#Taxa that are negatively associated with oil and positive with BJ
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           percentile_site==1)

#Taxa that are negatively associated with oil and positive with F
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==1 &
           percentile_site==10)

#What is the overlap of the oil and the two sites together?  What is common to the two sites?

#Can't answer this because oil may be the reason for the differences between the two sites.  And It's impossible for the sites to share anything, except in the sense of the middle 50th percentile.  That is, the sites had differential oiling, so we can't assume the oiling was a controlled event like winter or summer.

#Taxa that are positively associated with oil and Summer
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs==10 &
           percentile_season==10)

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

#Log ratio LM of marinobacter----

sb.mod %>%
  filter(Genus=="Marinobacter")

oil.ass.list <- sb.mod %>%
  filter(featureid=="777466") %>%
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
qqnorm(resid(mod))
abline(0,1)

#Log ratio LM of HC degraders----

#Rhizobiales has stronger relationship than Alteromonadales

sb.mod %>%
  filter(Order=="Desulfobacterales")

oil.ass.list <- sb.mod %>%
  filter(rare_dom=="dominant" & Order=="Desulfobacterales" & percentile_PAHs==10) %>%
  rownames()

oil.noass.list <- sb.mod %>%
  filter(Order=="Desulfobacterales" & percentile_PAHs %in% c(5,6)) %>%
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

#Marinobacter

sb.mod %>%
  filter(Genus=="Marinobacter")

oil.ass.list <- sb.mod %>%
  filter(rare_dom=="dominant" & Genus=="Marinobacter" & percentile_PAHs==10) %>%
  rownames()

oil.noass.list <- sb.mod %>%
  filter(Order=="Desulfobacterales" & percentile_PAHs %in% c(5,6)) %>%
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

#No site or season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(5,6))

#Fungi-----
         
         
#load Songbird output-----
         
         sb <- read.delim("songbird/logdir/new_final_2/fungi/S2_MG_fungi_site_season_oil/differentials.tsv")
         
         #rename columns
         names(sb) <- c("featureid", "Intercept", "Fourchon_rel_BJ", "Summer_rel_Wint", "scaled_PAHs")
         
         #Using not filtered because the biom used in songbird wasn't filtered this way 
         fung.2season.sub <- subset_taxa(fung.2season_with_outliers, rownames(tax_table(fung.2season_with_outliers)) %in% sb$featureid)
         
         #For fun subset to filtered version
         #fung.2season.sub.filtered <- subset_taxa(fung.2season.sub, rownames(tax_table(fung.2season.sub)) %in% rownames(tax_table(fung.2season)))
         
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
fid <- prune_taxa(sb.mod$featureid, fung.2season_with_outliers) 
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
  mutate(freq = n / sum(n))

sb.mod %>%
  group_by(Phylum, rare_dom) %>%
  summarise (n = sum(tax.sums)) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#Taxa Overlap----

#F no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(10) &
           percentile_season %in% c(5,6))

#BJ no season
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(1) &
           percentile_season %in% c(5,6))

#Winter no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(1))

#Summer no site
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_site %in% c(5,6) &
           percentile_season %in% c(10))

#PAH associated fungi----

#positive
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs %in% c(10)) %>%
  select(Genus)

#negative
sb.mod %>%
  filter(rare_dom=="dominant" &
           percentile_PAHs %in% c(1))


#Taxa that are positively associated with oil and BJ regardless of season
         sb.mod %>%
           filter(rare_dom=="dominant" &
                    percentile_PAHs==10 &
                    percentile_site==1 &
                    percentile_season %in% c(5,6))
         
#Taxa that are positively associated with oil and F regardless of season
         sb.mod %>%
           filter(rare_dom=="dominant" &
                    percentile_PAHs==10 &
                    percentile_site==10 &
                    percentile_season %in% c(5,6))
         
#Taxa that are positively associated with oil and Summer, regardless of site
         sb.mod %>%
           filter(rare_dom=="dominant" &
                    percentile_PAHs==10 &
                    percentile_season==10 &
                    percentile_site %in% c(5,6))

#Taxa that are positively associated with oil and winter, regardless of site
         sb.mod %>%
           filter(rare_dom=="dominant" &
                    percentile_PAHs==10 &
                    percentile_season==10 &
                    percentile_site %in% c(5,6))

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

