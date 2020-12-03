#Script for S2_MG Manuscript

#Figure 3 - beta diversity curves 
#Description:   NMDS of bacteria and fungi

#Last updated 16 Oct 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(tidyverse)
library(cowplot)

#Make NMDS----

#fungi  
voi <- fung.2season_with_outliers

fung.ord <- ordinate(voi, method = "NMDS", distance = "bray")

#put together data
plot.df <- data.frame(sample_data(voi), fung.ord$points)

#relabel levels
plot.df$site <- as.factor(plot.df$site)
levels(plot.df$site) <- plyr:::revalue(levels(plot.df$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

plot.df$season <- as.factor(plot.df$season)
levels(plot.df$season) <- plyr:::revalue(levels(plot.df$season), c("SUMMER" = "Summer", "WINTER" = "Winter"))

#plot
fung.NMDS.plot <- ggplot(data = plot.df,
                        aes(x = MDS1,
                            y = MDS2,
                            shape = site,
                          fill = season)) +
  geom_point(size = 3,
             stroke = 0.5,
             color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("lightgray", "black")) + 
  theme_bw(base_size = 14) +
  theme_bw() +
  labs(subtitle = "Fungi",
       x = "NMDS1",
       y = "NMDS2",
       color = "Season",
       shape = "Site",
       caption = paste0("Stress = ", round(fung.ord$stress, 2))) +
  theme(legend.text = element_text(size = 10),
        legend.position = "right") +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5, shape=21)))

fung.NMDS.plot

#PERMANOVA-----

set.seed(1)

voi <- fung.2season_with_outliers
#voi <- fung.2season_with_outliers_filtered

df <- voi %>%
  otu_table() %>%
  data.frame() %>%
  t()

df.env <- voi %>%
  sample_data() %>%
  data.frame()

#https://www.researchgate.net/post/How_do_I_know_how_many_permutations_to_use

set.seed(1)

fung.adonis.BC <- adonis(formula = df ~ site*season, data = df.env, method = "bray", permutations = 9999)

fung.adonis.BC <- adonis(formula = df ~ site*season + Total.relevant.PAHs, data = df.env, method = "bray", permutations = 9999)

fung.adonis.JC <- adonis(formula = df ~ site*season, data = df.env, method = "jaccard", permutations = 9999, binary = TRUE)

fung.adonis.JC <- adonis(formula = df ~ site*season + Total.relevant.PAHs, data = df.env, method = "jaccard", permutations = 9999, binary = TRUE)

#PERMDISP
A <- vegdist(x = df, method = "bray")
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

A <- vegdist(x = df, method = "jaccard", binary = TRUE)
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

#dbRDA with oil------

df <- voi %>%
  otu_table() %>%
  data.frame() %>%
  t()

df.env <- voi %>%
  sample_data() %>%
  data.frame() %>%
  select(SampleID,
         site,
         season,
         Total.C1.C3.Chrysenes,
         Total.C1.C4.naphthalenes,
         Total.C1.C4.phenanthrenes,
         Total.C1.C3.dibenzothiophenes,
         Total.relevant.PAHs) %>%
  na.omit()

df <- df[rownames(df) %in% df.env$SampleID, ]

fung.dbRDA.BC <- dbrda(formula = df ~ site*season + Total.relevant.PAHs, data = df.env, dist = "bray", permutations = 9999)

anova(fung.dbRDA.BC, by = "terms")
anova(fung.dbRDA.BC, by = "margin")

fung.dbRDA.JC <- dbrda(formula = df ~ site*season + Total.relevant.PAHs, data = df.env, dist = "jaccard", permutations = 9999, binary = TRUE)

anova(fung.dbRDA.JC, by = "margin")
anova(fung.dbRDA.JC, by = "terms")







