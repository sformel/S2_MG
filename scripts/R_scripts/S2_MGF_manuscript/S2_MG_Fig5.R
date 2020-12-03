#Script for S2_MG Manuscript

#Figure 3 - beta diversity curves 
#Description:   NMDS of bacteria and fungi

subtitle <- "S2_MG_Fig5.R"

#Last updated 2 Sep 2020 by Steve Formel

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(tidyverse)
library(cowplot)

#Make NMDS----

voi <- bac.2season_with_outliers

bac.ord <- ordinate(voi, method = "NMDS", distance = "bray")

#put together data
plot.df <- data.frame(sample_data(voi), bac.ord$points)

#relabel levels
plot.df$site <- as.factor(plot.df$site)
levels(plot.df$site) <- plyr:::revalue(levels(plot.df$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

plot.df$season <- as.factor(plot.df$season)
levels(plot.df$season) <- plyr:::revalue(levels(plot.df$season), c("SUMMER" = "Summer", "WINTER" = "Winter"))

#plot
bac.NMDS.plot <- ggplot(data = plot.df,
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
  labs(subtitle = "Bacteria",
       x = "NMDS1",
       y = "NMDS2",
       color = "Season",
       shape = "Site",
       caption = paste0("Stress = ", round(bac.ord$stress, 2))) +
  theme(legend.text = element_text(size = 10),
        legend.position = "right") +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5, shape=21)))
  
bac.NMDS.plot 

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

#final plot----
top <- plot_grid(bac.NMDS.plot + theme(legend.position = "none"), 
                 fung.NMDS.plot + theme(legend.position = "none"), nrow = 1,
                 labels = c("A", "B"))

fig5 <- plot_grid(top, 
                   get_legend(bac.NMDS.plot + theme(legend.position = "right")), 
                   nrow = 1, 
                   rel_widths = c(1,0.2))

fig5

ggsave("images/manuscript/S2_MG_Fig5_v1.png", width = 10, height = 6, units = "in")

#PERMANOVA-----

voi <- bac.2season_with_outliers
#voi <- bac.2season_with_outliers_filtered

df <- voi %>%
  otu_table() %>%
  data.frame() %>%
  t()

df.env <- voi %>%
  sample_data() %>%
  data.frame()

#https://www.researchgate.net/post/How_do_I_know_how_many_permutations_to_use

set.seed(1)

bac.adonis.BC <- adonis(formula = df ~ site*season, data = df.env, method = "bray", permutations = 9999)
bac.adonis.JC <- adonis(formula = df ~ site*season, data = df.env, method = "jaccard", permutations = 9999, binary = TRUE)

#PERMDISP
A <- vegdist(x = df, method = "bray")
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

A <- vegdist(x = df, method = "jaccard", binary = TRUE)
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))


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
fung.adonis.JC <- adonis(formula = df ~ site*season, data = df.env, method = "jaccard", permutations = 9999, binary = TRUE)

#PERMDISP
A <- vegdist(x = df, method = "bray")
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

A <- vegdist(x = df, method = "jaccard", binary = TRUE)
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))






