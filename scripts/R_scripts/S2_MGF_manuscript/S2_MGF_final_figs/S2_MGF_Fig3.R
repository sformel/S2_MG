#Script for S2_MG Manuscript

#Figure 3 - beta diversity curves 
#Description:   multSE, NMDS and PERMANOVA

#Last updated 11 Nov 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(tidyverse)
library(cowplot)
library(devtools)
library(compositions)

#color palette

cPAL <- c("#E69F00", "#0072B2")

#Custom function----

#https://jonlefcheck.net/2015/03/31/how-much-is-enough-a-new-technique-for-quantifying-precision-of-community-surveys/

#Load multSE function

multSE <- source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/multSE.R")[[1]]

Sys.time()
voi <- fung.2season_with_outliers
voi.dm <- distance(voi, method = "bray")
out.multiSE.fungi <- multSE(mat = voi.dm, 
                      group = sample_data(voi)$site_season, 
                      nresamp = 10000, 
                      permanova = FALSE)
Sys.time()

#took about 1 min

out.multiSE.fungi$site <- c(rep("Heavily Oiled", 39), rep("Lightly Oiled", 46))
out.multiSE.fungi$season <- c(rep("Winter", 20), rep("Summer", 42), rep("Winter", 23))

# Plot output----

multSE.plot <- ggplot(out.multiSE.fungi, aes(x = n.samp, 
                                             y = means,
                                             shape = site,
                                             fill = season)) +
  geom_errorbar(aes(ymax = upper.ci, 
                    ymin = lower.ci), 
                width = 0.2) +
  geom_point() + 
  theme(legend.text = element_text(size = 10) ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = cPAL) +
  scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5)) +
  #guides(shape = guide_legend(override.aes = list(size = 5))) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  facet_grid(rows = vars(season),
             cols = vars(site)) +
  labs(x = "Sample size (n)", 
       y = "Multivariate pseudo SE")
  #theme(legend.position = c(0.8, 0.8), 
  #      panel.grid.minor = element_blank(),
  #      legend.box.background = element_rect(colour = "black")) 

#make NMDS-----

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
  scale_fill_manual(values = cPAL) + 
  theme_bw(base_size = 14) +
  theme_bw() +
  labs(x = "NMDS1",
       y = "NMDS2",
       color = "Season",
       shape = "Site",
       caption = paste0("Stress = ", round(fung.ord$stress, 2))) +
  theme(legend.text = element_text(size = 10),
        legend.position = "right") +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5, shape=21)))


#plot both together
plot_grid(multSE.plot + theme(legend.position = "none"), fung.NMDS.plot, ncol = 1, labels = "AUTO")

ggsave("images/manuscript/S2_MGF_final/S2_MGF_Fig3_v1.png", width = 7, height = 8, units = "in" )

#Run again without PERMANOVA so I get group specific thresholds

#Load multSE function
minsamp <- source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/minsamp.R")[[1]]

minsamp.fung <- minsamp(out.multiSE.fungi, group = sample_data(voi)$site_season)
minsamp.unfiltered <- minsamp.fung

minsamp.unfiltered$filtered <- "NO"

#Same with filtered samples----

#Run again without PERMANOVA so I get group specific thresholds

Sys.time()
voi <- fung.2season_with_outliers_filtered
voi.dm <- distance(voi, method = "bray")
out.multiSE.fungi <- multSE(mat = voi.dm, 
                            group = sample_data(voi)$site_season, 
                            nresamp = 10000, 
                            permanova = FALSE)
Sys.time()

minsamp.fung <- minsamp(out.multiSE.fungi, group = sample_data(voi)$site_season)
minsamp.filtered <- minsamp.fung

minsamp.filtered$filtered <- "YES"

multSE.table <- rbind(minsamp.filtered, minsamp.unfiltered)

write.csv(multSE.table, file = "images/manuscript/S2_MGF_final/S2_MGF_fig3_minsamp_table.csv", row.names = FALSE)


#PERMANOVA and Betadisp-----

set.seed(1)

voi <- fung.2season_with_outliers
#voi <- fung.2season_with_outliers_filtered

df.full <- voi %>%
  otu_table() %>%
  data.frame() %>%
  t()

df.env <- voi %>%
  sample_data() %>%
  data.frame()

df.env.na <- data.frame(sample_data(voi)) %>%
  select(SampleID,
         site,
         season,
         Total.relevant.PAHs) %>%
  na.omit()

df.full.na<- df.full[rownames(df.full) %in% rownames(df.env.na),]

#https://www.researchgate.net/post/How_do_I_know_how_many_permutations_to_use

set.seed(1)

fung.adonis.BC <- adonis(formula = df.full ~ site*season, data = df.env, method = "bray", permutations = 9999)

fung.adonis.BC <- adonis(formula = df.full.na ~ site*season + Total.relevant.PAHs, data = df.env.na, method = "bray", permutations = 9999)

fung.adonis.JC <- adonis(formula = df.full ~ site*season, data = df.env, method = "jaccard", permutations = 9999, binary = TRUE)

fung.adonis.JC <- adonis(formula = df.full ~ site*season + Total.relevant.PAHs, data = df.env.na, method = "jaccard", permutations = 9999, binary = TRUE)

#Aitchison
clr.df <- clr(df.full)
clr.df <- as.data.frame(clr.df)

clr.df.na <- clr(df.full.na) %>%
  as.data.frame()

fung.adonis.AT <- adonis(formula = clr.df ~ site*season, data = df.env, method = "euclidean", permutations = 9999)

fung.adonis.AT <- adonis(formula = clr.df.na ~ site*season + Total.relevant.PAHs, data = df.env.na, method = "euclidean", permutations = 9999)

#PERMDISP----
A <- vegdist(x = df.full, method = "bray")
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

A <- vegdist(x = df.full, method = "jaccard", binary = TRUE)
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

A <- vegdist(x = clr.df, method = "euclidean")
anova(betadisper(d = A, group = df.env$site))
anova(betadisper(d = A, group = df.env$season))

#dbRDA with oil------

fung.dbRDA.BC <- dbrda(formula = df.full.na ~ site*season + Total.relevant.PAHs, data = df.env.na, dist = "bray", permutations = 9999)

anova(fung.dbRDA.BC, by = "terms")
anova(fung.dbRDA.BC, by = "margin")

fung.dbRDA.JC <- dbrda(formula = df.full.na ~ site*season + Total.relevant.PAHs, data = df.env.na, dist = "jaccard", permutations = 9999, binary = TRUE)

anova(fung.dbRDA.JC, by = "margin")
anova(fung.dbRDA.JC, by = "terms")

fung.dbRDA.AT <- dbrda(formula = clr.df.na ~ site*season + Total.relevant.PAHs, data = df.env.na, dist = "euclidean", permutations = 9999)

anova(fung.dbRDA.AT, by = "margin")
anova(fung.dbRDA.AT, by = "terms")
