#Script for S2_MG Manuscript

#Figure 1 stats:  model of PAHs between site and season
#Description:   Do our results match our experimental strategy? Are there gradients of oil within each site and a large difference between sites?

#Last updated 10 Nov 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(compositions)
library(ggplot)
library(cowplot)


#color palette

cPAL <- c("#E69F00", "#0072B2")

#notes-----

#CLR Transformations------

voi <- fung.2season_with_outliers

df <- data.frame(sample_data(voi)) %>%
  select(SampleID,
         site,
         season,
         Total.C1.C3.Chrysenes,
         Total.C1.C4.naphthalenes,
         Total.C1.C4.phenanthrenes,
         Total.C1.C3.dibenzothiophenes,
         Total.relevant.PAHs)

#relabel levels
df$site <- as.factor(df$site)
levels(df$site) <- plyr:::revalue(levels(df$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

df$season <- as.factor(df$season)
levels(df$season) <- plyr:::revalue(levels(df$season), c("WINTER" = "Winter", "SUMMER" = "Summer"))

#reorder
levels(df$season) <- c("Winter", "Summer")

#http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.876.3979&rep=rep1&type=pdf

#good explanations and tutorial
#https://cran.r-project.org/web/packages/compositions/vignettes/UsingCompositions.pdf

#This transforms each number into it's relative abundance in that sample
comps <- acomp(df[,c(4:7)])

clr.df <- cbind(comps, df)

clr.df.gathered <- clr.df %>%
  gather(key = "PAH", value = "Proportion", c(1:4))

clr.df.gathered %>%
  ggplot(aes(x = PAH,
             y = Proportion,
             color = PAH)) +
  geom_boxplot() +
  facet_grid(rows = vars(site),
             cols = vars(season))
  
#But as I understand it, it's inappropriate to do the stats on this because it's closed.  Need to translate to CLR.

# 7 Multivariate Methods
# The central idea of the package ??? following the coordinate approach of [Pawlowsky-Glahn(2003)]
# and [Pawlowsky-Glahn and Mateu-Figueras(2005)] ??? is to transform the data by one
# of transforms into a classical multivariate dataset, to apply classical multivariate
# statistics and to back transform or interpreted the results afterwards in the original
# space

clr.PAH <- clr(na.omit(df[,c(4:7)]))
pc <- princomp(x = clr.PAH)

pc$Loadings # The loadings as compositional vector
pc$loadings # The loadings in clr-space
df.pca <- pc$scores

aitchison.pca.plot <- cbind(na.omit(df), df.pca) %>%
  ggplot(aes(x = Comp.1,
             y = Comp.2,
             fill = season,
             shape = site)) +
  geom_point(size = 4,
             stroke = 1,
             alpha = 0.5,
             color = "black") +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  labs(fill = "Season",
       shape = "Site") +
  theme_bw() +
  theme(legend.position = "bottom", legend.box="vertical") +
  guides(fill=guide_legend(override.aes=list(shape=21)))


#PERMANOVA----
clr.PAH <- clr(na.omit(df[,c(4:7)]))
clr.PAH <- as.data.frame(clr.PAH)

df.na <- na.omit(df) %>%
  select(SampleID, site, season)

df.clr <- cbind(df.na, clr.PAH)

set.seed(1)
adonis(formula = clr.PAH ~ site*season, data = df.na, permutations = 9999, method = "euclidean")

CLR.boxplot <- df.clr %>%
  gather(key = "PAH", value = "clr_val", Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
  ggplot(aes(x = PAH,
         y = clr_val),
         color = "black") +
  geom_boxplot(aes(fill = season),
         alpha = 0.5) +
  facet_grid(rows = vars(site), cols = vars(season)) +
  theme_bw() +
  scale_fill_manual(values = cPAL) +
  scale_x_discrete(labels = c("Chrysenes", "Dibenzothiophenes", "Naphthalenes", "Phenanthrenes")) +
  labs(y = "Centered Log-Ratio") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))


#plot Fig S1
plot_grid(aitchison.pca.plot, CLR.boxplot, labels = "AUTO")

ggsave("images/manuscript/S2_MGF_final/S2_MG_S1_PAH_comp.png", width = 10, height = 6, units = "in")
