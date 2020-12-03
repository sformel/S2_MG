#Script for S2_MG Manuscript

#Figure 4 - beta diversity curves for fungi only

#Last updated 15 Oct. 2020

#Description:   How many samples are necessary to legetimately describe salt marsh soil microbial communities?

#MultSE explained: 

#https://www.biorxiv.org/content/10.1101/2020.03.19.996991v1.full.pdf
#https://jonlefcheck.net/2015/03/31/how-much-is-enough-a-new-technique-for-quantifying-precision-of-community-surveys/

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(tidyverse)
library(cowplot)
#library(metagMisc)  #investigate this package more! https://rdrr.io/github/vmikk/metagMisc/

#Custom function----

#https://jonlefcheck.net/2015/03/31/how-much-is-enough-a-new-technique-for-quantifying-precision-of-community-surveys/

#Load multSE function
library(devtools)
multSE <- source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/multSE.R")[[1]]

Sys.time()
voi <- fung.2season_with_outliers
voi.dm <- distance(voi, method = "bray")
out.multiSE.fungi <- multSE(mat = voi.dm, 
                      group = sample_data(voi)$site_season, 
                      nresamp = 10000, 
                      permanova = FALSE)
Sys.time()

#took about 20 min each one

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
  scale_fill_manual(values = c("lightgray", "black")) +
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


#plot both together
library(cowplot)

plot_grid(multSE.plot + theme(legend.position = "none"), fung.NMDS.plot, ncol = 1, labels = "AUTO")

ggsave("images/manuscript/S2_MG_fungi_Fig3_v1.png", width = 7, height = 8, units = "in" )

#Run again without PERMANOVA so I get group specific thresholds

#Load multSE function
library(devtools)
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

fig4 <- rbind(minsamp.filtered, minsamp.unfiltered)

write.csv(fig4, file = "images/manuscript/fig3_fungi_minsamp_table.csv", row.names = FALSE)
