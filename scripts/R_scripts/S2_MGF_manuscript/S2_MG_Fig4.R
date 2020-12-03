#Script for S2_MG Manuscript

#Figure 4 - beta diversity curves 
#Description:   How many samples are necessary to legetimately describe salt marsh soil microbial communities?

#MultSE explained: 

#https://www.biorxiv.org/content/10.1101/2020.03.19.996991v1.full.pdf
#https://jonlefcheck.net/2015/03/31/how-much-is-enough-a-new-technique-for-quantifying-precision-of-community-surveys/

subtitle <- "S2_MG_Fig3.R"

#Last updated 11 Aug 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(tidyverse)
library(cowplot)
#library(metagMisc)  #investigate this package more! https://rdrr.io/github/vmikk/metagMisc/

#Add Bacteria or Fungi to each pseq object metadata

voi.list <- c(bac.2season_with_outliers, fung.2season_with_outliers)

microbe.list <- c(rep("Bacteria", 1), rep("Fungi", 1))
for (i in seq_along(voi.list)) {
  sample_data(voi.list[[i]])$microbe.type <- microbe.list[i]
}
#Custom function----

#https://jonlefcheck.net/2015/03/31/how-much-is-enough-a-new-technique-for-quantifying-precision-of-community-surveys/

#Load multSE function
library(devtools)
multSE <- source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/multSE.R")[[1]]

#source("multSE.R")  #saved on computer

#First run to make one curve for bacteria and one for fungi-----

Sys.time()
voi <- bac.2season_with_outliers
voi.dm <- distance(voi, method = "bray")
out.multiSE.bacteria <- multSE(mat = voi.dm, 
                      group = sample_data(voi)$site_season, 
                      nresamp = 10000, 
                      permanova = TRUE)
Sys.time()
voi <- fung.2season_with_outliers
voi.dm <- distance(voi, method = "bray")
out.multiSE.fungi <- multSE(mat = voi.dm, 
                      group = sample_data(voi)$site_season, 
                      nresamp = 10000, 
                      permanova = TRUE)
Sys.time()

#took about 20 min each one

out.multiSE.bacteria$microbe.type <- rep("Bacteria", nrow(out.multiSE.bacteria))
out.multiSE.fungi$microbe.type <- rep("Fungi", nrow(out.multiSE.fungi))

out.multiSE.all <- rbind(out.multiSE.bacteria, out.multiSE.fungi)

# Plot output----

both.multSE.plot <- ggplot(out.multiSE.all, aes(x = n.samp, y = means, group = microbe.type)) +
  geom_hline(aes(yintercept = upper.ci[18])) +
  geom_hline(aes(yintercept = upper.ci[36]), linetype = 2) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.2)+
  geom_point(aes(shape = microbe.type), fill = "gray", size = 2) + 
  labs(x = "Number of Soil Cores",
       y = "Effective Number of Species (q = 1)",
       shape = "Microbe Type") +
  theme(legend.text = element_text(size = 10) ) +
  scale_shape_manual(values = c(24, 25)) +
  scale_x_continuous(breaks = c(1:20)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 14) +
  labs(x = "Sample size (n)", y = "Multivariate pseudo SE") +
  theme(legend.position = c(0.8, 0.8), 
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(colour = "black")) 

ggsave("images/manuscript/S2_MG_Fig4_v1.png", width = 7, height = 6, units = "in")

#Run again without PERMANOVA so I get group specific thresholds

Sys.time()
voi <- bac.2season_with_outliers
voi.dm <- distance(voi, method = "bray")
out.multiSE.bacteria <- multSE(mat = voi.dm, 
                               group = sample_data(voi)$site_season, 
                               nresamp = 10000, 
                               permanova = FALSE)
Sys.time()
voi <- fung.2season_with_outliers
voi.dm <- distance(voi, method = "bray")
out.multiSE.fungi <- multSE(mat = voi.dm, 
                            group = sample_data(voi)$site_season, 
                            nresamp = 10000, 
                            permanova = FALSE)
Sys.time()

#Load multSE function
library(devtools)
minsamp <- source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/minsamp.R")[[1]]

minsamp.bac <- minsamp(out.multiSE.bacteria, group = sample_data(voi)$site_season)
minsamp.fung <- minsamp(out.multiSE.fungi, group = sample_data(voi)$site_season)

minsamp.bac$microbe <- "bacteria"
minsamp.fung$microbe <- "fungi"

minsamp.unfiltered <- rbind(minsamp.bac,minsamp.fung)

minsamp.unfiltered$filtered <- "NO"

#Same with filtered samples----

#Run again without PERMANOVA so I get group specific thresholds

Sys.time()
voi <- bac.2season_with_outliers_filtered
voi.dm <- distance(voi, method = "bray")
out.multiSE.bacteria <- multSE(mat = voi.dm, 
                               group = sample_data(voi)$site_season, 
                               nresamp = 10000, 
                               permanova = FALSE)
Sys.time()
voi <- fung.2season_with_outliers_filtered
voi.dm <- distance(voi, method = "bray")
out.multiSE.fungi <- multSE(mat = voi.dm, 
                            group = sample_data(voi)$site_season, 
                            nresamp = 10000, 
                            permanova = FALSE)
Sys.time()

minsamp.bac <- minsamp(out.multiSE.bacteria, group = sample_data(voi)$site_season)
minsamp.fung <- minsamp(out.multiSE.fungi, group = sample_data(voi)$site_season)

minsamp.bac$microbe <- "bacteria"
minsamp.fung$microbe <- "fungi"

minsamp.filtered <- rbind(minsamp.bac,minsamp.fung)

minsamp.filtered$filtered <- "YES"

fig4 <- rbind(minsamp.filtered, minsamp.unfiltered)

write.csv(fig4, file = "images/manuscript/fig4_minsamp_table.csv", row.names = FALSE)
