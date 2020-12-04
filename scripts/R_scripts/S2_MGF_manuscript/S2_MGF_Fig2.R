#Script for S2_MG Manuscript

#Figure 2 - rarefaction curves of diversity and oil
#Description:   How many samples are necessary to legetimately describe salt marsh soil microbial communities?

#Last updated 30 Sep 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(cowplot)
library(tidyverse)
library(iNEXT)


# iNEXT: make incidence rarefaction curves with both bacteria and fungi unfiltered-----

#sampling effort for Shannon's Diversity

voi.list <- list(bac.2season.BJ.summer_with_outliers, 
                 bac.2season.BJ.winter_with_outliers, 
                 bac.2season.F.summer_with_outliers, 
                 bac.2season.F.winter_with_outliers, 
                 fung.2season.BJ.summer_with_outliers, 
                 fung.2season.BJ.winter_with_outliers, 
                 fung.2season.F.summer_with_outliers, 
                 fung.2season.F.winter_with_outliers)

#Add Bacteria or Fungi to each pseq object metadata

microbe.list <- c(rep("Bacteria", 4), rep("Fungi", 4))
for (i in seq_along(voi.list)) {
  sample_data(voi.list[[i]])$microbe.type <- microbe.list[i]
}

#Generate curve data

#make otu_table_list
voi.list.OT <- lapply(voi.list, function(x) {
  as.data.frame(otu_table(x))
  })

#add names
names(voi.list.OT) <- c("bac.2season.BJ.summer_with_outliers", 
                        "bac.2season.BJ.winter_with_outliers", 
                        "bac.2season.F.summer_with_outliers", 
                        "bac.2season.F.winter_with_outliers", 
                        "fung.2season.BJ.summer_with_outliers", 
                        "fung.2season.BJ.winter_with_outliers", 
                        "fung.2season.F.summer_with_outliers", 
                        "fung.2season.F.winter_with_outliers")


voi.list.OT <- lapply(voi.list.OT, function(x){
  x[x>0] <- 1
  return(x)  #important to remember this line of code
})

#Run iNEXT------

#Sys.time()
# voi.inext.df <- iNEXT(x = voi.list.OT, 
#                       q = c(0,1,2), 
#                       datatype = "incidence_raw", 
#                       size = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), 
#                       se = TRUE, 
#                       conf = 0.95, 
#                       nboot = 1000)  
# Sys.time() 

#took 25 min for 1000 bootstraps on new macbook pro. Only do once and reload below for editing

#write objects to rds to save time for the future
#write_rds(voi.inext.df, "images/manuscript/S2_MG_fig2_df.rds")

voi.inext.df <- read_rds("images/manuscript/fig2/S2_MG_fig2_df.rds")

#make into one list
voi.inext.all <- rbindlist(voi.inext.df$iNextEst)

#add names column
voi.inext.all$group <- rep(names(voi.inext.df$iNextEst), each = 90)

voi.inext.all$microbe.type <- c(rep("Bacteria", 90*4), rep("Fungi", 90*4))
voi.inext.all$site <- c(rep("Heavily Oiled", 90*2), rep("Lightly Oiled", 90*2),rep("Heavily Oiled", 90*2), rep("Lightly Oiled", 90*2))

voi.inext.all$season <- rep(c(rep("Summer", 90), rep("Winter", 90)), 4)

#grab max_cores confidence intervals
hline.div <- function(x) {
  x %>%
    filter(order == 1, method == "observed") %>%
    slice_max(t) %>%
    select(qD.LCL, qD.UCL)
}

hline.list <- lapply(voi.inext.df$iNextEst, hline.div)
hline.list.final <- rbindlist(hline.list)
hline.list.final$group <- c("bac.2season.BJ.summer", 
                        "bac.2season.BJ.winter", 
                        "bac.2season.F.summer", 
                        "bac.2season.F.winter", 
                        "fung.2season.BJ.summer", 
                        "fung.2season.BJ.winter", 
                        "fung.2season.F.summer", 
                        "fung.2season.F.winter")

hline.list.final$microbe.type <- c(rep("Bacteria", 4), rep("Fungi", 4))
hline.list.final$site <- c(rep("Heavily Oiled", 2), rep("Lightly Oiled", 2),rep("Heavily Oiled", 2), rep("Lightly Oiled", 2))

hline.list.final$season <- c(rep(c("Summer", "Winter"), 4))

#reorder season factors
hline.list.final$season <- factor(hline.list.final$season, levels = c("Winter", "Summer"))
voi.inext.all$season <- factor(voi.inext.all$season, levels = c("Winter", "Summer"))

#plot q = 1 only, not extrapolated points------
voi.inext.all %>%
  filter(order==1, method!="extrapolated") %>%
  ggplot(aes(x = t,
           y = qD,
           fill = season,
           shape = site)) +
  geom_point(size = 2, fill = "gray") +
  geom_rect(inherit.aes = FALSE,
            data = hline.list.final,
            aes(xmin = -Inf,
                xmax = Inf,
                ymin = qD.LCL,
                ymax = qD.UCL,
                fill = group),
            color = "black",
            fill = "gray",
            alpha=0.3) +
  geom_errorbar(aes(ymax = qD.UCL,
                    ymin = qD.LCL)) +
  facet_grid(rows = vars(microbe.type),
             cols = vars(season),  
             scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species (Hill order = 1)",
     shape = "Site") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(5, 8)) +
  guides(shape = guide_legend(override.aes = list(size = 5)))

#ggsave("images/manuscript/S2_MG_Fig2_v4.png", width = 8.5, height = 7, units = "in")

#calculate min samples for each group-----
voi.hill <- subset(voi.inext.all, method=="observed")

#split data frame into list
voi.inext.list <- split.data.frame(voi.inext.all, voi.inext.all$group)

minsamp.list.0 <- list()
minsamp.list.1 <- list()
minsamp.list.2 <- list()

voi.hill.0 <- subset(voi.hill, order==0)
voi.hill.1 <- subset(voi.hill, order==1)
voi.hill.2 <- subset(voi.hill, order==2)

for(i in 1:nrow(voi.hill.0)){
  
  minsamp.list.0[[i]] <- voi.inext.list[[voi.hill.0$group[i]]] %>%
    as.data.frame(value = names(voi.inext.list[[voi.hill.0$group[i]]])) %>%
    filter(order==0 & qD.UCL < voi.hill.0[i]$qD.LCL) %>%
    select(t, order, group) %>%
    slice(n())
}
  
for(i in 1:nrow(voi.hill.0)){
  
  minsamp.list.1[[i]] <- voi.inext.list[[voi.hill.1$group[i]]] %>%
    as.data.frame(value = names(voi.inext.list[[voi.hill.1$group[i]]])) %>%
    filter(order==1 & qD.UCL < voi.hill.1[i]$qD.LCL) %>%
    select(t, order, group) %>%
    slice(n())
}
 
for(i in 1:nrow(voi.hill.0)){
 minsamp.list.2[[i]] <- voi.inext.list[[voi.hill.2$group[i]]] %>%
    as.data.frame(value = names(voi.inext.list[[voi.hill.2$group[i]]])) %>%
    filter(order==2 & qD.UCL < voi.hill.2[i]$qD.LCL) %>%
    select(t, order, group) %>%
    slice(n())
}

#lists of threshold.  Remember that you need to add one to this because the code above asked to return the samples that had lower CI values
A <- rbindlist(minsamp.list.0)
B <- rbindlist(minsamp.list.1)
C <- rbindlist(minsamp.list.2)

#make into pretty table
voi.hill.table <- voi.hill %>%
  select(microbe.type,site,season, order,qD,qD.LCL, qD.UCL,t ) %>%
  arrange(order)


voi.hill.table$min.samples <- NA
voi.hill.table$min.samples[voi.hill.table$order==0] <- A$t
voi.hill.table$min.samples[voi.hill.table$order==1] <- B$t
voi.hill.table$min.samples[voi.hill.table$order==2] <- C$t
         
names(voi.hill.table) <- c("Microbe Type", "Site", "Season", "q (Hill order)", "Effective Number of Species", "Lower CI", "Upper CI","N", "Min. Samples")

write.csv(voi.hill.table, file = "S2_fig2_stats.csv", row.names = FALSE)

#plot all three orders, not extrapolated points------
voi.inext.all %>%
  filter(method!="extrapolated") %>%
  ggplot(aes(x = t,
           y = qD,
           fill = site,
           shape = factor(order))) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL, linetype = site),
               alpha = 0.5) +
  # geom_point(inherit.aes = FALSE, aes(x = t,
  #                                     y = qD),
  #            shape = 21, 
  #            fill = "white",
  #            color = "white",
  #            size = 2) +
  geom_point(size = 2) +
  facet_grid(rows = vars(microbe.type),
             cols = vars(season),  
             scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48,49,50)) +
  scale_fill_manual(values = c("black", "lightgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

#ggsave("images/manuscript/S2_MG_Fig2_v5.png", width = 8.5, height = 7, units = "in")

#plot same but faceted by site------

voi.inext.all %>%
  filter(method!="extrapolated") %>%
  ggplot(aes(x = t,
           y = qD,
           fill = season,
           shape = factor(order))) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL),
               alpha = 0.5) +
  # geom_point(inherit.aes = FALSE, aes(x = t,
  #                                     y = qD),
  #            shape = 21, 
  #            fill = "white",
  #            color = "white",
  #            size = 2) +
  geom_point(size = 2) +
  facet_grid(rows = vars(microbe.type),
             cols = vars(site),  
             scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48,49,50)) +
  scale_fill_manual(values = c("black", "lightgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

ggsave("images/manuscript/S2_MG_Fig2_v5_site.png", width = 8.5, height = 7, units = "in")


#Same comparison but with abundance-weighted Hill numbers------

#iNEXT just takes too darn long, so I'm going to compare observed estimates and see how it looks

p <- plot_richness(bac.2season_with_outliers)

#convert value to Hill numbers

p1 <- p$data

p1$value[which(p1$variable == "Observed")] <- p1$value[which(p1$variable == "Observed")]

p1$value[which(p1$variable == "Shannon")] <- exp(p1$value[which(p1$variable == "Shannon")])

p1$value[which(p1$variable == "Simpson")] <- 1/(1 - p1$value[which(p1$variable == "Simpson")])

p1$hill_order <- NA

p1$hill_order[which(p1$variable == "Observed")] <- 0
p1$hill_order[which(p1$variable == "Shannon")] <- 1
p1$hill_order[which(p1$variable == "Simpson")] <- 2

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

p1$season <- factor(p1$season, levels = c("WINTER", "SUMMER"))

library(Rmisc)

p1.summary <- summarySE(data = p1,
                        measurevar = "value",
                        groupvars = c("site", "season", "hill_order"))
#plot all three orders------

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = site,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ season) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = site,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_bac_unfiltered_abund.png", width = 8.5, height = 7, units = "in")

#plot same by site-----

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = season,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ site) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = season,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_bac_unfiltered_abund_site.png", width = 8.5, height = 7, units = "in")

#fungi-----
p <- plot_richness(fung.2season_with_outliers)

#convert value to Hill numbers

p1 <- p$data

p1$value[which(p1$variable == "Observed")] <- p1$value[which(p1$variable == "Observed")]

p1$value[which(p1$variable == "Shannon")] <- exp(p1$value[which(p1$variable == "Shannon")])

p1$value[which(p1$variable == "Simpson")] <- 1/(1 - p1$value[which(p1$variable == "Simpson")])

p1$hill_order <- NA

p1$hill_order[which(p1$variable == "Observed")] <- 0
p1$hill_order[which(p1$variable == "Shannon")] <- 1
p1$hill_order[which(p1$variable == "Simpson")] <- 2

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

p1$season <- factor(p1$season, levels = c("WINTER", "SUMMER"))

library(Rmisc)

p1.summary <- summarySE(data = p1,
                        measurevar = "value",
                        groupvars = c("site", "season", "hill_order"))
#plot all three orders------

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = site,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ season) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = site,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_fung_unfiltered_abund.png", width = 8.5, height = 7, units = "in")

#plot same by site-----

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = season,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ site) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = season,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_fung_unfiltered_abund_site.png", width = 8.5, height = 7, units = "in")

#Same but with filtered species------

# Subset "2season" by Site using filtered taxa------------

bac.2season.BJ_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#Subset to "2season" by WINTER-------

bac.2season.BJ.winter_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.winter_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.winter_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("WINTER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Subset to "2season" by SUMMER-------

bac.2season.BJ.summer_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

bac.2season.F.summer_with_outliers_filtered <- bac.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.summer_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("BJ") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer_with_outliers_filtered <- fung.2season_with_outliers_filtered %>%
  subset_samples( . , site %in% c("F") & season %in% c("SUMMER")) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

#sampling effort for Shannon's Diversity

voi.list <- list(bac.2season.BJ.summer_with_outliers_filtered, 
                 bac.2season.BJ.winter_with_outliers_filtered, 
                 bac.2season.F.summer_with_outliers_filtered, 
                 bac.2season.F.winter_with_outliers_filtered, 
                 fung.2season.BJ.summer_with_outliers_filtered, 
                 fung.2season.BJ.winter_with_outliers_filtered, 
                 fung.2season.F.summer_with_outliers_filtered, 
                 fung.2season.F.winter_with_outliers_filtered)

#Add Bacteria or Fungi to each pseq object metadata

microbe.list <- c(rep("Bacteria", 4), rep("Fungi", 4))
for (i in seq_along(voi.list)) {
  sample_data(voi.list[[i]])$microbe.type <- microbe.list[i]
}

#Generate curve data

#make otu_table_list
voi.list.OT <- lapply(voi.list, function(x) {
  as.data.frame(otu_table(x))
  })

#add names
names(voi.list.OT) <- c("bac.2season.BJ.summer_with_outliers_filtered", 
                        "bac.2season.BJ.winter_with_outliers_filtered", 
                        "bac.2season.F.summer_with_outliers_filtered", 
                        "bac.2season.F.winter_with_outliers_filtered", 
                        "fung.2season.BJ.summer_with_outliers_filtered", 
                        "fung.2season.BJ.winter_with_outliers_filtered", 
                        "fung.2season.F.summer_with_outliers_filtered", 
                        "fung.2season.F.winter_with_outliers_filtered")


voi.list.OT <- lapply(voi.list.OT, function(x){
  x[x>0] <- 1
  return(x)  #important to remember this line of code
})

#Run iNEXT------

Sys.time()
 voi.inext.df <- iNEXT(x = voi.list.OT, 
                      q = c(0,1,2),
                      datatype = "incidence_raw",
                      size = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),
                      se = TRUE,
                      conf = 0.95,
                      nboot = 1000)
Sys.time()

#took 2 min for 1000 bootstraps on new macbook pro. Clearly the hyperdiversity of the unfiltered set makes a big difference

#make into one list
voi.inext.all <- rbindlist(voi.inext.df$iNextEst)

#add names column
voi.inext.all$group <- rep(names(voi.inext.df$iNextEst), each = 90)

voi.inext.all$microbe.type <- c(rep("Bacteria", 90*4), rep("Fungi", 90*4))
voi.inext.all$site <- c(rep("Heavily Oiled", 90*2), rep("Lightly Oiled", 90*2),rep("Heavily Oiled", 90*2), rep("Lightly Oiled", 90*2))

voi.inext.all$season <- rep(c(rep("Summer", 90), rep("Winter", 90)), 4)

hline.div <- function(x) {
  x %>%
    filter(order == 1, method == "observed") %>%
    slice_max(t) %>%
    select(qD.LCL, qD.UCL)
}

hline.list <- lapply(voi.inext.df$iNextEst, hline.div)
hline.list.final <- rbindlist(hline.list)
hline.list.final$group <- c("bac.2season.BJ.summer_filtered", 
                        "bac.2season.BJ.winter_filtered", 
                        "bac.2season.F.summer_filtered", 
                        "bac.2season.F.winter_filtered", 
                        "fung.2season.BJ.summer_filtered", 
                        "fung.2season.BJ.winter_filtered", 
                        "fung.2season.F.summer_filtered", 
                        "fung.2season.F.winter_filtered")

hline.list.final$microbe.type <- c(rep("Bacteria", 4), rep("Fungi", 4))
hline.list.final$site <- c(rep("Heavily Oiled", 2), rep("Lightly Oiled", 2),rep("Heavily Oiled", 2), rep("Lightly Oiled", 2))

hline.list.final$season <- c(rep(c("Summer", "Winter"), 4))

#reorder season factors
hline.list.final$season <- factor(hline.list.final$season, levels = c("Winter", "Summer"))
voi.inext.all$season <- factor(voi.inext.all$season, levels = c("Winter", "Summer"))

#plot all three orders, not extrapolated points------

voi.inext.all %>%
  filter(method!="extrapolated") %>%
  ggplot(aes(x = t,
           y = qD,
           fill = site,
           shape = factor(order))) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL, linetype = site),
               alpha = 0.5) +
  # geom_point(inherit.aes = FALSE, aes(x = t,
  #                                     y = qD),
  #            shape = 21, 
  #            fill = "white",
  #            color = "white",
  #            size = 2) +
  geom_point(size = 2) +
  facet_grid(rows = vars(microbe.type),
             cols = vars(season),  
             scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48,49,50)) +
  scale_fill_manual(values = c("black", "lightgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

ggsave("images/manuscript/S2_MG_Fig2_v5_filtered.png", width = 8.5, height = 7, units = "in")

voi.inext.all %>%
  filter(method!="extrapolated") %>%
  ggplot(aes(x = t,
           y = qD,
           fill = season,
           shape = factor(order))) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL, linetype = site),
               alpha = 0.5) +
  # geom_point(inherit.aes = FALSE, aes(x = t,
  #                                     y = qD),
  #            shape = 21, 
  #            fill = "white",
  #            color = "white",
  #            size = 2) +
  geom_point(size = 2) +
  facet_grid(rows = vars(microbe.type),
             cols = vars(site),  
             scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48,49,50)) +
  scale_fill_manual(values = c("black", "lightgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

ggsave("images/manuscript/S2_MG_Fig2_v5_filtered_site.png", width = 8.5, height = 7, units = "in")

#calculate min samples for each group-----
voi.hill <- subset(voi.inext.all, method=="observed")

#split data frame into list
voi.inext.list <- split.data.frame(voi.inext.all, voi.inext.all$group)

minsamp.list.0 <- list()
minsamp.list.1 <- list()
minsamp.list.2 <- list()

voi.hill.0 <- subset(voi.hill, order==0)
voi.hill.1 <- subset(voi.hill, order==1)
voi.hill.2 <- subset(voi.hill, order==2)

for(i in 1:nrow(voi.hill.0)){
  
  minsamp.list.0[[i]] <- voi.inext.list[[voi.hill.0$group[i]]] %>%
    as.data.frame(value = names(voi.inext.list[[voi.hill.0$group[i]]])) %>%
    filter(order==0 & qD.UCL < voi.hill.0[i]$qD.LCL) %>%
    select(t, order, group) %>%
    slice(n())
}
  
for(i in 1:nrow(voi.hill.0)){
  
  minsamp.list.1[[i]] <- voi.inext.list[[voi.hill.1$group[i]]] %>%
    as.data.frame(value = names(voi.inext.list[[voi.hill.1$group[i]]])) %>%
    filter(order==1 & qD.UCL < voi.hill.1[i]$qD.LCL) %>%
    select(t, order, group) %>%
    slice(n())
}
 
for(i in 1:nrow(voi.hill.0)){
 minsamp.list.2[[i]] <- voi.inext.list[[voi.hill.2$group[i]]] %>%
    as.data.frame(value = names(voi.inext.list[[voi.hill.2$group[i]]])) %>%
    filter(order==2 & qD.UCL < voi.hill.2[i]$qD.LCL) %>%
    select(t, order, group) %>%
    slice(n())
}

#lists of threshold.  Remember that you need to add one to this because the code above asked to return the samples that had lower CI values
A <- rbindlist(minsamp.list.0)
B <- rbindlist(minsamp.list.1)
C <- rbindlist(minsamp.list.2)

#make into pretty table
voi.hill.table <- voi.hill %>%
  select(microbe.type,site,season, order,qD,qD.LCL, qD.UCL,t ) %>%
  arrange(order)


voi.hill.table$min.samples <- NA
voi.hill.table$min.samples[voi.hill.table$order==0] <- A$t
voi.hill.table$min.samples[voi.hill.table$order==1] <- B$t
voi.hill.table$min.samples[voi.hill.table$order==2] <- C$t
         
names(voi.hill.table) <- c("Microbe Type", "Site", "Season", "q (Hill order)", "Effective Number of Species", "Lower CI", "Upper CI","N", "Min. Samples")

#write.csv(voi.hill.table, file = "S2_fig2_filtered_stats.csv", row.names = FALSE)

#Same comparison but with abundance-weighted Hill numbers------

#iNEXT just takes too darn long, so I'm going to compare observed estimates and see how it looks

p <- plot_richness(bac.2season_with_outliers_filtered)

#convert value to Hill numbers

p1 <- p$data

p1$value[which(p1$variable == "Observed")] <- p1$value[which(p1$variable == "Observed")]

p1$value[which(p1$variable == "Shannon")] <- exp(p1$value[which(p1$variable == "Shannon")])

p1$value[which(p1$variable == "Simpson")] <- 1/(1 - p1$value[which(p1$variable == "Simpson")])

p1$hill_order <- NA

p1$hill_order[which(p1$variable == "Observed")] <- 0
p1$hill_order[which(p1$variable == "Shannon")] <- 1
p1$hill_order[which(p1$variable == "Simpson")] <- 2

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

p1$season <- factor(p1$season, levels = c("WINTER", "SUMMER"))

library(Rmisc)

p1.summary <- summarySE(data = p1,
                        measurevar = "value",
                        groupvars = c("site", "season", "hill_order"))
#plot all three orders------

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = site,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ season) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = site,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_bac_filtered_abund.png", width = 8.5, height = 7, units = "in")

#plot same by site-----

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = season,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ site) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = season,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_bac_filtered_abund_site.png", width = 8.5, height = 7, units = "in")

#fungi-----

p <- plot_richness(fung.2season_with_outliers_filtered)

#convert value to Hill numbers

p1 <- p$data

p1$value[which(p1$variable == "Observed")] <- p1$value[which(p1$variable == "Observed")]

p1$value[which(p1$variable == "Shannon")] <- exp(p1$value[which(p1$variable == "Shannon")])

p1$value[which(p1$variable == "Simpson")] <- 1/(1 - p1$value[which(p1$variable == "Simpson")])

p1$hill_order <- NA

p1$hill_order[which(p1$variable == "Observed")] <- 0
p1$hill_order[which(p1$variable == "Shannon")] <- 1
p1$hill_order[which(p1$variable == "Simpson")] <- 2

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

p1$season <- factor(p1$season, levels = c("WINTER", "SUMMER"))

library(Rmisc)

p1.summary <- summarySE(data = p1,
                        measurevar = "value",
                        groupvars = c("site", "season", "hill_order"))
#plot all three orders------

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = site,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ season) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = site,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


ggsave("images/manuscript/S2_MG_Fig2_v5_fung_filtered_abund.png", width = 8.5, height = 7, units = "in")

#plot same by site----

p1.summary %>%
  filter(hill_order %in% c(0,1,2)) %>%
  ggplot(aes(x = season,
           y = value,
           shape = factor(hill_order))) +
  facet_wrap(~ site) +
  geom_errorbar(aes(ymax = value + ci,
                    ymin = value-ci),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = season,
                 y = value),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,49,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))


#ggsave("images/manuscript/S2_MG_Fig2_v5_fung_filtered_abund_site.png", width = 8.5, height = 7, units = "in")
