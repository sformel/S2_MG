#Script for S2_MG Manuscript

#Figure 2 - rarefaction curves of diversity and oil
#Description:   How many samples are necessary to legetimately describe salt marsh soil microbial communities?

#Last updated 15 Oct 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(cowplot)
library(tidyverse)
library(iNEXT)



# iNEXT: make abundance rarefaction curves with fungi unfiltered-----

#Run iNEXT abundance------

# melt phyloseq data to data frame
mdf <- psmelt(fung.2season_with_outliers)

# compute taxa sum according to sample type
sampletype_abund <- dplyr::group_by(mdf, OTU, site,season) %>% 
  dplyr::summarise(abundance = sum(Abundance))

df.taxasums <- sampletype_abund %>%
pivot_wider(names_from = c(site,season), names_sep = ".", values_from = c(abundance)) %>%
  as.data.frame() 

df.taxasums1 <- df.taxasums[,-1]
rownames(df.taxasums1) <- df.taxasums[,1]

#Sys.time()
#iNEXT(x = df.taxasums1, 
#                        q = c(0,2),
#                       datatype = "abundance",
#                       se = TRUE,
#                       conf = 0.95,
#                       nboot = 1000,
#                       endpoint = 600000)
# Sys.time()

#10 seconds

source("iNEXT_parallel.R")
Sys.time()
abund.out <- iNEXT_parallel(x = df.taxasums1, 
                       q = c(0,1,2),
                      datatype = "abundance",
                      se = TRUE,
                      conf = 0.95,
                      nboot = 1000,
                      endpoint = 600000)
Sys.time() 

#under a minute with nboot = 1000, endpoint = 1000
#about 4 minutes with nboot = 1000, endpoint = 10000
#about 27 minutes with nboot = 1000, endpoint = 100000
#about 45-60 minutes with nboot = 1000, endpoint = 600000

#Throws a warning that can be ignored.


#make into one list
## Extract results (list with data.frames for each sample)
res <- plyr::llply(.data = abund.out, .fun = function(z){ z$iNextEst })

voi.inext.all.abund <- rbindlist(res)

#add names column
voi.inext.all.abund$group <- rep(names(voi.inext.df.abund$iNextEst), each = 120)
voi.inext.all.abund$site <- c(rep("Heavily Oiled", 240), 
                        rep("Lightly Oiled", 240))

voi.inext.all.abund$season <- rep(c(rep("Summer", 120), rep("Winter", 120)), 2)

#adjust factors
voi.inext.all.abund$season <- factor(voi.inext.all.abund$season, levels = c("Winter", "Summer"))
voi.inext.all.abund$order <- factor(voi.inext.all.abund$order)

#plot------

#Just 0 and 2
voi.inext.all.abund <- voi.inext.all.abund %>%
  filter(order %in% c(0,2), method=="observed")

abund.plot <- voi.inext.all.abund %>%
  ggplot(aes(x = site,
           y = qD,
           shape = factor(order))) +
  facet_wrap(~ season) +
  geom_errorbar(aes(ymax = qD.UCL,
                    ymin = qD.LCL),
                width = 0.2) + 
  geom_point(inherit.aes = FALSE, 
             aes(x = site,
                 y = qD),
              shape = 21, 
              fill = "white",
              color = "black",
              size = 3,
             stroke = 0.5) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4))) +
  facet_grid(cols = vars(season),
              rows = vars(order),
              scales = "free") 

# abund.plot  <- voi.inext.all.abund %>%
#   ggplot(aes(x = m,
#            y = qD,
#            fill = site,
#            shape = order)) +
#   geom_ribbon(aes(ymax = qD.UCL,
#                     ymin = qD.LCL),
#                alpha = 0.5) +
#   geom_line(inherit.aes = FALSE, data = voi.inext.all.abund %>%
#               filter(method=="extrapolated"),
#              aes(x = m,
#                  y = qD,
#                  group = site),
#              color = "black",
#              size = 0.5,
#             linetype = "dashed") +
#   geom_point(inherit.aes = FALSE, data = voi.inext.all.abund %>%
#               filter(method!="extrapolated"), 
#              aes(x = m,
#                  y = qD,
#                  fill = site,
#                  shape = order),
#            size = 2, fill = "gray") +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank()) +
#   labs(x = "Number of Reads",
#      y = "Effective Number of Species",
#      shape = "Hill Order",
#      color = "Site") +
#   theme(legend.text = element_text(size = 10),
#         legend.position = "bottom") +
#   scale_shape_manual(values = c(48,50)) +
#   scale_fill_manual(values = c("black", "darkgray")) +
#   scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
#   guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4))) +
#   facet_grid(cols = vars(season),
#              rows = vars(order),
#              scales = "free") 

#Run iNEXT incidence----
voi.list <- list(fung.2season.BJ.summer_with_outliers, 
                 fung.2season.BJ.winter_with_outliers, 
                 fung.2season.F.summer_with_outliers, 
                 fung.2season.F.winter_with_outliers)

#make otu_table_list
voi.list.OT <- lapply(voi.list, function(x) {
  as.data.frame(otu_table(x))
  })

#add names
names(voi.list.OT) <- c("fung.2season.BJ.summer_with_outliers", 
                        "fung.2season.BJ.winter_with_outliers", 
                        "fung.2season.F.summer_with_outliers", 
                        "fung.2season.F.winter_with_outliers")


voi.list.OT <- lapply(voi.list.OT, function(x){
  x[x>0] <- 1
  return(x)  #important to remember this line of code
})


Sys.time()
voi.inext.df.inc <- iNEXT(x = voi.list.OT,
                      q = c(0,1,2),
                      datatype = "incidence_raw",
                      size = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),
                      se = TRUE,
                      conf = 0.95,
                      nboot = 1000)
Sys.time()

#took about 30 seconds

#make into one list
voi.inext.all <- rbindlist(voi.inext.df.inc$iNextEst)

#add names column
voi.inext.all$group <- rep(names(voi.inext.df.inc$iNextEst), each = 90)
voi.inext.all$site <- c(rep("Heavily Oiled", 180), 
                        rep("Lightly Oiled", 180))

voi.inext.all$season <- rep(c(rep("Summer", 90), rep("Winter", 90)), 2)

#adjust factors
voi.inext.all$season <- factor(voi.inext.all$season, levels = c("Winter", "Summer"))
voi.inext.all$order <- factor(voi.inext.all$order)

#plot------

voi.inext.all <- voi.inext.all %>%
  filter(order %in% c(0,2))

inc.plot <- voi.inext.all %>%
  ggplot(aes(x = t,
           y = qD,
           fill = site,
           shape = order)) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL),
               alpha = 0.5) +
  geom_line(inherit.aes = FALSE, data = voi.inext.all %>%
              filter(method=="extrapolated"),
             aes(x = t,
                 y = qD,
                 group = interaction(site, order)),
             color = "black",
             size = 0.5,
            linetype = "dashed")  +
  geom_point(inherit.aes = FALSE, data = voi.inext.all %>%
              filter(method!="extrapolated"), 
             aes(x = t,
                 y = qD,
                 fill = site,
                 shape = order),
           size = 2, fill = "gray") +
  facet_wrap(~ season) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species",
     shape = "Hill Order",
     fill = "Site") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48,50)) +
  scale_fill_manual(values = c("black", "darkgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

#plot of avg observed Hill orders by abundance-----

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
p1$season <- plyr:::revalue(p1$season, replace = c("WINTER" = "Winter", "SUMMER" = "Summer"))

library(Rmisc)

p1.summary <- summarySE(data = p1,
                        measurevar = "value",
                        groupvars = c("site", "season", "hill_order"))

#plot by site----

avg.abund.plot <- p1.summary %>%
  filter(hill_order %in% c(0,2)) %>%
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
              size = 3,
             stroke = 0.5) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Site",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  scale_shape_manual(values = c(48,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4))) +
  facet_grid(cols = vars(season),
              rows = vars(hill_order),
              scales = "free") 

#plot all 3 plots----

top <- plot_grid(avg.abund.plot + theme(legend.position = "none"), abund.plot + theme(legend.position = "none"), labels = c("A", "B"))
plot_grid(top, inc.plot  + theme(legend.position = "right"), nrow = 2, labels = c('','C'))

ggsave("images/manuscript/S2_MGF_Fig2_v2.png", width = 10, height = 8, units = "in")

#FILTERED-------

#Run iNEXT abundance------

# melt phyloseq data to data frame
mdf <- psmelt(fung.2season_with_outliers_filtered)

# compute taxa sum according to sample type
sampletype_abund <- dplyr::group_by(mdf, OTU, site,season) %>% 
  dplyr::summarise(abundance = sum(Abundance))

df.taxasums <- sampletype_abund %>%
pivot_wider(names_from = c(site,season), names_sep = ".", values_from = c(abundance)) %>%
  as.data.frame() 

df.taxasums1 <- df.taxasums[,-1]
rownames(df.taxasums1) <- df.taxasums[,1]

Sys.time()
 voi.inext.df.abund <- iNEXT(x = df.taxasums1, 
                       q = c(0,2),
                      datatype = "abundance",
                      se = TRUE,
                      conf = 0.95,
                      nboot = 100,
                      endpoint = 600000)
Sys.time()

#Took about 30 seconds on macbook pro

#make into one list
voi.inext.all.abund <- rbindlist(voi.inext.df.abund$iNextEst)

#add names column
voi.inext.all.abund$group <- rep(names(voi.inext.df.abund$iNextEst), each = 80)
voi.inext.all.abund$site <- c(rep("Heavily Oiled", 160), 
                        rep("Lightly Oiled", 160))

voi.inext.all.abund$season <- rep(c(rep("Summer", 80), rep("Winter", 80)), 2)

#adjust factors
voi.inext.all.abund$season <- factor(voi.inext.all.abund$season, levels = c("Winter", "Summer"))
voi.inext.all.abund$order <- factor(voi.inext.all.abund$order)

#plot------
abund <- voi.inext.all.abund %>%
  ggplot(aes(x = m,
           y = qD,
           fill = site,
           shape = order)) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL),
               alpha = 0.5) +
  geom_point(inherit.aes = FALSE, data = voi.inext.all.abund %>%
              filter(method=="extrapolated"),
             aes(x = m,
                 y = qD),
             shape = 21,
             fill = "white",
             color = "white",
             size = 3) +
  geom_point(size = 2, fill = "gray") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Number of Reads",
     y = "Effective Number of Species",
     shape = "Hill Order",
     color = "Site") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48, 50)) +
  scale_fill_manual(values = c("black", "lightgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4))) +
  facet_grid(cols = vars(season),
             rows = vars(order),
             scales = "free") 

#Run iNEXT incidence----
voi.list <- list(fung.2season.BJ.summer_with_outliers_filtered, 
                 fung.2season.BJ.winter_with_outliers_filtered, 
                 fung.2season.F.summer_with_outliers_filtered, 
                 fung.2season.F.winter_with_outliers_filtered)

#make otu_table_list
voi.list.OT <- lapply(voi.list, function(x) {
  as.data.frame(otu_table(x))
  })

#add names
names(voi.list.OT) <- c("fung.2season.BJ.summer_with_outliers_filtered", 
                        "fung.2season.BJ.winter_with_outliers_filtered", 
                        "fung.2season.F.summer_with_outliers_filtered", 
                        "fung.2season.F.winter_with_outliers_filtered")


voi.list.OT <- lapply(voi.list.OT, function(x){
  x[x>0] <- 1
  return(x)  #important to remember this line of code
})


Sys.time()
voi.inext.df.inc <- iNEXT(x = voi.list.OT,
                      q = c(0,2),
                      datatype = "incidence_raw",
                      size = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),
                      se = TRUE,
                      conf = 0.95,
                      nboot = 1000)
Sys.time()

#took about 30 seconds

#make into one list
voi.inext.all <- rbindlist(voi.inext.df.inc$iNextEst)

#add names column
voi.inext.all$group <- rep(names(voi.inext.df.inc$iNextEst), each = 60)
voi.inext.all$site <- c(rep("Heavily Oiled", 120), 
                        rep("Lightly Oiled", 120))

voi.inext.all$season <- rep(c(rep("Summer", 60), rep("Winter", 60)), 2)

#adjust factors
voi.inext.all$season <- factor(voi.inext.all$season, levels = c("Winter", "Summer"))
voi.inext.all$order <- factor(voi.inext.all$order)

#plot------

inc.plot <- voi.inext.all %>%
  ggplot(aes(x = t,
           y = qD,
           fill = season,
           shape = order)) +
  geom_ribbon(aes(ymax = qD.UCL,
                    ymin = qD.LCL),
               alpha = 0.5) +
  geom_point(inherit.aes = FALSE, data = voi.inext.all %>%
              filter(method=="extrapolated"),
             aes(x = t,
                 y = qD),
             shape = 21,
             fill = "white",
             color = "white",
             size = 3) +
  geom_point(size = 2, fill = "gray") +
  facet_wrap(~ site,  
             scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5)) +
  labs(x = "Number of Soil Cores",
     y = "Effective Number of Species",
     shape = "Hill Order") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_shape_manual(values = c(48, 50)) +
  scale_fill_manual(values = c("black", "lightgray")) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

#ggsave("images/manuscript/S2_.png", width = 8.5, height = 7, units = "in")

#plot of observed Hill orders by abundance-----

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

#plot by site----

avg_abund <- p1.summary %>%
  filter(hill_order %in% c(0,2)) %>%
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
  scale_shape_manual(values = c(48,50)) +
  guides(shape = guide_legend(label = FALSE, override.aes = list(size = 4)))

#plot all 3 plots----

plot_grid(abund, inc.plot, avg_abund, labels = c("Abund", "Incid", "Avg Abun"))

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
  select(site,season, order,qD,qD.LCL, qD.UCL,t ) %>%
  arrange(order)


voi.hill.table$min.samples <- NA
voi.hill.table$min.samples[voi.hill.table$order==0] <- A$t +1
voi.hill.table$min.samples[voi.hill.table$order==1] <- B$t + 1
voi.hill.table$min.samples[voi.hill.table$order==2] <- C$t + 1
         
names(voi.hill.table) <- c("Site", "Season", "q (Hill order)", "Effective Number of Species", "Lower CI", "Upper CI","N", "Min. Samples")

write.csv(voi.hill.table, file = "S2_fig2_fungi_stats.csv", row.names = FALSE)


#Extrapolated version

voi.hill <- subset(voi.inext.all, method=="extrapolated") %>%
  filter(t==30)

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
  select(site,season, order,qD,qD.LCL, qD.UCL,t ) %>%
  arrange(order)


voi.hill.table$min.samples <- NA
voi.hill.table$min.samples[voi.hill.table$order==0] <- A$t + 1
voi.hill.table$min.samples[voi.hill.table$order==1] <- B$t + 1
voi.hill.table$min.samples[voi.hill.table$order==2] <- C$t + 1
         
names(voi.hill.table) <- c("Site", "Season", "q (Hill order)", "Effective Number of Species", "Lower CI", "Upper CI","N", "Min. Samples")

write.csv(voi.hill.table, file = "S2_fig2_fungi_stats_extrapoltaed.csv", row.names = FALSE)
