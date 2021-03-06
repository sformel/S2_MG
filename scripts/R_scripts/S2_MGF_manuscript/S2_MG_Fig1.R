#Script for S2_MG Manuscript

#Figure 1 - two season map of oil distribution
#Description:   Do our results match our experimental strategy? Are there gradients of oil within each site and a large difference between sites?

#Last updated 12 Sep 2020 by Steve Formel

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

#load libraries-----

library(rgdal)
library(ggsn)
library(cowplot)
library(tidyverse)
library(plyr)

#notes-----

#followed this tutorial: http://zevross.com/blog/2014/07/16/mapping-in-r-using-the-ggplot2-package/


#make map of total PAH distribution over sites-------

#define variable of interest (VOI) - doesn't make a difference if I use bacteria or fungi here, they both should include all samples from the 4 dates being considered.

voi <- bac.2season_with_outliers

#extract oil data and omit any samples without oil data (NA)

df <- t(data.frame(otu_table(voi))) %>%
  data.frame()

df.env <- data.frame(sample_data(voi))

df.env.na <- df.env %>%
  select(SampleID,
         site,
         season,
         Total.C1.C3.Chrysenes,
         Total.C1.C4.naphthalenes,
         Total.C1.C4.phenanthrenes,
         Total.C1.C4.phenanthrenes,
         Total.relevant.PAHs,
         Longitude,
         Latitude) %>%
  na.omit()

df.na <- df[df.env.na$SampleID,]

#read in BJ shapefile
BJ.shp <- readOGR("GIS/ArcGIS/BJ_site.shp")

#check projection, should be WGS84
proj4string(BJ.shp)

# Next the shapefile has to be converted to a dataframe for use in ggplot2
BJ.shp_df <- fortify(BJ.shp)

#make our data object unique so we don't have to back transform it after this
df.env.na.maps <- df.env.na

#make it spatial points
coordinates(df.env.na.maps) <- ~ Longitude + Latitude

#we know that the coordinate system is WGS84 (google maps) so we can manually tell R what the coordinate system is
proj4string(df.env.na.maps)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# now we can use the spTransform function to project. We will project the mapdata and for coordinate reference system (CRS) we will
# assign the projection from counties

df.env.na.maps<-spTransform(df.env.na.maps, CRS(proj4string(BJ.shp)))

# double check that they match
identical(proj4string(df.env.na.maps),proj4string(BJ.shp))

# ggplot can't deal with a SpatialPointsDataFrame so we can convert back to a data.frame
df.env.na.maps_df <- data.frame(df.env.na.maps)

#subset for BJ
df.env.na.maps_df.BJ <- df.env.na.maps_df[df.env.na.maps_df$site=="BJ",]

#make Heavily Oiled (Bay Jimmy) plot-----

map <- ggplot() +
  geom_polygon(data = BJ.shp_df, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = "white", size = .2)

map.styled <- map +
  geom_point(data=df.env.na.maps_df.BJ, 
             aes(x=Longitude, 
                 y=Latitude, 
                 size= sqrt(Total.relevant.PAHs), 
                 fill = season), 
             shape = 21, 
             color = "black", 
             position=position_jitter(width = 0.00004, height = 0.00004)) +
  scale_size(range = c(rel(0), rel(10))) +
  theme_void() +
  scale_fill_manual(values = c("grey", "white")) +
  theme(legend.position = "none")

#Add north arrow and scale bar

BJ.map.no.inset <- map.styled +
  north(BJ.shp_df, 
        location = "bottomleft", 
        anchor = c(x = -89.8909, y = 29.443)) +
  scalebar(BJ.shp_df, 
           dist = 0.1, 
           dist_unit = "km", 
           transform = TRUE, 
           model = 'WGS84', 
           location = "bottomleft", 
           st.size = 2.5, 
           border.size = 0.1)

#Add coordinates if desired
#coord_cartesian(xlim = c(-89.8905, -89.88664), ylim = c(29.44503, 29.443))

#add inset of southern louisiana

#read in BJ shapefile
LA.shp <- readOGR("GIS/S6_ArcGIS/states.shp")

#now we can use the spTransform function to project. We will project the mapdata and for coordinate reference system (CRS) we will assign the projection from counties

LA.shp.tran <-spTransform(LA.shp, CRS(proj4string(BJ.shp)))

#subset to LA
LA.shp.tran <- LA.shp.tran[LA.shp.tran$STATE_NAME=="Louisiana",]

#check projection on our data
proj4string(LA.shp.tran)

# double check that they match
identical(proj4string(LA.shp.tran),proj4string(BJ.shp))

# Next the shapefile has to be converted to a dataframe for use in ggplot2
LA.shp_df <- fortify(LA.shp.tran)

#subset to southern LA
LA.shp_df.south <- LA.shp_df[LA.shp_df$long>-91.88 & LA.shp_df$lat<31.03,]

#add data point for New Orleans, Fourchon
NOLA <- data.frame("Long" = c(-90, -90.14535), "Lat" = c(30, 29.13347))


LA.inset <- ggplotGrob(ggplot() +
                         geom_polygon(data=LA.shp_df.south, 
                                      aes(x = long, y = lat, group = group), 
                                      colour="grey10",fill="white", size = 0.2) +
                         geom_point(data = df.env.na.maps_df.BJ, 
                                    aes(x = Longitude[1], y = Latitude[1]), 
                                    size = 1) +
                         geom_text(data = df.env.na.maps_df.BJ, 
                                   aes(x = Longitude[1], y = Latitude[1]), 
                                   label = "Heavily Oiled", 
                                   vjust = 0, 
                                   hjust = 1, 
                                   nudge_x = -0.05, 
                                   nudge_y = 0.05, 
                                   size = 4, 
                                   check_overlap = TRUE) +
                         geom_point(data = NOLA, 
                                    aes(x = Long, y = Lat), 
                                    size = 1) +
                         geom_text(data = NOLA, 
                                   aes(x = Long[1], y = Lat[1]), 
                                   label = c("New Orleans"), 
                                   vjust = 1, 
                                   hjust = 1, 
                                   nudge_x = -0.05, 
                                   nudge_y = -0.00, 
                                   size = 4) +
                         geom_text(data = NOLA, 
                                   aes(x = Long[2], y = Lat[2]), 
                                   label = c("Lightly Oiled"), 
                                   vjust = 1, hjust = 1, 
                                   nudge_x = -0.05, 
                                   nudge_y = -0.05, 
                                   size = 4) +
                         theme_void() +
                         theme(panel.background = element_rect(fill = "gray"), 
                               panel.border = element_rect(color = "gray", fill = NA, size = 1)))

#Heavily Oiled site map
HO.site <- BJ.map.no.inset +
  annotation_custom(grob = LA.inset, xmin = -89.8852, xmax = -89.8878, ymin = 29.4440, ymax = 29.44575)

#view
HO.site

#Make Lightly Oiled (Fourchon) site-----

#read in shapefile
F.shp <- readOGR("GIS/ArcGIS/Fourchon_site.shp")

# Next the shapefile has to be converted to a dataframe for use in ggplot2
F.shp_df <- fortify(F.shp)
F.shp_df.trans <- F.shp_df

#subset to southern LA - these coordinates were fine tuned through trial and error
F.shp_df <- F.shp_df.trans[F.shp_df.trans$long<(-90.1452) & F.shp_df.trans$lat<29.133538,]

#make this object unique so we don't have to back transform it after this
df.env.na.maps <- df.env.na

#make it spatial points
coordinates(df.env.na.maps) <- ~ Longitude + Latitude

# we know that the coordinate system is WGS84 (google maps) so we can manually tell R what the coordinate system is
proj4string(df.env.na.maps)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# now we can use the spTransform function to project. We will project the mapdata and for coordinate reference system (CRS) we will assign the projection from counties
df.env.na.maps<-spTransform(df.env.na.maps, CRS(proj4string(F.shp)))

# double check that they match
identical(proj4string(df.env.na.maps),proj4string(F.shp))

# ggplot can't deal with a SpatialPointsDataFrame so we can convert back to a data.frame
df.env.na.maps_df <- data.frame(df.env.na.maps)

#subset for Fourchon
df.env.na.maps_df.F <- df.env.na.maps_df[df.env.na.maps_df$site=="F",]

#make Fourchon map
map <- ggplot() +
  geom_polygon(data = F.shp_df, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = "white", size = .2)

#style map
map.styled <- map +
  geom_point(data=df.env.na.maps_df.F, 
             aes(x=Longitude, y=Latitude, 
                 size= sqrt(Total.relevant.PAHs), 
                 fill = season), 
             shape = 21, 
             color = "black", 
             position=position_jitter(width = 0.00002, height = 0.00002)) +
  scale_size(range = c(rel(0), rel(10))) +
  theme_void() + 
  scale_fill_manual(values = c("grey", "white")) +
  theme(legend.position = "none")

#Add north arrow and scale bar

F.map.no.inset <- map.styled +
  north(F.shp_df, 
        location = "bottomleft", 
        anchor = c(x = -90.14567, y = 29.13329)) +
  scalebar(F.shp_df, 
           dist = 0.01, 
           dist_unit = "km", 
           transform = TRUE, 
           model = 'WGS84', 
           location = "bottomleft", 
           st.size = 2.5, 
           border.size = 0.1) 


#make inset - not needed in ultimate figure, but it was done
# LA.inset.F <- ggplotGrob(ggplot() +
#                            geom_polygon(data=LA.shp_df.south, 
#                                         aes(x = long, y = lat, group = group), 
#                                         colour="grey10",fill="white", size = 0.2) +
#                            geom_point(data = NOLA, 
#                                       aes(x = Long, y = Lat), size = 1) +
#                            geom_text(data = NOLA, 
#                                      aes(x = Long, y = Lat), 
#                                      label = c("New Orleans", "Lightly Oiled"), vjust = 0, hjust = 1, nudge_x = -0.05, nudge_y = 0.00, size = 4) +
#                            theme_void())
# 
# 
# 
# #LO.site <- F.map.no.inset +
#   annotation_custom(grob = LA.inset.F, xmin = -90.145297, xmax = -90.14509, ymin = 29.1334, ymax = 29.13364)


#make strip plot of PAH classes----

df <- sample_data(bac.2season_with_outliers)

df.PAHs <- df %>%
  gather(key = "Class", 
         value = "Abundance", 
         Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes)

#make Heavily Oiled (BJ) plot
df.PAHs.plot.BJ <- ggplot(data = df.PAHs[df.PAHs$site=="BJ",], 
                          aes(x = Class, 
                              y = Abundance, 
                              color = season, 
                              fill = season, 
                              shape = site)) + 
  geom_jitter(size = 4,
              width = 0.2, 
              stroke = 1, 
              alpha = 3/5) +
  scale_color_manual(values = c(rep("black", 4)), 
                     name = "Season",
                     label = c("Heavily Oiled", "Lightly Oiled")) +
  scale_fill_manual(values = c("grey", "white"), 
                    name = "Season",
                    label = c("Heavily Oiled", "Lightly Oiled")) +
  scale_shape_manual(values = c(21,24), 
                     name = "Site",
                     label = c("Heavily Oiled", "Lightly Oiled")) +
  theme_classic() +
  theme(panel.border = element_rect(color="black", fill=NA), 
        axis.title.x = element_blank(), 
        text=element_text(size=14), 
        legend.position = "none", 
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = paste0("PAHs (", "\U03BC", "g/g)")) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_discrete(labels = c("Total.C1.C3.Chrysenes" = "C1-C3 Chrysenes", 
                              "Total.C1.C4.naphthalenes" =  "C1-C4 Naphthalenes", 
                              "Total.C1.C4.phenanthrenes" = "C1-C4 Phenanthrenes", 
                              "Total.C1.C3.dibenzothiophenes" = "C1-C3 Dibenzothiophenes") ) +
  ggtitle("Heavily Oiled")

#Get legend for final figure
legend <- get_legend(ggplot(data = df.PAHs[df.PAHs$site=="BJ",], 
                            aes(x = Class, 
                                y = Abundance, 
                                color = season, 
                                fill = season)) + 
                       geom_jitter(size = 4,
                                   width = 0.2, 
                                   stroke = 1, 
                                   alpha = 3/5, 
                                   shape = 21) +
                       scale_shape_manual(values= c(21, 21), 
                                          name = "Season",
                                          label = c("Winter", "Summer")) +
                       scale_color_manual(values = c(rep("black", 4)), 
                                          name = "Season",
                                          label = c("Winter", "Summer")) +
                       scale_fill_manual(values = c("grey", "white"), 
                                         name = "Season",
                                         label = c("Winter", "Summer")) +
                       theme(legend.position="bottom", 
                             legend.direction="horizontal", 
                             legend.justification="center"))

#make Lightly Oiled plot (Fourchon)
df.PAHs.plot.F <- ggplot(data = df.PAHs[df.PAHs$site=="F",], 
                          aes(x = Class, 
                              y = Abundance, 
                              color = season, 
                              fill = season, 
                              shape = site)) + 
  geom_jitter(size = 4,
              width = 0.2, 
              stroke = 1, 
              alpha = 3/5) +
  scale_color_manual(values = c(rep("black", 4)), 
                     name = "Season",
                     label = c("Heavily Oiled", "Lightly Oiled")) +
  scale_fill_manual(values = c("grey", "white"), 
                    name = "Season",
                    label = c("Heavily Oiled", "Lightly Oiled")) +
  scale_shape_manual(values = c(21,24), 
                     name = "Site",
                     label = c("Heavily Oiled", "Lightly Oiled")) +
  theme_classic() +
  theme(panel.border = element_rect(color="black", fill=NA), 
        axis.title.x = element_blank(), 
        text=element_text(size=14), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(y = paste0("PAHs (", "\U03BC", "g/g)")) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  scale_x_discrete(labels = c("Total.C1.C3.Chrysenes" = "C1-C3 Chrysenes", 
                              "Total.C1.C4.naphthalenes" =  "C1-C4 Naphthalenes", 
                              "Total.C1.C4.phenanthrenes" = "C1-C4 Phenanthrenes", 
                              "Total.C1.C3.dibenzothiophenes" = "C1-C3 Dibenzothiophenes") ) +
  ggtitle("Lightly Oiled")

#Build Final figure-----  

left <- plot_grid(df.PAHs.plot.BJ, 
                  df.PAHs.plot.F, 
                  legend, 
                  nrow = 3, 
                  rel_heights = c(0.6, 1, 0.1), 
                  labels = c("A", "B"), 
                  align = "V")

right <- plot_grid(HO.site, 
                   F.map.no.inset, 
                   nrow = 2, 
                   labels = c("C", "D"))

#final plot

plot_grid(left, 
          right, 
          ncol = 2, 
          rel_widths = c(1,0.75))

#save
ggsave(filename = "images/manuscript/S2_MG_Fig1_V1.png", width = 8.5, height = 6, units = "in")

