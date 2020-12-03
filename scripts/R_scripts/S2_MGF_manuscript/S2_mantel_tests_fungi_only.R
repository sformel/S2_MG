#Script for S2_MG Manuscript
#Mantel test

#Description: Is there spatial autocorrelation with oil or community composition?

#Last updated 16, Oct 2020 by Steve Formel

#Load and Clean data-----

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries and notes-------

library(vegan)
library(ecodist)
library(geosphere)
library(compositions)
#Mantel tests------

#what about distance regardless of oil content?

#Subset unfiltered

fung.2season.BJ.summer.unfiltered <- fung.2season.not_filtered %>%
  subset_samples( . , (site %in% c("BJ") & season %in% c("SUMMER"))) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.BJ.winter.unfiltered <- fung.2season.not_filtered %>%
  subset_samples( . , (site %in% c("BJ") & season %in% c("WINTER"))) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.summer.unfiltered <- fung.2season.not_filtered %>%
  subset_samples( . , (site %in% c("F") & season %in% c("SUMMER"))) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}

fung.2season.F.winter.unfiltered <- fung.2season.not_filtered %>%
  subset_samples( . , (site %in% c("F") & season %in% c("WINTER"))) %>%
  {prune_taxa(taxa_sums(.) > 0 , .)}


#Get list of all bacterial and fungal phyloseq objects
varnames <- c("fung.2season.F.summer_with_outliers", 
              "fung.2season.F.winter_with_outliers", 
              "fung.2season.BJ.summer_with_outliers", 
              "fung.2season.BJ.winter_with_outliers")

# varnames <- c("bac.2season.F.summer.unfiltered", 
#               "bac.2season.F.winter.unfiltered", 
#               "bac.2season.BJ.summer.unfiltered", 
#               "bac.2season.BJ.summer.unfiltered",
#               "fung.2season.F.summer.unfiltered", 
#               "fung.2season.F.winter.unfiltered", 
#               "fung.2season.BJ.summer.unfiltered", 
#               "fung.2season.BJ.summer.unfiltered")
# 
# varnames <- c("bac.2season.F.summer", 
#               "bac.2season.F.winter", 
#               "bac.2season.BJ.summer", 
#               "bac.2season.BJ.summer",
#               "fung.2season.F.summer", 
#               "fung.2season.F.winter", 
#               "fung.2season.BJ.summer", 
#               "fung.2season.BJ.summer")

#make empty matrix
B <- matrix(ncol = 4, nrow = length(varnames))

for (i in 1:length(varnames)) {
  
  voi <- get(varnames[i])
  df <- t(data.frame(otu_table(voi)))
  df <- data.frame(df)
  
  df.env.na <- voi %>%
    sample_data() %>%
    data.frame() %>%
    select(SampleID,
           site,
           season,
           Latitude,
           Longitude) %>%
    na.omit()
  df.na <- df[df.env.na$SampleID,]
  
  dist.df <- vegdist(df.na, method = "bray")
  
  #make distance metric based on lat long
  mtx <- as.matrix(data.frame(df.env.na$Longitude, df.env.na$Latitude))
  
  # create distance matrix
  mat <- as.dist(distm(mtx, fun=distVincentyEllipsoid))
  
  man.results <- mantel(dist.df ~ mat)  #note that this is calling the eocdist version, not the vegan version
  
  man.results[7] <- ifelse(man.results[4] <= 0.05, "significant", "not signficant")
    
  B[i,] <- c("object_name" = varnames[i], man.results[c(1,4)], "result" = man.results[c(7)])
  
}

#add names and make into data frame
colnames(B) <- c("object_name", "mantel.stat(r)", "p-value", "results")
C <- as.data.frame(B) 

#view results
C[(C$results=="significant"),] 

#sort by signficance
C <- C[order(C$results=="significant", decreasing = TRUE),]

#print to copy and paste into notebook
print(C, row.names = FALSE)

write.csv(x = C,file = "images/manuscript/S2_MG_mantel_comm_comp_results.csv", row.names = FALSE)

#Mantel test on oiliness-----

#I collapse PAHs into a principle component like I've done for other tests and then made a distance matrix based on euclidean distance.

#make empty matrix
B <- matrix(ncol = 4, nrow = length(varnames))

for (i in 1:length(varnames)) {
  
  voi <- get(varnames[i])
  df.env.na <- data.frame(sample_data(voi)) %>%
  select(SampleID,
         site,
         season,
         Total.C1.C3.Chrysenes,
         Total.C1.C4.naphthalenes,
         Total.C1.C4.phenanthrenes,
         Total.C1.C3.dibenzothiophenes,
         Latitude,
         Longitude) %>%
    na.omit()
  
  comps <- acomp(df.env.na[,c(4:7)])
  pc <- princomp.acomp(x = comps)
  df.pca <- as.data.frame(pc$scores)

  dist.df <- vegdist(df.pca, method = "euclidean")
  
  #make distance metric based on lat long
  mtx <- as.matrix(data.frame(df.env.na$Longitude, df.env.na$Latitude))
  
  # create distance matrix
  mat <- as.dist(distm(mtx, fun=distVincentyEllipsoid))
  
  man.results <- mantel(dist.df ~ mat)  #note that this is calling the eocdist version, not the vegan version
  
  man.results[7] <- ifelse(man.results[4] <= 0.05, "significant", "not signficant")
  
  B[i,] <- c("object_name" = varnames[i], man.results[c(1,4)], "result" = man.results[c(7)])
  
}

#add names and make into data frame
colnames(B) <- c("object_name", "mantel.stat(r)", "p-value", "results")
C <- as.data.frame(B) 

#view results
C[(C$results=="significant"),] 

#sort by signficance
C <- C[order(C$results=="significant", decreasing = TRUE),]

print(C, row.names = FALSE)

write.csv(x = C,file = "results/S2_MG_mantel_oil_results.csv", row.names = FALSE)
