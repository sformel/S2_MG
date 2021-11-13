#Make exploratory plots to visualize how OTUs correspond to a explanatory variable
#last update Nov 8 2021
#by Steve Formel

library(vegan)

 #load community matrix
data("dune")

#load environmental data as an object
data(dune.env)

#Add explanatory variable of interest to data frame

df <- dune
df$A1 <- dune.env$A1
#add column of samples names
df$sampleID <- row.names(df)

library(tidyr)

#make long version of table - Achimill and Callcusp are species names (column names)
df_long <- gather(df, species, abundance,  Achimill:Callcusp)

#plot each species
ggplot(data = df_long, aes(x = A1, y = abundance, group = species)) +
  geom_point() +
  geom_smooth(method = lm) + facet_wrap( ~ species)

