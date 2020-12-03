#Script for S2_MG Manuscript

#Figure 6 - dbRDA of bacteria and fungi
#Description:  Is there a linear relationship between community composition and oil abundance?

subtitle <- "S2_MG_Fig6.R"

#Last updated 16 Sep 2020 by Steve Formel

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

set.seed(1)
#load libraries-----

library(tidyverse)
library(cowplot)

#Make dbRDA----

#dbRDA of bacteria as a function of oil----

voi <- bac.2season_with_outliers

#load environmental data as an object
df.env<- data.frame(sample_data(voi))

#What is correlated?
library(ggcorrplot)

#how bad are the correlations among explanatory factors?
ggcorrplot(cor(df.env[,17:20], use = "complete.obs"), hc.order = TRUE, type = "lower",
           lab = TRUE)

#Generate PC1 for Oil----

#get rid of NA values
df.env.na <- df.env %>%
  select(Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
  na.omit(.)

df.env.na.scaled <- scale(df.env.na, center = TRUE)

PC.oil <- rda(df.env.na.scaled)

#PCA of oil to see if PC1 explains most of the variation.  If yes, then use PC1 as explanatory variable.  
round(cumsum(100*PC.oil$CA$eig / sum(PC.oil$CA$eig)),2)

#PC1 % variance explained
#80.34%

PCaxes <- scores(PC.oil)
PCaxes <- as.data.frame(PCaxes$sites)
plot(PCaxes$PC1, na.omit(df.env$Total.relevant.PAHs))

df.env.na.scaled <- as.data.frame(df.env.na.scaled)
df.env.na.scaled$PC1 <- -(PCaxes$PC1)  #flip the sign to make the relationship more intuitive

#add backinto data
sample_data(bac.2season_with_outliers)$PC1 <- NA
sample_data(fung.2season_with_outliers)$PC1 <- NA

sample_data(bac.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(bac.2season_with_outliers)), rownames(df.env.na.scaled))]

sample_data(fung.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(fung.2season_with_outliers)), rownames(df.env.na.scaled))]

#make dbRDA model-----
voi <- bac.2season_with_outliers
df.env <- data.frame(sample_data(voi))

cm <- data.frame(t(otu_table(voi)))
cm <- cm[rownames(cm) %in% rownames(df.env.na.scaled),]

exp.df <- df.env[rownames(df.env) %in% rownames(df.env.na.scaled),]

#Does scaling PAHs make a difference?  No.
#exp.df$Total.relevant.PAHs <- scale(exp.df$Total.relevant.PAHs, center = TRUE)

#relabel levels
exp.df$site <- as.factor(exp.df$site)
levels(exp.df$site) <- plyr:::revalue(levels(exp.df$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

exp.df$season <- as.factor(exp.df$season)
levels(exp.df$season) <- plyr:::revalue(levels(exp.df$season), c("SUMMER" = "Summer", "WINTER" = "Winter"))

PAH.mod <- dbrda(cm ~ site + season + Total.relevant.PAHs, 
           data = exp.df, 
           distance = "bray") 

PC1.mod <- dbrda(cm ~ site + season + PC1, 
                 data = exp.df, 
                 distance = "bray") 

set.seed(1)
anova(PAH.mod, by = "margin", permutations = 9999)
anova(PC1.mod, by = "margin", permutations = 9999)

anova(PAH.mod, by="terms", perm.max=500)
anova(PC1.mod, by="terms", perm.max=500)

#plot PAH model----

A <- PAH.mod
plot(A)
A.1 <- scores(A)
A.2 <- A.1$sites
A.3 <- cbind(A.2, exp.df)

#scores for arows
A.4 <- data.frame(scores(A, display = "bp"))
A.4 <- A.4[3,]

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = dbRDA1, yend = dbRDA2, x = 0, y = 0, shape = NULL, color = NULL, fill = NULL)
label_map <- aes(x = 1.5*dbRDA1, y = 1.5*dbRDA2, label = row.names(A.4), shape = NULL, color = NULL, fill = NULL)
arrowhead = arrow(length = unit(0.02, "npc"))


p <- ggplot(data = A.3, aes(x = dbRDA1, y = dbRDA2))

p.dbrda <- p +
  geom_point(data = A.3, 
             color = "black", 
             alpha = 5/5, 
             aes(size = log(Total.relevant.PAHs), 
                 fill = season, 
                 shape = site), 
             stroke = 0.5) +
  theme_bw() + 
  theme(plot.margin = unit(c(0,0.5,0,0), "cm"), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 8), 
        legend.position = "right") +
  geom_segment(arrow_map, size = rel(1), data = A.4, color = "black", arrow = arrowhead) +
  geom_text(label_map, size = rel(4), data = A.4, show.legend = FALSE, label = "Total PAHs") +
  labs(subtitle = "Bacteria",
       shape = "Site",
       fill = "Season") +
  theme_bw() + 
  scale_shape_manual(values=c(21, 24)) +
  scale_fill_manual(values = c("grey", "white")) +
  scale_size_continuous(name = "log(Total PAHs)", breaks = c(-0.5,1,2,4)) +
  guides(
    shape = guide_legend(override.aes = list(size= rel(4)), order = 1),
    fill = guide_legend(override.aes = list(size= rel(4), shape = 22), order = 3),
    size = guide_legend(override.aes = list(shape = 22, fill = "black"), order = 4))

p.dbrda

ggsave("images/manuscript/S2_MG_Fig6.png", width = 8, height = 6, units = "in")

#dbRDA of fungi as a function of oil----

voi <- fung.2season_with_outliers

#load environmental data as an object
df.env<- data.frame(sample_data(voi))

#What is correlated?
library(ggcorrplot)

#how bad are the correlations among explanatory factors?
ggcorrplot(cor(df.env[,17:20], use = "complete.obs"), hc.order = TRUE, type = "lower",
           lab = TRUE)

#Generate PC1 for Oil----

#get rid of NA values
df.env.na <- df.env %>%
  select(Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
  na.omit(.)

df.env.na.scaled <- scale(df.env.na, center = TRUE)

PC.oil <- rda(df.env.na.scaled)

#PCA of oil to see if PC1 explains most of the variation.  If yes, then use PC1 as explanatory variable.  
round(cumsum(100*PC.oil$CA$eig / sum(PC.oil$CA$eig)),2)

#PC1 % variance explained
#80.34%

PCaxes <- scores(PC.oil)
PCaxes <- as.data.frame(PCaxes$sites)
plot(PCaxes$PC1, na.omit(df.env$Total.relevant.PAHs))

df.env.na.scaled <- as.data.frame(df.env.na.scaled)
df.env.na.scaled$PC1 <- -(PCaxes$PC1)  #flip the sign to make the relationship more intuitive

#add backinto data
sample_data(bac.2season_with_outliers)$PC1 <- NA
sample_data(fung.2season_with_outliers)$PC1 <- NA

sample_data(bac.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(bac.2season_with_outliers)), rownames(df.env.na.scaled))]

sample_data(fung.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(fung.2season_with_outliers)), rownames(df.env.na.scaled))]

#make dbRDA model-----
voi <- fung.2season_with_outliers
df.env <- data.frame(sample_data(voi))

cm <- data.frame(t(otu_table(voi)))
cm <- cm[rownames(cm) %in% rownames(df.env.na.scaled),]

exp.df <- df.env[rownames(df.env) %in% rownames(df.env.na.scaled),]

#Does scaling PAHs make a difference?  No.
#exp.df$Total.relevant.PAHs <- scale(exp.df$Total.relevant.PAHs, center = TRUE)

PAH.mod <- dbrda(cm ~ site + season + Total.relevant.PAHs, 
                 data = exp.df, 
                 distance = "bray") 

PC1.mod <- dbrda(cm ~ site + season + PC1, 
                 data = exp.df, 
                 distance = "bray") 

set.seed(1)
anova(PAH.mod, by = "margin", permutations = 9999)
anova(PC1.mod, by = "margin", permutations = 9999)

anova(PAH.mod, by="terms", perm.max=500)
anova(PC1.mod, by="terms", perm.max=500)

#plot PAH model----

A <- PAH.mod
plot(A)
A.1 <- scores(A)
A.2 <- A.1$sites
A.3 <- cbind(A.2, exp.df)

#scores for arows
A.4 <- data.frame(scores(A, display = "bp"))
A.4 <- A.4[3,]

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = dbRDA1, yend = dbRDA2, x = 0, y = 0, shape = NULL, color = NULL, fill = NULL)
label_map <- aes(x = 1.5*dbRDA1, y = 1.5*dbRDA2, label = row.names(A.4), shape = NULL, color = NULL, fill = NULL)
arrowhead = arrow(length = unit(0.02, "npc"))


p <- ggplot(data = A.3, aes(x = dbRDA1, y = dbRDA2))

p.dbrda <- p +
  geom_point(data = A.3, 
             color = "black", 
             alpha = 5/5, 
             aes(size = log(Total.relevant.PAHs), 
                 fill = season, 
                 shape = site), 
             stroke = 0.5) +
  theme_bw() + 
  theme(plot.margin = unit(c(0,0.5,0,0), "cm"), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 8), 
        legend.position = "right") +
  geom_segment(arrow_map, size = rel(1), data = A.4, color = "black", arrow = arrowhead) +
  geom_text(label_map, size = rel(4), data = A.4, show.legend = FALSE, label = "Oiliness") +
  ggtitle("Fungi") +
  theme_bw() + 
  scale_shape_manual(values=c(21, 24)) +
  scale_fill_manual(values = c("grey", "white")) +
  scale_size_continuous(name = "Total PAHs:", breaks = c(-0.5,1,2,4)) +
  guides(
    shape = guide_legend(override.aes = list(size= rel(4)), order = 1),
    fill = guide_legend(override.aes = list(size= rel(4), shape = 22), order = 3),
    size = guide_legend(override.aes = list(shape = 22, fill = "black"), order = 4))

p.dbrda

