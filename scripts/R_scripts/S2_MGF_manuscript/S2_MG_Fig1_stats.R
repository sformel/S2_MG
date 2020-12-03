#Script for S2_MG Manuscript

#Figure 1 stats:  model of PAHs between site and season
#Description:   Do our results match our experimental strategy? Are there gradients of oil within each site and a large difference between sites?

#Last updated 5 Oct 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(tidyverse)
library(brms)

#notes-----

#followed this tutorial: http://zevross.com/blog/2014/07/16/mapping-in-r-using-the-ggplot2-package/


#doesn't make a difference if I use bacteria or fungi here, they both should include all samples from the 4 dates being considered.

voi <- bac.2season_with_outliers

df <- data.frame(sample_data(voi)) %>%
  select(SampleID,
         site,
         season,
         Total.C1.C3.Chrysenes,
         Total.C1.C4.naphthalenes,
         Total.C1.C4.phenanthrenes,
         Total.C1.C3.dibenzothiophenes,
         Total.relevant.PAHs)

#df[,4:8] <- scale(df[,4:8], center = TRUE)

#Run bayesian lm model------

#Total PAHs
M.total <- brm(log(Total.relevant.PAHs) ~ 1 + site*season, 
          data = df,
          family = "skew_normal",
          chains = 4, cores = 4,
          seed = 1,
          iter = 5000, 
          control = list(max_treedepth = 15))

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- M.total
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
as.data.frame(posterior_summary(M))

#R2
as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#Total PAHs - not transformed
M.total <- brm(Total.relevant.PAHs ~ 1 + site*season, 
          data = df,
          family = "gamma",
          chains = 4, cores = 4,
          seed = 1,
          iter = 5000, 
          control = list(max_treedepth = 15))

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- M.total
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
as.data.frame(posterior_summary(M))

#R2
as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#CLR Transformations------

library("compositions")

voi <- bac.2season_with_outliers

df <- data.frame(sample_data(voi)) %>%
  select(SampleID,
         site,
         season,
         Total.C1.C3.Chrysenes,
         Total.C1.C4.naphthalenes,
         Total.C1.C4.phenanthrenes,
         Total.C1.C3.dibenzothiophenes,
         Total.relevant.PAHs)

#http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.876.3979&rep=rep1&type=pdf

#good explanations and tutorial
#https://cran.r-project.org/web/packages/compositions/vignettes/UsingCompositions.pdf

comps <- acomp(df[,c(4:7)])
plot(comps) #ternary diagram

mean(comps)
var(comps)

plot(mean(comps),pch=20,add=T,col="red") # The geometric mean

plot(acomp(comps),margin="Total.C1.C3.Chrysenes")

barplot(comps)
boxplot(comps)

totals(comps)

# 7 Multivariate Methods
# The central idea of the package ??? following the coordinate approach of [Pawlowsky-Glahn(2003)]
# and [Pawlowsky-Glahn and Mateu-Figueras(2005)] ??? is to transform the data by one
# of transforms into a classical multivariate dataset, to apply classical multivariate
# statistics and to back transform or interpreted the results afterwards in the original
# space

pc <- princomp.acomp(x = comps)

pc$Loadings # The loadings as compositional vector
pc$loadings # The loadings in clr-space
df.pca <- pc$scores

cbind(na.omit(df), df.pca) %>%
  ggplot(aes(x = Comp.1,
             y = Comp.2,
             color = season,
             fill = season,
             shape = site)) +
  geom_point(size = 4,
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
  theme_classic()

ggsave("images/manuscript/S2_MG_Fig1_supp.png")

#CLR----
clr.PAH <- clr(na.omit(df[,c(4:7)]))
clr.PAH <- as.data.frame(clr.PAH)

plot(clr.PAH)

df.clr <- cbind(na.omit(df), clr.PAH)

adonis(formula = clr.PAH ~ site*season, data = na.omit(df), permutations = 9999, method = "euclidean")

df.clr %>%
  gather(key = "PAH", value = "clr_val", Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
  ggplot(aes(x = PAH,
         y = clr_val,
         color = PAH)) +
  geom_boxplot() +
  facet_grid(rows = vars(site), cols = vars(season))

#ILR----
#great explanation of ILR
#https://stats.stackexchange.com/questions/244118/how-to-use-isometric-logratio-ilr-from-a-package-compositions

#http://www.compositionaldata.com/codawork2015/images/ProceedingsBook.pdf

#https://stats.stackexchange.com/questions/259208/how-to-perform-isometric-log-ratio-transformation
#I believe the below should use Chrysenes as the reference level for the other PAHs 
# 
# ilr.PAH <- ilr(na.omit(df[,c(4:7)]), ) 
# 
# ilr.PAH <- as.data.frame(ilr.PAH)
# 
# library(vegan)
# 
# adonis(formula = ilr.PAH ~ site*season, data = na.omit(df), permutations = 9999, method = "euclidean")
# 
# pc <- princomp(x = ilr.PAH)
# pc <- princomp(x = df[,c(4:7)])
# 
# pc$Loadings # The loadings as compositional vector
# pc$loadings # The loadings in clr-space
# df.pca <- pc$scores
# 
# cbind(na.omit(df), df.pca) %>%
#   ggplot(aes(x = Comp.1,
#              y = Comp.2,
#              color = site,
#              shape = season)) +
#   geom_point()
