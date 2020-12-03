#Script for S2_MG Manuscript

#Figure 3 - linear fit of diversity and oil
#Description:   Is there a relationship between diversity of salt marsh soil microbial communities and oil abundance?

#Last updated 11 Nov 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(cowplot)
library(tidyverse)
library(brms)

#Generate PC1 for Oil----

voi <- fung.2season_with_outliers
df<-t(data.frame(otu_table(voi)))

#load environmental data as an object
df.env<- data.frame(sample_data(voi))

#get rid of NA values
df.env.na <- df.env %>%
  select(Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
    na.omit(.)

df.env.na.scaled <- scale(df.env.na, center = TRUE)

PC.oil <- rda(df.env.na.scaled)

PCaxes <- scores(PC.oil)
PCaxes <- as.data.frame(PCaxes$sites)
plot(PCaxes$PC1, na.omit(df.env$Total.relevant.PAHs))

df.env.na.scaled <- as.data.frame(df.env.na.scaled)
df.env.na.scaled$PC1 <- -(PCaxes$PC1)  #flip the sign to make the relationship more intuitive

#add backinto data
sample_data(fung.2season_with_outliers)$PC1 <- NA

sample_data(fung.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(fung.2season_with_outliers)), rownames(df.env.na.scaled))]


#Fungi bayesian LM of div against total relevant PAHs bacteria-----
p <- plot_richness(fung.2season_with_outliers)

#convert value to Hill numbers

p1 <- p$data
p1$hill_value <- NA

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

#total PAHs----
p1 <- p1 %>%
  filter(variable %in% c("Observed", "Shannon", "Simpson")) %>%
  droplevels()

p.list <- split.data.frame(x = p1, f = p1$variable)

#also for plotting ultimate figure
p.fung <-p.list

#Make robust linear bayesian model
#based on  https://solomonkurz.netlify.app/post/robust-linear-regression-with-the-robust-student-s-t-distribution/

#Fungi Richness----

f0 <- brm(value ~ Total.relevant.PAHs, 
          data = p.list[[1]],
          family = "gaussian",
          chains = 4, cores = 4,
          seed = 1)

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- f0
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
summ <- as.data.frame(posterior_summary(M))

#R2
r2 <- as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#get parameters for plotting
df.plot <- rbind(summ[2,], r2[1,])
df.plot$hill_order <- "0"
df.plot$param <- rownames(df.plot)

df.plot.done <- df.plot

#Fungi Shannon-----

f1 <- brm(value ~ Total.relevant.PAHs, 
          data = p.list[[2]],
          family = "skew_normal",
          chains = 4, cores = 4,
          seed = 1)

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- f1
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
summ <- as.data.frame(posterior_summary(M))

#R2
r2 <- as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#get parameters for plotting
df.plot <- rbind(summ[2,], r2[1,])
df.plot$hill_order <- "1"
df.plot$param <- rownames(df.plot)

df.plot.done <- rbind(df.plot.done, df.plot)
#Fungi Simpson-----

f2 <- brm(value ~ Total.relevant.PAHs, 
                data = p.list[[3]],
                family = "skew_normal",
                chains = 4, cores = 4,
                seed = 1,
          control = list(adapt_delta = 0.99))

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- f2
pp_check(object = M, type = "dens_overlay", nsamples = 100)
pp_check(M, type = "stat", stat = 'median', nsamples = 100)
pp_check(M, type = "stat", stat = 'mean', nsamples = 100)
pp_check(M,type = 'intervals')

plot(M)
print(M)
summ <- as.data.frame(posterior_summary(M))

#R2
r2 <- as.data.frame(bayes_R2(object = M,resp = NULL,summary = TRUE,robust = FALSE,probs = c(0.025, 0.975)))

#get parameters for plotting
df.plot <- rbind(summ[2,], r2[1,])
df.plot$hill_order <- "2"
df.plot$param <- rownames(df.plot)

df.plot.done <- rbind(df.plot.done, df.plot)
df.plot.done$microbe <- "Fungi"
fungi.plot.done <- df.plot.done

plot.params.df <- fungi.plot.done

plot.params.df

write.csv(plot.params.df, "images/manuscript/S2_MGF_final/S2_MGF_oil_div_results.csv")
