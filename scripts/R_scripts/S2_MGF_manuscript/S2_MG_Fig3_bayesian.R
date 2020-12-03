#Script for S2_MG Manuscript

#Figure 3 - linear fit of diversity and oil
#Description:   Is there a relationship between diversity of salt marsh soil microbial communities and oil abundance?

subtitle <- "S2_MG_Fig3_bayesian.R"

#Last updated 26 Aug 2020 by Steve Formel

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(cowplot)
library(tidyverse)
library(brms)

#Generate PC1 for Oil----

voi <- bac.2season_with_outliers
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
sample_data(bac.2season_with_outliers)$PC1 <- NA
sample_data(fung.2season_with_outliers)$PC1 <- NA

sample_data(bac.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(bac.2season_with_outliers)), rownames(df.env.na.scaled))]

sample_data(fung.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(fung.2season_with_outliers)), rownames(df.env.na.scaled))]

#Bacteria bayesian LM of div against total relevant PAHs bacteria-----
p <- plot_richness(bac.2season_with_outliers)

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

#also for plotting ultimate figure

p.list <- split.data.frame(x = p1, f = p1$variable)

p.bac <-p.list

#good explanation of checking fit
#https://towardsdatascience.com/introduction-to-bayesian-linear-regression-e66e60791ea7

#Bac Richness----

b0 <- brm(value ~ 1 + PC1, 
          data = p.list[[1]],
          family = "gaussian",
          chains = 4, cores = 4,
          seed = 1)

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- b0
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

#Bac Shannon-----

b1 <- brm(value ~ 1 + PC1, 
          data = p.list[[2]],
          family = "skew_normal",
          chains = 4, cores = 4,
          seed = 1)

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- b1
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

#Bac Simpson-----

b2 <- brm(value ~ 1 + PC1, 
          data = p.list[[3]],
          family = "skew_normal",
          chains = 4, cores = 4,
          prior(normal(200, 50), class = b),
          seed = 1)

#https://bayesed-and-confused.netlify.app/post/model-fit-checks/

M <- b2
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
df.plot$hill_order <- "2"
df.plot$param <- rownames(df.plot)

df.plot.done <- rbind(df.plot.done, df.plot)
df.plot.done$microbe <- "Bacteria"
bac.plot.done <- df.plot.done

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

f0 <- brm(value ~ 1 + PC1, 
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

f1 <- brm(value ~ 1 + PC1, 
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

f2 <- brm(value ~ 1 + PC1, 
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

plot.params.df <- rbind(fungi.plot.done, bac.plot.done)

#plot Figure 3----

#Bacteria
bac.scat <- ggplot(data = p.bac[[2]], 
                   aes(x=PC1, 
                       y=value)) +
  geom_point(aes(shape = site,
                 fill = season), 
             size = 3,
             stroke = 0.5,
             color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("lightgray", "black")) + 
  theme_bw(base_size = 14) +
  #scale_x_continuous(trans = "log10", 
   #                  labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  labs(x = "Total PAHs (\u03BCg/g)",
       y = "Effective Number \nof Species (q = 1)",
       shape = "Site",
       fill = "Season",
       subtitle = "Bacteria") +
  theme(legend.text = element_text(size = 10),
        legend.position = "right") +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5, shape=21))) 

bac.scat

#Fungi
fung.scat <- ggplot(data = p.fung[[2]], aes(x=PC1, y=value)) +
  geom_point(aes(shape = site,
                 fill = season), 
             size = 3,
             stroke = 0.5,
             color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("lightgray", "black")) + 
  theme_bw(base_size = 14) +
  #scale_x_continuous(trans = "log10", 
   #                  labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  labs(x = "Total PAHs (\u03BCg/g)",
       y = "Effective Number \nof Species (q = 1)",
       shape = "Site",
       fill = "Season",
       subtitle = "Fungi") +
  theme(legend.text = element_text(size = 10),
        legend.position = "right") +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5, shape=21))) 

fung.scat

#plot param estimates

plot.params.df$param <- plyr:::revalue(plot.params.df$param, c(b_PC1 = "\U03B2 (Slope)", R2 = "R\U00B2"))

params.left <- ggplot(data = subset(plot.params.df, plot.params.df$param=="R\U00B2"),
       aes(x = hill_order,
           shape = microbe)) +
  geom_errorbar(aes(ymin = Q2.5,
                    ymax = Q97.5),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_point(aes(y = Estimate),
             size = 3,
             fill = "gray",
             position = position_dodge(width=0.5)) +
  scale_shape_manual(values = c(5, 8)) +
  theme_bw() +
  labs(x = "Hill Order",
       subtitle = "R\U00B2") +
  theme(legend.position = "right")

params.right <- ggplot(data = subset(plot.params.df, plot.params.df$param!="R\U00B2"),
                      aes(x = hill_order,
                          shape = microbe)) +
  geom_errorbar(aes(ymin = Q2.5,
                    ymax = Q97.5),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_point(aes(y = Estimate),
             size = 3,
             fill = "gray",
             position = position_dodge(width=0.5)) +
  scale_shape_manual(values = c(5, 8)) +
  theme_bw() +
  labs(x = "Hill Order",
       subtitle = "\U03B2 (Slope)") +
  theme(legend.position = "bottom")

#put together

top <- plot_grid(bac.scat + theme(legend.position = "none"), fung.scat + theme(legend.position = "none"), labels = c("A", "B"))

top.legend <- get_legend(bac.scat )
top.done <- plot_grid(top, top.legend, nrow = 1, rel_widths = c(1, 0.2))
bottom <- plot_grid(params.right + theme(legend.position = "none"), params.left + theme(legend.position = "none"), labels = c("C", "D"), nrow = 1)
bottom.legend <- get_legend(params.left)
bottom.done <- plot_grid(bottom, bottom.legend, nrow = 1, rel_widths = c(1, 0.2))

fig3.done <- plot_grid(top.done, bottom.done, ncol = 1)
fig3.done

#ggsave("images/manuscript/S2_MG_Fig3_v3.png", height = 6, width = 8, units = "in" )