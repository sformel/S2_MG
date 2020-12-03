#Script for S2_MG Manuscript

#Figure 3 - linear fit of diversity and oil
#Description:   Is there a relationship between diversity of salt marsh soil microbial communities and oil abundance?

#Last updated 16 Oct 2020 by Steve Formel

#Load and Clean data------

source("scripts/R_scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(cowplot)
library(tidyverse)

#Generate PC1 for Oil----

voi <- fung.2season_with_outliers
df<-t(data.frame(otu_table(voi)))

#load environmental data as an object
df.env<- data.frame(sample_data(voi))

#get rid of NA values
df.env.na <- df.env %>%
  select(Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
    na.omit(.)

#df.env.na.scaled <- scale(df.env.na, center = TRUE)
df.env.na.log <- log(df.env.na)    
df.env.na.scaled <- scale(df.env.na.log, center = TRUE)

PC.oil <- rda(df.env.na.scaled)

PCaxes <- scores(PC.oil)
PCaxes <- as.data.frame(PCaxes$sites)
plot(PCaxes$PC1, na.omit(df.env$Total.relevant.PAHs))

df.env.na.scaled <- as.data.frame(df.env.na.scaled)
df.env.na.scaled$PC1 <- -(PCaxes$PC1)  #flip the sign to make the relationship more intuitive

#add backinto data
sample_data(fung.2season_with_outliers)$PC1 <- NA

sample_data(fung.2season_with_outliers)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(fung.2season_with_outliers)), rownames(df.env.na.scaled))]

#make table of linear model results-------

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

p1 <- p1 %>%
  filter(variable %in% c("Observed", "Shannon", "Simpson")) %>%
  droplevels()

p.list <- split.data.frame(x = p1, f = p1$variable)

#Total PAHs-----

p.list.lm <- list()

#Simpson's doesn't need to be transformed
p.list.lm[1] <- lapply(p.list[c(1)], lm, formula = value ~ log(Total.relevant.PAHs))
p.list.lm[c(2:3)] <- lapply(p.list[c(2:3)], lm, formula = sqrt(value) ~ log(Total.relevant.PAHs))

p.lm.summ <- lapply(p.list.lm, summary)

p.list.trans <- c("none", "square root", "square root")
hist(resid(p.list.lm[[1]]))
hist(resid(p.list.lm[[2]]))
hist(resid(p.list.lm[[3]]))

#list of p-values and R2 values

pval.list <- lapply(p.lm.summ, function(x){
  round(x$coefficients[2,4], 3)
})

r2.list <- lapply(p.lm.summ, function(x){
  round(x$r.squared, 3)
})

shapiro.list <- lapply(p.lm.summ, function(x){
  shapiro.test(resid(x))
})

shapiro.list <- lapply(shapiro.list, function(x){
  round(x$p.value, 3)
})

# #generate correlation stats
# corr.list[c(1)] <- lapply(p.list[c(1)], function(x){
#   cor.test(sqrt(x$value), log(x$Total.relevant.PAHs), method = 'spearman', exact = FALSE)
# })
# 
# corr.list[c(2:3)] <- lapply(p.list[c(2:3)], function(x){
#   cor.test(sqrt(x$value), log(x$Total.relevant.PAHs), method = 'spearman', exact = FALSE)
# })
# 
# 
# corr.list <- lapply(corr.list, function(x){
#   round(x$estimate, 3)
# })

lm.results.fungi.PAH <- data.frame("pval" = unlist(pval.list), 
                                   "R2" = unlist(r2.list),
                                   "Shapiro_pval" = unlist(shapiro.list),
                                   "hill_order" = c(0,1,2),
                                   "Microbe" = "Fungi",
                                   "Independent" = "Total_PAHs",
                                   "transformation" = p.list.trans)

#PC1----
p.list.lm[1] <- lapply(p.list[1], lm, formula = value ~ PC1)
p.list.lm[c(2:3)] <- lapply(p.list[c(2:3)], lm, formula = sqrt(value) ~ PC1)

p.lm.summ <- lapply(p.list.lm, summary)

#histograms
hist(resid(p.list.lm[[1]]), breaks = 20)
hist(resid(p.list.lm[[2]]), breaks = 20)
hist(resid(p.list.lm[[3]]), breaks = 20)

#list of p-values and R2 values

pval.list <- lapply(p.lm.summ, function(x){
  round(x$coefficients[2,4], 3)
})

r2.list <- lapply(p.lm.summ, function(x){
  round(x$r.squared, 3)
})

shapiro.list <- lapply(p.lm.summ, function(x){
  shapiro.test(resid(x))
})

shapiro.list <- lapply(shapiro.list, function(x){
  round(x$p.value, 3)
})

# #generate correlation stats
# corr.list[c(1)] <- lapply(p.list[c(1)], function(x){
#   cor.test(x$value^2, x$PC1, method = 'spearman', exact = FALSE)
# })
# 
# corr.list[c(2:3)] <- lapply(p.list[c(2:3)], function(x){
#   cor.test(sqrt(x$value), x$PC1, method = 'spearman', exact = FALSE)
# })
# 
# corr.list <- lapply(corr.list, function(x){
#   round(x$estimate, 3)
# })

lm.results.fungi.PC1 <- data.frame("pval" = unlist(pval.list), 
                                   "R2" = unlist(r2.list),
                                   "Shapiro_pval" = unlist(shapiro.list),
                                   "hill_order" = c(0,1,2),
                                   "Microbe" = "Fungi",
                                   "Independent" = "PC1",
                                   "transformation" = p.list.trans)


lm.results.fungi <- rbind(lm.results.fungi.PAH, lm.results.fungi.PC1)

lm.results.fungi

#Combine-----

lm.plot.df <- rbind(lm.results.bacteria, lm.results.fungi)

#What are the patterns?
  lm.plot.df[order(lm.plot.df$hill_order),]

#So it looks like fungal significance increases with order, but only for PC1
#Bacteria does the opposite.  But this could be due to the fact that fungi are less diverse and each unit change would be proportionally larger.  To check, I'll conpare the slopes of q = 0 to q =1

#plot richness vs Shannon for fungi and bac-----

p.bac <- plot_richness(bac.2season_with_outliers)
p.fung <- plot_richness(fung.2season_with_outliers)

p.bac <- p.bac$data
p.fung <- p.fung$data

plot(exp(p.bac$value[which(p.bac$variable == "Shannon")]) ~ p.bac$value[which(p.bac$variable == "Observed")])
summary(lm(exp(p.bac$value[which(p.bac$variable == "Shannon")]) ~ p.bac$value[which(p.bac$variable == "Observed")]))

plot(exp(p.fung$value[which(p.fung$variable == "Shannon")]) ~ p.fung$value[which(p.fung$variable == "Observed")])
summary(lm(exp(p.fung$value[which(p.fung$variable == "Shannon")]) ~ p.fung$value[which(p.fung$variable == "Observed")]))

#interpreting transformations
#https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faqhow-do-i-interpret-a-regression-model-when-some-variables-are-log-transformed/

#plot of stats----
lm.plot.df$hill_order <- as.factor(lm.plot.df$hill_order)

lm.plot.df.gathered <- lm.plot.df %>%
  gather(key = "Statistic", value = "stat.value", pval:Spear_corr_rho)

ggplot(data = lm.plot.df.gathered,
       aes(x = hill_order,
           y = stat.value)) +
  geom_point(aes(color = Independent,
                 shape = Microbe),
              size = 3) +
  facet_wrap(~ Statistic, scales = "free")

#They only become linear once I log transform Total PAHs
summary(lm(PC1 ~ Total.relevant.PAHs, data = p1)) #R2 = 0.5

#plot LM of shannon's div against total relevant PAHs bacteria-----
p <- plot_richness(bac.2season_with_outliers)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))



#get R2 values for below plot
bac.lm <- lm(p1$exp_Shannon^2 ~ p1$Total.relevant.PAHs) 
lm.sum <- summary(bac.lm)
hist(resid(bac.lm))
qqnorm(resid(bac.lm), pch = 1, frame = FALSE)
qqline(resid(bac.lm), col = "steelblue", lwd = 2)
shapiro.test(resid(bac.lm))

#plot
p3.1 <- ggplot(data = p1, aes(x=Total.relevant.PAHs, y=exp_Shannon))

bac.plot <- p3.1 +
  geom_point(aes(shape = site,
                 color = season), 
             size = 3,
             stroke = 1) +
  scale_shape_manual(values = c(21, 8)) +
  scale_color_manual(values = c("darkgray", "black")) + 
  theme_bw(base_size = 14) +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  labs(x = "Total PAHs (\U03BCg/g)",
       y = "Effective Number of Species (q = 1)",
       shape = "Site",
       color = "Season") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("text", 
           hjust = 1, 
           x = 1.5, 
           y = 225, 
           label = paste("R� = ", 
                         round(lm.sum$r.squared, 3), 
                         "\np = ", 
                         round(lm.sum$coefficients[2,4], 3))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = rel(1.5))

bac.plot

#plot LM of shannon's div against total relevant PAHs fungi-----
p <- plot_richness(fung.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

#get R2 values for below plot
fung.lm <- lm(p1$exp_Shannon ~ p1$PC1)
lm.sum <- summary(fung.lm)
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#plot
p3.1 <- ggplot(data = p1, aes(x=PC1, y=exp_Shannon))

fung.plot <- p3.1 +
  geom_point(aes(shape = site,
                 color = season), 
             size = 3,
             stroke = 1) +
  scale_shape_manual(values = c(21, 8)) +
  scale_color_manual(values = c("darkgray", "black")) + 
  theme_bw(base_size = 14) +
  #scale_x_continuous(trans = "log10") +
  theme_bw() +
  labs(x = "Relative Oil Abundance (PC1)",
       y = "Effective Number of Species (q = 1)",
       shape = "Site",
       color = "Season") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("text", 
           hjust = 1, 
           x = 1.5, 
           y = 10, 
           label = paste("R� = ", 
                         round(lm.sum$r.squared, 3), 
                         "\np = ", 
                         round(lm.sum$coefficients[2,4], 3))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = rel(1.5))

fung.plot

#plot them together-----

plots <- plot_grid(bac.plot +theme(legend.position = "none"), fung.plot + theme(legend.position = "none") + labs(y = element_blank()), nrow = 1, labels = c("A", "B"))

plot_grid(plots, get_legend(bac.plot), ncol = 1, rel_heights = c(1,0.2))

ggsave("images/manuscript/S2_MG_Fig3_v1.png", width = 9, height = 5, units = "in")

#Total PAHs-----

#It occurred to me that Total PAHs is easier to read than PC1 and also gives the reader a reference point that shows oil was there in reasonable amounts, we're not hiding behind the term relative.

#plot LM of shannon's div against total relevant PAHs bacteria-----
p <- plot_richness(bac.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))



#get R2 values for below plot
bac.lm <- lm(p1$exp_Shannon ~ (p1$Total.relevant.PAHs)) 
lm.sum <- summary(bac.lm)
hist(resid(bac.lm))
qqnorm(resid(bac.lm), pch = 1, frame = FALSE)
qqline(resid(bac.lm), col = "steelblue", lwd = 2)
shapiro.test(resid(bac.lm))

#squaring the Shannon's exp makes it normal, but the results are essentially the same, so I feel comfortable delivering the results this way.

#plot
p3.1 <- ggplot(data = p1, aes(x=Total.relevant.PAHs, y=exp_Shannon))

bac.plot <- p3.1 +
  geom_point(aes(shape = site,
                 color = season), 
             size = 3,
             stroke = 1) +
  scale_shape_manual(values = c(21, 8)) +
  scale_color_manual(values = c("darkgray", "black")) + 
  theme_bw(base_size = 14) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  labs(x = "Total PAHs (\U03BCg/g)",
       y = "Effective Number of Species (q = 1)",
       shape = "Site",
       color = "Season",
       subtitle = "Bacteria") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("text", 
           hjust = 1, 
           x = 10, 
           y = 225, 
           label = paste("R� = ", 
                         round(lm.sum$r.squared, 3), 
                         "\np = ", 
                         round(lm.sum$coefficients[2,4], 3))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = rel(1.5))

bac.plot

#plot LM of shannon's div against total relevant PAHs fungi-----
p <- plot_richness(fung.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

#get R2 values for below plot
fung.lm <- lm(p1$exp_Shannon ~ p1$Total.relevant.PAHs)
lm.sum <- summary(lm(p1$exp_Shannon ~ p1$Total.relevant.PAHs))
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#plot
p3.1 <- ggplot(data = p1, aes(x=Total.relevant.PAHs, y=exp_Shannon))

fung.plot <- p3.1 +
  geom_point(aes(shape = site,
                 color = season), 
             size = 3,
             stroke = 1) +
  scale_shape_manual(values = c(21, 8)) +
  scale_color_manual(values = c("darkgray", "black")) + 
  theme_bw(base_size = 14) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  labs(x = "Total PAHs (\U03BCg/g)",
       y = "Effective Number of Species (q = 1)",
       shape = "Site",
       color = "Season",
       subtitle = "Fungi") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("text", 
           hjust = 1, 
           x = 10, 
           y = 10, 
           label = paste("R� = ", 
                         round(lm.sum$r.squared, 3), 
                         "\np = ", 
                         round(lm.sum$coefficients[2,4], 3))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = rel(1.5))

fung.plot

#plot them together-----

plots <- plot_grid(bac.plot +theme(legend.position = "none"), fung.plot + theme(legend.position = "none") + labs(y = element_blank()), nrow = 1, labels = c("A", "B"))

plot_grid(plots, get_legend(bac.plot), ncol = 1, rel_heights = c(1,0.2))

ggsave("images/manuscript/S2_MG_Fig3_v2.png", width = 9, height = 5, units = "in")


#Answer questions about LM------

#So none of these models meet the normality assumptions unless you square the diversity
#So are any of the other orders of diversity strictly "significant"? 

#plot LM of richness against oil fungi-----
p <- plot_richness(fung.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Observed")

#LM model
fung.lm <- lm(p1$value ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#richness is not significant but is normal for Total PAHs

#LM model
fung.lm <- lm(p1$value ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#same for PC1

#plot LM of q = 1 against oil fungi-----
p <- plot_richness(fung.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#LM model
fung.lm <- lm(p1$exp_Shannon ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#q = 1 is not significant and is not normal for Total PAHs

#LM model
fung.lm <- lm(p1$exp_Shannon ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#q = 1 is significant for PC1 but not normal

#Square root exp shannons

#LM model
fung.lm <- lm(sqrt(p1$exp_Shannon) ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#q = 1 is significant for PC1 and normal if sqrt transformed

#LM model
fung.lm <- lm(sqrt(p1$exp_Shannon) ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#q = 1 is not significant for Total PAHs and normal if sqrt transformed

#plot LM of q = 2 against oil fungi-----
p <- plot_richness(fung.2season)

#subset to Simpson
p1 <- subset(p$data, p$data$variable=="Simpson")

#revalue to Hill q = 2
p1$simpson2 <- 1/(1-p1$value)

#LM model
fung.lm <- lm(p1$simpson2 ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#significant but not normal

#LM model
fung.lm <- lm(p1$simpson2 ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#not significant or normal

#Transformable?

#For neither PAHs or PC1 using Simpson's diveristy or Hill q =2 could I transform the data to be normal.

#Bacteria------

#plot LM of richness against oil fungi-----
p <- plot_richness(bac.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Observed")

#LM model
fung.lm <- lm(p1$value ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#richness is not significant and not normal for Total PAHs
#square will transform it but does not make significant

#LM model
fung.lm <- lm(p1$value ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#same for PC1

#plot LM of q = 1 against oil fungi-----
p <- plot_richness(bac.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#LM model
fung.lm <- lm(p1$exp_Shannon ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#LM model
fung.lm <- lm(p1$exp_Shannon ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#plot LM of q = 2 against oil fungi-----
p <- plot_richness(bac.2season)

#subset to Simpson
p1 <- subset(p$data, p$data$variable=="Simpson")

#revalue to Hill q = 2
p1$simpson2 <- 1/(1-p1$value)

#LM model
fung.lm <- lm(p1$simpson2^2 ~ p1$PC1)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#significant but not normal

#LM model
fung.lm <- lm(p1$simpson2^2 ~ p1$Total.relevant.PAHs)
lm.sum <- summary(fung.lm)
lm.sum
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#not significant or normal

#Transformable?

#For neither PAHs or PC1 using Simpson's diveristy or Hill q =2 could I transform the data to be normal.

#Plot points and p-values ~ Hill order------

#Load and Clean data------

source("scripts/S2_load_packages_and_clean_data.R")

#load libraries-----
library(cowplot)
library(tidyverse)

#Generate PC1 for Oil----

voi <- bac.2season
df<-t(data.frame(otu_table(voi)))

#load environmental data as an object
df.env<- data.frame(sample_data(voi))

#get rid of NA values
df.env.na <- df.env %>%
  select(Total.C1.C3.Chrysenes:Total.C1.C3.dibenzothiophenes) %>%
  na.omit(.)

#df.env.na.scaled <- scale(df.env.na, center = TRUE)
df.env.na.log <- log(df.env.na)    
df.env.na.scaled <- scale(df.env.na.log, center = TRUE)

PC.oil <- rda(df.env.na.scaled)

PCaxes <- scores(PC.oil)
PCaxes <- as.data.frame(PCaxes$sites)
plot(PCaxes$PC1, na.omit(df.env$Total.relevant.PAHs))

df.env.na.scaled <- as.data.frame(df.env.na.scaled)
df.env.na.scaled$PC1 <- -(PCaxes$PC1)  #flip the sign to make the relationship more intuitive

#add backinto data
sample_data(bac.2season)$PC1 <- NA
sample_data(fung.2season)$PC1 <- NA

sample_data(bac.2season)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(bac.2season)), rownames(df.env.na.scaled))]

sample_data(fung.2season)$PC1 <- df.env.na.scaled$PC1[match(rownames(sample_data(fung.2season)), rownames(df.env.na.scaled))]

#Bacteria-----
p <- plot_richness(bac.2season)

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

#make lm for each combination----
p1 <- p1 %>%
  filter(variable %in% c("Observed", "Shannon", "Simpson")) %>%
  droplevels()

p.list <- split.data.frame(x = p1, f = p1$variable)
p.list.lm <- lapply(p.list, lm, formula = value ~ Total.relevant.PAHs)
p.lm.summ <- lapply(p.list.lm, summary)

#list of p-values and R2 values

pval.list <- lapply(p.lm.summ, function(x){
  round(x$coefficients[2,4], 3)
  })

r2.list <- lapply(p.lm.summ, function(x){
  round(x$r.squared, 3)
})

shapiro.list <- lapply(p.lm.summ, function(x){
  shapiro.test(resid(x))
})

shapiro.list <- lapply(shapiro.list, function(x){
  round(x$p.value, 3)
})

lm.results.bacteria.PAH <- data.frame("pval" = unlist(pval.list), 
                         "R2" = unlist(r2.list),
                         "Shapiro_pval" = unlist(shapiro.list),
                         "hill_order" = c(0,1,2),
                         "Microbe" = "Bacteria",
                         "Independent" = "Total_PAHs")

p.list.lm <- lapply(p.list, lm, formula = value ~ PC1)

p.lm.summ <- lapply(p.list.lm, summary)

#list of p-values and R2 values

pval.list <- lapply(p.lm.summ, function(x){
  round(x$coefficients[2,4], 3)
})

r2.list <- lapply(p.lm.summ, function(x){
  round(x$r.squared, 3)
})

shapiro.list <- lapply(p.lm.summ, function(x){
  shapiro.test(resid(x))
})

shapiro.list <- lapply(shapiro.list, function(x){
  round(x$p.value, 3)
})

lm.results.bacteria.PC1 <- data.frame("pval" = unlist(pval.list), 
                                  "R2" = unlist(r2.list),
                                  "Shapiro_pval" = unlist(shapiro.list),
                                  "hill_order" = c(0,1,2),
                                  "Microbe" = "Bacteria",
                                  "Independent" = "PC1")


lm.results.bacteria <- rbind(lm.results.bacteria.PAH, lm.results.bacteria.PC1)

#Fungi----

p <- plot_richness(fung.2season)

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

#make lm for each combination----
p1 <- p1 %>%
  filter(variable %in% c("Observed", "Shannon", "Simpson")) %>%
  droplevels()

p.list <- split.data.frame(x = p1, f = p1$variable)
p.list.lm <- lapply(p.list, lm, formula = value ~ Total.relevant.PAHs)
p.lm.summ <- lapply(p.list.lm, summary)

#list of p-values and R2 values

pval.list <- lapply(p.lm.summ, function(x){
  round(x$coefficients[2,4], 3)
})

r2.list <- lapply(p.lm.summ, function(x){
  round(x$r.squared, 3)
})

shapiro.list <- lapply(p.lm.summ, function(x){
  shapiro.test(resid(x))
})

shapiro.list <- lapply(shapiro.list, function(x){
  round(x$p.value, 3)
})

lm.results.fungi.PAH <- data.frame("pval" = unlist(pval.list), 
                                      "R2" = unlist(r2.list),
                                      "Shapiro_pval" = unlist(shapiro.list),
                                      "hill_order" = c(0,1,2),
                                      "Microbe" = "Fungi",
                                      "Independent" = "Total_PAHs")

p.list.lm <- lapply(p.list, lm, formula = value ~ PC1)

p.lm.summ <- lapply(p.list.lm, summary)

#histograms
hist(resid(p.list.lm[[1]]), breaks = 20)
hist(resid(p.list.lm[[2]]), breaks = 20)
hist(resid(p.list.lm[[3]]), breaks = 20)

#list of p-values and R2 values

pval.list <- lapply(p.lm.summ, function(x){
  round(x$coefficients[2,4], 3)
})

r2.list <- lapply(p.lm.summ, function(x){
  round(x$r.squared, 3)
})

shapiro.list <- lapply(p.lm.summ, function(x){
  shapiro.test(resid(x))
})

shapiro.list <- lapply(shapiro.list, function(x){
  round(x$p.value, 3)
})

#What does a KS test for normality look like?
KS.list <- lapply(p.lm.summ, function(x){
  ks.test(resid(x), "pnorm", mean = 1, sd = 2)
})

lm.results.fungi.PC1 <- data.frame("pval" = unlist(pval.list), 
                                      "R2" = unlist(r2.list),
                                      "Shapiro_pval" = unlist(shapiro.list),
                                      "hill_order" = c(0,1,2),
                                      "Microbe" = "Fungi",
                                      "Independent" = "PC1")


lm.results.fungi <- rbind(lm.results.fungi.PAH, lm.results.fungi.PC1)

#Combine

lm.plot.df <- rbind(lm.results.bacteria, lm.results.fungi)

lm.plot.df[order(lm.plot.df$Shapiro_pval),]

#plot
p3.1 <- ggplot(data = p1, aes(x=PC1, y=exp_Shannon))

bac.plot <- p3.1 +
  geom_point(aes(shape = site,
                 color = season), 
             size = 3,
             stroke = 1) +
  scale_shape_manual(values = c(21, 8)) +
  scale_color_manual(values = c("darkgray", "black")) + 
  theme_bw(base_size = 14) +
  #scale_x_continuous(trans = "log10") +
  theme_bw() +
  labs(x = "Relative Oil Abundance (PC1)",
       y = "Effective Number of Species (q = 1)",
       shape = "Site",
       color = "Season") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("text", 
           hjust = 1, 
           x = 1.5, 
           y = 225, 
           label = paste("R� = ", 
                         round(lm.sum$r.squared, 3), 
                         "\np = ", 
                         round(lm.sum$coefficients[2,4], 3))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = rel(1.5))

bac.plot

#plot LM of shannon's div against total relevant PAHs fungi-----
p <- plot_richness(fung.2season)

#subset to Shannon
p1 <- subset(p$data, p$data$variable=="Shannon")

#revalue Shannon to Hill q = 1
p1$exp_Shannon <- exp(p1$value)

#rename levels for site
p1$site <- as.factor(p1$site)
levels(p1$site) <- plyr:::revalue(levels(p1$site), c("BJ" = "Heavily Oiled", "F" = "Lightly Oiled"))

#get R2 values for below plot
fung.lm <- lm(p1$exp_Shannon ~ p1$PC1)
lm.sum <- summary(fung.lm)
hist(resid(fung.lm))
shapiro.test(resid(fung.lm))

#plot
p3.1 <- ggplot(data = p1, aes(x=PC1, y=exp_Shannon))

fung.plot <- p3.1 +
  geom_point(aes(shape = site,
                 color = season), 
             size = 3,
             stroke = 1) +
  scale_shape_manual(values = c(21, 8)) +
  scale_color_manual(values = c("darkgray", "black")) + 
  theme_bw(base_size = 14) +
  #scale_x_continuous(trans = "log10") +
  theme_bw() +
  labs(x = "Relative Oil Abundance (PC1)",
       y = "Effective Number of Species (q = 1)",
       shape = "Site",
       color = "Season") +
  theme(legend.text = element_text(size = 10),
        legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("text", 
           hjust = 1, 
           x = 1.5, 
           y = 10, 
           label = paste("R� = ", 
                         round(lm.sum$r.squared, 3), 
                         "\np = ", 
                         round(lm.sum$coefficients[2,4], 3))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = rel(1.5))

fung.plot

#plot them together-----

plots <- plot_grid(bac.plot +theme(legend.position = "none"), fung.plot + theme(legend.position = "none") + labs(y = element_blank()), nrow = 1, labels = c("A", "B"))

plot_grid(plots, get_legend(bac.plot), ncol = 1, rel_heights = c(1,0.2))

ggsave("images/manuscript/S2_MG_Fig3_v1.png", width = 9, height = 5, units = "in")