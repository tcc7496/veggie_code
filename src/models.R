

#### Libraries ----
library(nlme)
library(lme4)
library(glmm)
library(ggplot2)
library(mgcv)
library(MuMIn)
library(dplyr)


#### Working Dir ----

# set working directory
wd <- "/Users/taracunningham/projects/dissertation/sen2processing/processing/"
setwd(wd)


#### Functions ----

source('../src/texture_functions.R')
# this should contain only functions cause it runs the whole script when you enter it

# read in results of hmg per landuse
hmg_fp <- 'texture_metrics/from_evi/hmg/'
inputfile <- paste0(hmg_fp, 'hmg_per_landuse_normalised.csv')

hmg_df <- read.csv2(inputfile)

# sample dataframe
n = c(1000, 1000)

hmg_per_landuse_mini_sample <- 
  sample.df(hmg_df, n_sample = n, seed_n = 42) %>%
  convert.dtpyes.for.plot()

hmg_per_landuse_mini_sample$year <-
  hmg_per_landuse_mini_sample$year - min(hmg_per_landuse_mini_sample$year)

# check distribution of response variable
hist(hmg_df$value_norm) # poisson I think

# creates grid of plots
attach(hmg_per_landuse_mini_sample)
pairs(hmg_per_landuse_mini_sample[,c(1,3,4,6)], lower.panel = NULL, col = LandUse)

# basic linear model
model.lm <- lm(value_norm ~ LandUse + year, data = hmg_per_landuse_mini_sample)
summary(model.lm)

# include pixel id as an explanatory variable
model.lm2 <- lm(value_norm ~ LandUse + year + pixel_id, data = hmg_per_landuse_mini_sample)
# don't want this - want to control for the variation
# ref: https://ourcodingclub.github.io/tutorials/mixed-models/  Section 5

# LME model
model.lme3 <- nlme::lme(value_norm ~ year + LandUse,
                       random = ~ 1 + year | pixel_id,
                       correlation = corAR1(form = ~ year),
                       data = hmg_per_landuse_mini_sample)

pred.lme3 <- cbind(hmg_per_landuse_mini_sample, pred = predict(model.lme3))

model.lme2 <- nlme::lme(value_norm ~ year * LandUse,
                        random = ~ year | LandUse/pixel_id, 
                        data = hmg_per_landuse_mini_sample)
# ~ 1 + year | LandUse/pixel_id
# ~ 1 | pixel_id   # random intercept for each pixel
r.squaredGLMM(model.lme)

anova(model.lme, model.lme2)


# GAM model
set.seed(42) # for reproducibility

model.gam1 <- gam(value_norm ~ LandUse + s(year, k = 6), data = hmg_per_landuse_mini_sample)
model.gam2 <- gam(value_norm ~ LandUse * year, data = hmg_per_landuse_mini_sample)
model.gamm.pois <- gamm(value_norm ~ LandUse + s(year, k = 6, by = LandUse),
                   random = list(pixel_id = ~ year),
                   correlation = corAR1(form = ~ year),
                   family = c('poisson'),
                   data = hmg_per_landuse_mini_sample)

pred.gamm.pois <- cbind(hmg_per_landuse_mini_sample, pred = predict(model.gamm.pois))
#random=list(year=~1, year=~nf))

anova(model.lm, model.gam1, model.gam2, model.lme, test = 'F')
AIC(model.lm, model.gam1, model.gam2, model.lme)

ggplot(hmg_per_landuse_mini_sample, 
             aes(x = year, y = value_norm)) + 
  geom_line(aes(color = LandUse, group = pixel_id), alpha = 0.1, size = 0.1) +
  geom_line(data = pred.gamm, aes(y = pred, color = LandUse), size = 1) +
  geom_line(data = pred.gamm.pois, aes(y = pred, color = LandUse), size = 1) +
  theme_bw()
  #geom_line(data = pred.lme3, aes(y = pred), size = 1, col = 'red')
  #geom_smooth(method = 'lm', se = F, col="blue") +
  #geom_smooth(method = 'gam', formula = y ~ s(x, k = 6), se = F, col = 'blue')
  #geom_line(aes(color = pixel_id), alpha = 0.2) +
  #facet_wrap(~ LandUse) +
  #theme(legend.position = 'none')

tic('Simple ANOVA run against LandUse')
hmg_anova <- aov(value ~ LandUse + year, data = values_per_landuse_mini_sample)
toc()
summary(hmg_anova)

tic('Linear Mixed Model')
hmg.lmer <- lmer(value ~ LandUse + (1|year), data = hmg_values_per_landuse_tidy)
toc()
summary(hmg.lmer)

## GLMM

hmg.lme <- nlme::lme(value_norm ~ year + LandUse,
                     random = list(~ 1 + year | LandUse/pixel_id),
                     data = hmg_per_landuse_mini_sample)
