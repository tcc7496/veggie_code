#
# A script to model homogeneity texture over time
#

# Libraries ----

library(mgcv)
library(dplyr)
library(itsadug)
library(vcd)
library(tidyr)

# Functions ----

# load functions
source('../src/texture_functions.R')

# GAMM ----

## Data Pre-processing

# import homogeneity data
hmg_fp <- 'texture_metrics/from_evi/hmg/norm/'
inputfile <- paste0(hmg_fp, 'hmg_norm_same_pixels_F.csv')

hmg_df <- read.csv2(inputfile)

# investigate data distribution

# plot histogram
hist(hmg_df$value)

# test fit to poisson distribution
fit_pois <- goodfit(hmg_df$value_int, type='poisson')
summary(fit_pois)
rootogram(fit_pois)

# test fit to negative binomial distribution
fit_nb <- goodfit(hmg_df$value_int, type='nbinomial')
summary(fit_nb)
rootogram(fit_nb)

# take sample of data to reduce running time of model
# ensures only pixels that appear in all years are sampled, and an equal number from each land use
n = c(1000, 1000)    # 1000 pixels per year per land use

hmg_df_sample <- 
  sample.df(hmg_df, n_sample = n, seed_n = 42)

# prepare data for putting into model

hmg_df_sample <- hmg_df_sample %>%
  # convert pixel_id and LandUse to factors
  convert.dtpyes.for.plot()
# rescale year to start at 0
dplyr::mutate(year_zero = as.numeric(hmg_df_sample$year - min(hmg_df_sample$year)),
              # rescale normalised homogeneity to integer scale
              value_int = as.integer(rescale(value_norm, to = c(1, 1000))))

# make land use an ordered factor
hmg_df_sample$LandUse_ordered <- as.ordered(hmg_df_sample$LandUse)

# change contrast to treatment coding (difference curves)
contrasts(hmg_df_sample$LandUse_ordered) <- 'contr.treatment'

# inspect contrasts
contrasts(hmg_df_sample$LandUse_ordered)


## Model

# run model
mod.bam.OF <- bam(value_int ~ LandUse_ordered                 # parametric term
                  + s(year_zero, k = 6)                   # reference smooth
                  + s(year_zero, k = 6, by = LandUse_ordered) 
                  + s(pixel_id, bs = "re")                # random intercept}
                  + s(pixel_id, year_zero, bs = "re"),    # random slope}
                  method = "fREML",
                  family = nb(),                          # negative binomial distribution
                  data = hmg_df_sample)

# look at model summary
summary(mod.bam.OF)

# perform model checks
gam.check(mod.bam.OF)

# save model
filename <- 'models/mod_bam_OF.rda'
save(mod.bam.OF, file = filename)

# visualise model smooths
plot_smooth(
  x = mod.bam.OF,
  view = 'year_zero',
  cond = list(LandUse_ordered = "Conservancy"), rm.ranef=TRUE,
  col = 'red',
  transform = exp,
  ylim = c(100, 400)
)
plot_smooth(
  x = mod.bam.OF,
  view = 'year_zero',
  cond = list(LandUse_ordered = "Livestock Rearing Area"), rm.ranef=TRUE,
  col = 'blue',
  transform = exp,
  add = TRUE
)



# Linear Model ----

## Data pre-processing

inputdir_norm <- 'texture_metrics/from_evi/hmg/norm/'
filename_same_pix <- 'hmg_norm_same_pixels_T.csv'
filename_mean <- 'landuse_normalised_values.csv'

# read in mean values across land uses
hmg_mean <- as.data.frame(read.csv2(paste0(inputdir_norm, filename_mean)))

# Read in homogeneity values for all pixels across all years
hmg_same_pix_df <- as.data.frame(read.csv2(paste0(inputdir_norm, filename_same_pix)))

# calculate statistics per land use
stats_same_pix_df <- hmg_all_pix_df %>%
  dplyr::group_by(year, LandUse) %>%
  dplyr::summarise(sd = sd(value, na.rm = T),
                   sd_norm = sd(value_norm, na.rm = T),
                   mean = mean(value, na.rm = T),
                   mean_norm = mean(value_norm, na.rm = T),
                   range = range(value, na.rm = T),
                   range_norm = range(value_norm, na.rm = T),
                   cv = cv(value, na.rm = T),
                   cv_norm = cv(value_norm, na.rm = T),
                   min = min(value, na.rm = T),
                   min_norm = min(value_norm, na.rm = T),
                   max = max(value, na.rm = T),
                   max_norm = max(value_norm, na.rm = T),
                   med = median(value, na.rm = T),
                   med_norm = median(value_norm, na.rm = T))


# convert year to number and zero it
stats_same_pix_df$year_zero <- as.numeric(stats_same_pix_df$year) - min(as.numeric(stats_same_pix_df$year))

# create dataframe for each land use
stats_same_pix_df_cons <- stats_same_pix_df[stats_same_pix_df$LandUse == 'Conservancy',]
stats_same_pix_df_lra <- stats_same_pix_df[stats_same_pix_df$LandUse == 'Livestock Rearing Area',]

# run linear models for homogeneity against time for each land use
lm.cons <- lm(med_norm ~ year_zero, data = stats_same_pix_df_cons)
lm.lra <- lm(med_norm ~ year_zero, data = stats_same_pix_df_lra)

# look at model summaries
summary(lm.cons)
summary(lm.lra)

# look at model plots
plot(lm.cons)
plot(lm.lra)

