#
# A script to create plots
#

# Libraries ----

library(mgcv)
library(dplyr)
library(itsadug)
library(vcd)
library(tidyr)
library(ggplot2)
library(scales)

# Functions ----

# load functions
source('../src/texture_functions.R')

# Data ----

# import homogeneity data
hmg_fp <- 'texture_metrics/from_evi/hmg/norm/'
inputfile <- paste0(hmg_fp, 'hmg_norm_same_pixels_F.csv')

hmg_df <- read.csv2(inputfile)

###################################

# Histogram of data ----

hist_p <- ggplot(aes(value), data = hmg_df) +
  geom_histogram(bins = 50, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  labs(x = 'Homogeneity', y = 'Count') +
  theme_light() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

###################################

# Homogeneity raster plots ----

inputdir_abs <- 'texture_metrics/from_evi/hmg/abs/wo_outliers/'
study_area_fp <- '../../other_data/study_area_shapefile/study_area_buffer_aggr_w_conservancy.shp'

# get list of input homogeneity tifs and order them by year
inputfilelist_abs <- Sys.glob(file.path(inputdir_abs, '*wo_outliers.tif')) %>%
  str_sort(numeric = T)

# stack input files
s_abs <- raster::stack(inputfilelist_abs)

# find maximum homogeneity value for colour scale
max(s_abs)
max_hmg_abs <- 0.53 # actual max is 0.5278

# read in study area for overlaying
study_area_sf <- st_read(study_area_fp)

# define theme for maps
map_theme <- theme_light() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# create plot and extract legend in order to only have one legend in final plot for all maps
p2016_legend_abs <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[1]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  my_theme_top +
  labs(fill = 'Homogeneity') +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = margin(0, r = 10, 0, 10, "pt"))

legend_abs <- cowplot::get_legend(p2016_legend_abs)
# visualise legend to check
ggdraw(legend_norm)

# create homogeneity maps for all years
p2016 <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[1]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  map_theme +
  theme(legend.position = "none")

p2017 <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[2]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  map_theme +
  theme(legend.position = "none")

p2018 <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[3]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  map_theme +
  theme(legend.position = "none")

# bottom row
p2019 <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[4]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  map_theme +
  theme(legend.position = "none")

p2020 <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[5]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  map_theme +
  theme(legend.position = "none")

p2021 <- ggplot() +
  layer_spatial(study_area_sf, size = 0.2, col = "grey70", fill = NA) +
  layer_spatial(s_abs[[6]]) +
  scale_fill_viridis(limits = c(0, max_hmg_abs), na.value = NA) +
  map_theme +
  theme(legend.position = "none")

# put all years together in grid
hmg_maps_abs <- cowplot::plot_grid(
  plot_grid(p2016, p2017, p2018, p2019, p2020, p2021,
            labels=c("2016", "2017", "2018", "2019", "2020", "2021"),
            label_size = 10,
            label_fontface = 'plain',
            ncol = 3, nrow = 2),
  plot_grid(NULL, legend_abs, NULL, ncol = 1),
  ncol = 2,
  rel_widths=c(1, 0.14))

# visualise to check
hmg_maps_abs

# save image
save_plot(filename = 'hmg_time_series_maps.png', plot = hmg_maps_abs, base_asp = 0.8)

###################################

# Pixel value plots ----

inputdir_norm <- 'texture_metrics/from_evi/hmg/norm/'
filename_same_pix <- 'hmg_norm_same_pixels_T.csv'

# Read in homogeneity values for same pixels across all years
hmg_same_pix_df <- as.data.frame(read.csv2(paste0(inputdir_norm, filename_same_pix)))

# convert LandUse to a factor for plotting
hmg_same_pix_df$LandUse <- as.factor(hmg_same_pix_df$LandUse)
# convert year to character to divide boxes in plot between years
hmg_same_pix_df$year <- as.character(hmg_same_pix_df$year)

# take sample of dataframe to increase visibility of individual pixel lines
hmg_sample_df <- sample.df(hmg_same_pix_df, n_sample = c(1000, 1000), seed_n = 123)

# create plot
pixel_p <- ggplot(hmg_sample_df, aes(x = year, y = value_norm)) +
  geom_line(aes(color = LandUse, group = pixel_id), alpha = 0.1, size = 0.1) +
  facet_wrap(~ LandUse) +
  labs(x = 'Year', y = 'Normalised Homogeneity', colour = 'Land Use') +
  guides(colour = guide_legend(override.aes = list(size = 1, alpha = 0.9))) +
  theme_light() +
  theme(legend.position = 'bottom',
        legend.margin=margin(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12))

# visualise
pixel_p

# save plot
ggsave('pixels.png', plot = pixel_p, width = 7, height = 3.5, units = 'in', dpi = 300)

###################################

# Homogeneity Boxplots ----

inputdir_norm <- 'texture_metrics/from_evi/hmg/norm/'
filename_all_pix <- 'hmg_norm_same_pixels_F.csv'
filename_mean <- 'landuse_mean_values.csv'

# read in mean values across land uses
hmg_mean <- as.data.frame(read.csv2(paste0(inputdir_norm, filename_mean)))

# Read in homogeneity values for all pixels across all years
hmg_all_pix_df <- as.data.frame(read.csv2(paste0(inputdir_norm, filename_all_pix)))

# convert LandUse to a factor for plotting
hmg_all_pix_df$LandUse <- as.factor(hmg_all_pix_df$LandUse)
# convert year to character to divide boxes in plot between years
hmg_all_pix_df$year <- as.character(hmg_all_pix_df$year)

# calculate statistics per land use
stats_all_pix_df <- hmg_all_pix_df %>%
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
                   med_norm = median(value_norm, na.rm = T)) %>%
  ungroup()

# write out file
write.csv2(stats_all_pix_df, file = 'texture_metrics/from_evi/hmg/summary_stats_per_landuse_year.csv')

# calculate number of pixels per landuse per year when using all pixels
num_pixels <- hmg_all_pix_df %>%
  dplyr::group_by(year, LandUse) %>%
  na.omit() %>%
  summarise(count = n()) %>%
  ungroup()

# write out to csv file
write.csv2(num_pixels, file = 'texture_metrics/from_evi/hmg/norm/pixels_per_landuse.csv')

# create boxplot
p2 <- ggplot(hmg_all_pix_df, aes(x = year, y = value, fill = LandUse)) +
  geom_boxplot(alpha = 0.8,
               outlier.shape = 1,
               outlier.size = 1.5,
               outlier.stroke = 0.05,
  ) +
  ylab('Homogeneity') +
  xlab('Year') +
  labs(fill='Land use') +
  scale_fill_manual(values = c("#fbb04e", '#a2ceaa')) +
  geom_hline(aes(yintercept = mean, colour = LandUse), linetype = 'dashed', 
             data = hmg_mean, show.legend = T) +
  scale_colour_manual(name = 'Mean', 
                      values = c("#ff800e", '#638b66'), 
                      guide = guide_legend(override.aes = list(color = c("#ff800e", '#638b66')))) +
  theme_light() +
  theme(legend.position = 'bottom',
        legend.box="vertical",
        legend.margin=margin())

# visualise
p2

# save image
ggsave('hmg_box_all_pix_abs.png', plot = p2, width = 6, height = 4, units = 'in', dpi = 300)


###################################

# Linear Model Results ----

# create plot
lm_p <- ggplot(aes(year_zero + 2016, med_norm), data = stats_same_pix_df) +
  geom_point(aes(colour = LandUse)) +
  stat_smooth(method = 'lm', aes(fill = LandUse, colour = LandUse)) +
  labs(x = 'Year', y = 'Median Normalised Homogeneity', fill = 'Land Use', colour = 'Land Use') +
  scale_color_manual(values = c("#4E79A7", '#cf3e53')) +
  scale_fill_manual(values = c("#4E79A7", '#cf3e53')) +
  facet_wrap(~ LandUse) +
  theme_light() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# visualise
lm_p

# save image
ggsave('lm_plot.png', plot = lm_p, width = 7, height = 4, units = 'in', dpi = 300)


###################################

# Median Homogeneity with rainfall plot ----

# read in rainfall data
rain_fp <- 'rainfall/rainfall_sum_mar_to_june.csv'
rain_df <- as.data.frame(read.csv2(rain_fp))


# Median homogeneity with rainfall

# define scale coefficient for rainfall
# want 700 mm of rain to correspond to homogeneity = 0.25
coef2 <- 700/0.25

# change year to character
rain_df$year <- as.character(rain_df$year)

# create plot
hmg_rain_p <- ggplot(aes(x = year, y = mean, color = LandUse), data = stats_all_pix_df) +
  geom_bar(aes(x = year, y = rain/coef2, fill = LandUse), data = rain_df,
           stat = 'identity', position = position_dodge(), size = 0, width = 0.5) +
  geom_point(aes(shape = LandUse, color = LandUse), size = 2, position = position_dodge(width = 0.5)) +
  geom_line(aes(group = LandUse, color = LandUse), position = position_dodge(width = 0.5)) +
  scale_y_continuous(sec.axis = sec_axis(~.*coef2, name = "Rainfall (mm)")) +
  labs(x = 'Year', y = 'Absolute Mean Homogeneity', fill = 'Precipitation', col = 'Homogeneity', shape = 'Homogeneity') +
  scale_fill_manual(values = c("#4E79A7", '#A0CBE8')) +
  scale_color_manual(values = c("#cf3e53", '#fc7d0b')) +
  theme_light() +
  theme(legend.position = 'bottom',
        legend.box="vertical",
        legend.margin=margin(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

# visualise
hmg_rain_p

# save image
ggsave('hmg_rainfall.png', plot = hmg_rain_p, width = 7, height = 5, units = 'in', dpi = 300)


###################################

# GAMM plots ----

# load model
load('models/mod.bam.OF.rda')


# Smooths Plot

# create sequence of x values to make model predictions
x_values <- seq(0, 5, length.out = 101)
# make model predictions excluding random effects. outputs tibble
preds <- predict_gam(mod.bam.OF,
                     exclude_terms = c("s(pixel_id)", "s(pixel_id,year_zero)"),
                     values = list(year_zero = x_values))

# convert from log to original exp function. integer homogeneity
preds$exp_values <- exp(preds$fit)
preds$exp_CIupper <- exp(preds$fit + 2*preds$se.fit)
preds$exp_CIlower <- exp(preds$fit - 2*preds$se.fit)

# convert year_zero back to year
preds$year <- preds$year_zero + 2016

# create plot
smooths_p_log <- ggplot(aes(year, fit), data = preds) +
  geom_ribbon(aes(ymin = fit - 2*se.fit, ymax = fit + 2*se.fit, fill = LandUse, group = LandUse),
              alpha = 0.3) +
  geom_line(aes(colour = LandUse)) +
  labs(x = 'Year', y = 'log(Normalised Homogeneity)', fill = 'Land Use', col = 'Land Use') +
  scale_color_manual(values = c("#4E79A7", '#cf3e53')) +
  scale_fill_manual(values = c("#4E79A7", '#cf3e53')) +
  theme_light() +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# visualise
smooths_p_log

# save image
ggsave('smooths_log.png', plot = smooths_p_log, width = 7, height = 4.5, units = 'in', dpi = 300)


# Heterogeneity with rainfall plot

# create heterogeneity = inverse homogeneity column
preds$het <- 1000 - preds$exp_values
preds$CIupper_het <- 1000 - preds$exp_CIlower
preds$CIlower_het <- 1000 - preds$exp_CIupper

# convert year to number
rain_df$year <- as.numeric(rain_df$year)

# normalise rainfall to average per land use
rain_df <- rain_df %>%
  dplyr::group_by(LandUse) %>%
  dplyr::mutate(rain_norm = rain/sum(rain)) %>%
  dplyr::ungroup()

# create coefficient for second axis
# max value of heterogeneity = 826 => use 850 as maximum scale value
coef <- 850/0.34
# scale heterogeneity for plotting
preds$het_scaled <- (preds$het - 300)

# create plot
het_rain_p <- ggplot(aes(year, het_scaled), data = preds) +
  geom_bar(data = rain_df, aes(x = year, y = rain_norm*coef, fill = LandUse),
           stat = 'identity', position = position_dodge(), width = 0.5) +
  geom_line(aes(colour = LandUse), size = 1) +
  scale_fill_manual(values = c("#4E79A7", '#A0CBE8')) +
  scale_color_manual(values = c("#cf3e53", '#fc7d0b')) +
  scale_y_continuous(sec.axis = sec_axis(~.*coef, name = "Normalised Rainfall",
                                         labels = scales::scientific)) +
  scale_x_continuous(breaks = seq(2016, 2021, by = 1)) +
  labs(x = 'Year', y = 'Normalised Heterogeneity', fill = 'Precipitation', col = 'Heterogeneity') +
  theme_light() +
  theme(legend.position = 'bottom',
        legend.box="vertical",
        legend.margin=margin(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# visualise
het_rain_p

# save image
ggsave('het_rainfall.png', plot = het_rain_p, width = 7, height = 5, units = 'in', dpi = 300)


# Difference smooth plot

# get difference smooth using plot_smooths package
mod.diff <- get_smooths_difference(mod.bam13,
                                   year_zero,
                                   list(LandUse = c("Conservancy", "Livestock Rearing Area")))

# create plot
diff_smooth_p <- ggplot(aes(year_zero + 2016, difference), data = mod.diff) +
  geom_hline(aes(yintercept = 0)) +
  geom_ribbon(aes(ymin = -Inf, ymax = Inf, fill = sig_diff, group = group), alpha = 0.2) +
  scale_fill_manual(values = c( "#cf3e53", "#A0CBE8")) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.3) +
  geom_line(size = 1) +
  labs(x = 'Year',
       y = 'log(Normalised Homogeneity) Difference \n (Conserved - Livestock Rearing Area)',
       fill = 'Significant') +
  theme_light() +
  theme(legend.position = 'bottom',
        #legend.box="vertical",
        legend.margin=margin(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

#visualise
diff_smooth_p

# save plot
ggsave('diff_smooth_log.png', plot = diff_smooth_p, width = 7, height = 4.5, units = 'in', dpi = 300)

###################################

# Model check plots ----

# convert model to mgcViz object
modviz <- getViz(mod.bam.OF)

# run gam check function which outputs plots automatically
mod.ck_p <- check.gamViz(modviz13,
                        a.qq = list(method = "auto", 
                                    a.cipoly = list(fill = "lightblue")), 
                        a.respoi = list(size = 0.5), 
                        a.hist = list(bins = 50))

# visualise
mod.ck_p

###################################