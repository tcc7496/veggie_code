# Script to run in background while still playing in exploring_textures ntbk

#### Libraries ----
library(glcm)
library(raster)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(sp)
library(tictoc)     # for timing operations
library(reticulate) # for python integration
library(viridis)    # for plotting
library(tidyr)
library(ggplot2)
library(mapproj)
library(proj4)
library(ggalt)
library(sf)
library(exactextractr)
library(glue)
library(purrr)
library(stringr)
library(lme4)
library(nlme)
library(glmm)
library(tibble)
library(patchwork)


#### Working Dir ----

# set working directory
wd <- "/Users/taracunningham/projects/dissertation/sen2processing/processing/"
setwd(wd)


#### Functions ----

source('../src/texture_functions.R')


#### Calc Texture Metrics and Zonal Stats from EVI ####

## Texture Metrics ----

# define filepaths
evi_fp <- file.path('veg_indices', 'evi') # location of evi files to process
treemask_file <- 'tree_mask/tree_mask_species_map_4_inverse_buffered.tif'
swamp_file <- 'swamp_shapefiles/swamp_clipped.shp'
study_area <- '../../other_data/study_area_shapefile/study_area_buffer_aggr_w_conservancy.shp'
evi_outdir <- file.path('texture_metrics', 'from_evi')
evi_outdir_perc <- file.path('texture_metrics', 'from_evi_perc')

# define other inputs for glcm calculation
# textures = c('hmg')
# n_grey = 256
# 
# # call function to do calculations
# calc.textures.from.evi.batch(
#   textures = textures,
#   n_grey = 256,
#   evi_fp = evi_fp,
#   aoi_fp = study_area,
#   treemask_fp = treemask_file,
#   swamp_fp = swamp_file,
#   outdir = evi_outdir
#   )

## Zonal Stats ----

#hmg_fp <- 'texture_metrics/from_evi/hmg/'
#outfile <- paste0(inputdir, 'zs_from_evi.gpkg')

# call function to calculate zonal stats
# zs_sf <- calc.zonal.stats.batch(inputdir = hmg_fp,
#                        outfile = outfile,
#                        aoi = study_area
#                        )



## Rescale outliers in hmg tifs and read out ----

inputdir = 'texture_metrics/from_evi/hmg/abs/'
probs = c(0.02, 0.98)

# function call
remove.hmg.outliers(inputdir = inputdir, probs = probs)


## Normalise values to mean per polygon in raster ----

# to observe any trend over time visually between landuses

inputdir = 'texture_metrics/from_evi/hmg/abs/'
outputdir = 'texture_metrics/from_evi/hmg/norm/'
statsfile <- 'zs_from_evi_abs.gpkg'
meanfile <- '../landuse_normalised_values.csv'

# function call
calc.normalised.hmg.raster(inputdir, outputdir, statsfile, meanfile)


## Extract all raster values per polygon ----

# Do 2 runs, one with the only the same pixels across all the years, one without
# If normalise = T, output csv has normal values too
inputdir = 'texture_metrics/from_evi/hmg/abs/'
outputdir_abs = inputdir
outputdir_norm = 'texture_metrics/from_evi/hmg/norm/'

# Run 1: same_pixels = T, normalise = T
outfile1 = paste0(outputdir_norm, 'hmg_norm_same_pixels_T.csv')

# call function
values_per_landuse_norm <- extract.values.by.polygon.batch(
  inputdir = inputdir,
  outfile = outfile1,
  aoi = study_area,
  same_pixels = T,
  normalise = T
  )

# Run 2: same_pixels = F, normalise = T
outfile2 = paste0(outputdir_norm, 'hmg_norm_same_pixels_F.csv')

# call function
values_per_landuse_norm <- extract.values.by.polygon.batch(
  inputdir = inputdir,
  outfile = outfile2,
  aoi = study_area,
  same_pixels = F,
  normalise = T
)



# write file manually
# write.csv2(values_per_landuse_norm, file = outfile)

# check number of pixels per landuse
# pixels_per_landuse <- values_per_landuse_norm %>%
#   dplyr::group_by(LandUse, year) %>%
#   dplyr::summarize(count = dplyr::n())

# Conservancy: 19996
# Livestock Rearing Area: 445817
# Loads! 2 km2 of pixels... maybe not so much


## Rainfall data ----

# inputdir <- '../original/rainfall/monthly/tifs/'
# outfile <- 'rainfall/rainfall_sum_mar_to_june.csv'
# 
# rainfall_df <- calc.rainfall(inputdir, outfile)


## NDVI ----
# calculate average NDVI across land use per year

ndvi_fp = 'veg_indices/ndvi/'
year_list

# obtain list of ndvi files
ndvi_file_list = Sys.glob(file.path(ndvi_fp, "????????_ndvi.tif")) %>%
  str_sort(numeric = T)

# open swamp shapefile and study area
swamp_shp <- st_read(swamp_file)
study_area_sf <- st_read(study_area_fp) %>%
  dplyr::select('LandUse', 'geometry')

# create stack of ndvi files
ndvi_stack <- raster::stack(ndvi_file_list) %>%
  # mask swamp to get remove main water
  mask.raster.with.shp(swamp_shp)

# remove values less than zero (corresponding to water)
ndvi_stack[ndvi_stack < 0] = NA

stats_list <- c('mean', 'stdev')
# extract mean NDVI value per landuse
# calculate mean within whole study area for each year
ndvi_zs <- calc.zonal.stats(ndvi_stack, study_area_sf, stats = stats_list)
colnames(ndvi_zs) <- c('LandUse', paste0(stats_list[1], '_', year_list),
                       paste0(stats_list[2], '_', year_list), 'geometry')



ndvi_zs <- ndvi_zs %>%
  st_drop_geometry()

# find difference in average NDVIs
diff <- as.data.frame(ndvi_zs[1,2:7] - ndvi_zs[2,2:7])
colnames(diff) <- str_replace(colnames(diff), 'mean', 'diff')

ndvi_zs_test <- cbind(ndvi_zs, diff)
ndvi_zs_test2 <- ndvi_zs_test %>%
  tidyr::pivot_longer(cols = mean_2016:diff_2021,
                      names_to = c("stat", "year"),
                      names_sep = "_",
                      names_repair = "unique",
                      values_to = "value") %>%
  tidyr::pivot_wider(names_from = stat,
                      values_from = value)

# set second round of difference to NA so it doesn't plot
ndvi_zs_test2$diff[ndvi_zs_test2$LandUse == 'Conservancy'] <- NA
  
# plot bar chart of NDVI for each land Use
ggplot(data = ndvi_zs_test2,
       aes(x = year, y = mean, fill = LandUse)) +
  geom_bar(stat = 'identity', position = position_dodge())

# plot bar chart of difference
ggplot(data = ndvi_zs_test2,
       aes(x = year, y = diff)) +
  geom_bar(stat = 'identity')

# plot line graph of 1 - mean to compare shape to smooth in models
ggplot(data = ndvi_zs_test2,
       aes(x = year, y = mean, group = LandUse, colour = LandUse)) +
  geom_line()

# significant differences in NDVI across years



## Sample dataframe for plotting ----

# number of sample pixels per landuse per year.
# total number of data points will be 6 times the sum
# to use all datapoints in conservancy, use n = 19996
# n = c(100, 100)
# 
# hmg_per_landuse_norm_sample1 <- 
#   sample.df(values_per_landuse_norm, n_sample = n, seed_n = 42) %>%
#   convert.dtpyes.for.plot()
# hmg_per_landuse_norm_sample2 <- 
#   sample.df(values_per_landuse_norm, n_sample = n, seed_n = 123) %>%
#   convert.dtpyes.for.plot()
# hmg_per_landuse_norm_sample3 <-
#   sample.df(values_per_landuse_norm, n_sample = n, seed_n = 74) %>%
#   convert.dtpyes.for.plot()
# hmg_per_landuse_norm_sample4 <-
#   sample.df(values_per_landuse_norm, n_sample = n, seed_n = 4834) %>%
#   convert.dtpyes.for.plot()

# check the number of points per landuse is indeed equal
# pixels_per_landuse_sample <- hmg_per_landuse_norm_sample %>%
#    dplyr::group_by(LandUse) %>%
#    dplyr::summarize(count = dplyr::n())
# it is indeed


## Plot data ----

# rescale year to start at 0
#hmg_per_landuse_norm_sample$year <-
#  hmg_per_landuse_norm_sample$year - min(hmg_per_landuse_norm_sample$year)

# p1 <- ggplot(hmg_per_landuse_norm_sample1, 
#        aes(x = year, y = value_norm,
#            color = LandUse, shape = LandUse, group = pixel_id)) + 
#   #geom_point(alpha = 0.5) +
#   geom_line(aes(color = LandUse), alpha = 0.05, size = 0.1) +
#   facet_wrap(~ LandUse)
#   #theme_bw()
# p1
# p2 <- ggplot(hmg_per_landuse_norm_sample2, 
#              aes(x = year, y = value_norm,
#                  color = LandUse, shape = LandUse, group = pixel_id)) + 
#   #geom_point(alpha = 0.5) +
#   geom_line(aes(color = pixel_id), alpha = 0.2) +
#   facet_wrap(~ LandUse) +
#   theme_bw() +
#   theme(legend.position = 'none')
# p2
# p3 <- ggplot(hmg_per_landuse_norm_sample3, 
#              aes(x = year, y = value_norm,
#                  color = LandUse, shape = LandUse, group = pixel_id)) + 
#   geom_point(alpha = 0.5) +
#   geom_line(aes(color = LandUse), alpha = 0.2) +
#   theme_bw()
# p4 <- ggplot(hmg_per_landuse_norm_sample4, 
#              aes(x = year, y = value_norm,
#                  color = LandUse, shape = LandUse, group = pixel_id)) + 
#   geom_point(alpha = 0.5) +
#   geom_line(aes(color = LandUse), alpha = 0.2) +
#   theme_bw()
# 
# combined <- p1 + p2 + p3 + p4 +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")
# combined
# 
# # save image
# ggsave('hmg_per_landuse_10000.png', width = 8.82, height = 3.33, units = 'in', dpi = 300)
# dev.off()

## Modelling


## Linear Mixed Effects Model - NLME package
# incorporates temporal autocorrelation

# hmg.lme <- nlme::lme(value_norm ~ year + LandUse,
#                      random = list(~ 1 + year | LandUse/pixel_id),
#                      data = values_per_landuse_norm)
# 
# 
# 
# 
# ## GLMM
# 
# # convert values column to natural numbers for compatability with package
# values_per_landuse_glmm <- hmg_values_per_landuse_tidy %>%
#   dplyr::mutate(values_nat = as.integer(value * 1000))
# 
# tic('GLMM')
# hmg.glmm <- glmm(value_nat ~ LandUse, random = list(0 ~ year),
#                  varcomps.names = c("year"),
#                  data = values_per_landuse_glmm,
#                  family.glmm = poisson.glmm,
#                  m = 100,
#                  debug = T)
# toc()
# summary(hmg.glmm)
# 
# glmm(Mate ~ 0 + Cross, random = list(~ 0 + Female,
#                                      ~ 0 + Male), varcomps.names = c("F", "M"), data = salamander,
#      family.glmm = bernoulli.glmm, m = 10^4, debug = TRUE)
# 
# 


# this is the useless bit cause you need to divide by the average per landuse across all years
# not just per year.

## Analysing Homogeneity Zonal Stats
# 
# # filepath to folder where zonal stats are stored
# zs_fp <- 'texture_metrics/from_evi/hmg/abs/'
# zs_filename <- 'zs_from_evi_abs.gpkg'
# zs_outfile_fp <- 'texture_metrics/from_evi/hmg/norm/'
# zs_outfilename <- 'zs_from_evi_norm.gpkg'
# 
# 
# ## Divide stats by the mean of each year per LandUse
# 
# # read in zonal stats
# zs_sf <- st_read(file.path(zs_fp, zs_filename))
# 
# # convert LandUse to factor for plotting
# zs_sf$LandUse <- as.factor(zs_sf$LandUse)
# # convert year to numeric
# zs_sf$year <- as.numeric(zs_sf$year)
# 
# zs_sf_mean_norm <- zs_sf %>%
#   dplyr::group_by(LandUse, year) %>%
#   dplyr::mutate(min_norm = min / mean,
#                 max_norm = max / mean,
#                 mean_norm = mean / mean,
#                 median_norm = median / mean,
#                 stdev_norm = stdev / mean,
#                 q25_norm = q25 / mean,
#                 q75_norm = q75 / mean) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-min, -max, -mean, -median, -stdev, -q25, -q75) %>%
#   dplyr::relocate(geom, .after = q75_norm)
# 
# st_write(zs_sf_mean_norm, paste0(zs_outfile_fp, zs_outfilename), append = F)
# 
# 
# # plot normalised stats
# ggplot(zs_sf_mean_norm,
#        aes(x = year,
#            ymin = min_norm,
#            lower = q25_norm, middle = median_norm, upper = q75_norm,
#            ymax = max_norm,
#            fill = LandUse)) +
#   geom_boxplot(stat = 'identity')
#   geom_point(aes(x = year, y = mean, fill = LandUse))



