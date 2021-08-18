# Script to run in background while still playing in exploring_textures ntbk

#### Libraries ####
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


#### Working Dir ####

# set working directory
wd <- "/Users/taracunningham/projects/dissertation/sen2processing/processing/"
setwd(wd)


#### Functions ####

# open raster
open.band <- function(band) {
  band <- raster(band)
  # calc and save min and max values of raster to object
  band <- setMinMax(band)
  return(band)
}

# mask raster with shapefile
mask.raster.with.shp <- function(r, shp) {
  # values within polygon are set to NA
  # polygon must be contained within raster
  r_masked <- mask(r, mask = shp, inverse = TRUE)
  return(r_masked)
}

# mask raster with another raster
mask.raster.with.raster <- function(r, rmask) {
  # values with data in rmask are set to NA in r
  r[!is.na(rmask)] <- NA
  return(r)
}

# prepare the data for glcm calcuation to make sure it's the same each time
rescale.outliers <- function(r) {
  # set values less than 0 to zero
  # ref: Farwell et al. pg.3
  r[r < 0] <- 0
  
  # set values > 1 to NA because meaningless - cloud or water. Negligible pixel count.
  # don't have a reference for this really.
  r[r > 1] <- NA
  return(r)
}

# function to rescale outlying percentage of data to that range 
rescale.outliers.probs <- function(x, probs = c(0.02, 0.98)) {
  ranges <- quantile(x, probs = probs, na.rm = T)
  x[x < ranges[1]] <- ranges[[1]]
  x[x > ranges[2]] <- ranges[[2]]
  print(ranges) # for testing
  return(x)
}

# rescale veg index data and change to integers
rescale.to.integers <- function(r) {
  r <- as.integer(r*100)
  return(r)
}

# rescale veg index data to 0-255 range and convert to integers
# assumes range of veg. index is [0,1]
rescale.to.eightbit <- function(r, rmin = 0, rmax = 1) {
  r <- as.integer(round(r * 255.0))
  return(r)
}

# calculate texture metrics function with defaults
calc.textures <- function(r, stats = c('homogeneity', 'variance'), n_grey = 101) {
  r_textures <- glcm(r[[1]], window = c(3,3), n_grey = n_grey,
                   shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                   statistics = stats,
                   min_x = 0, max_x = (n_grey - 1),
                   na_opt = 'any', na_val = NA)
  return(r_textures)
}

# function to difference two rasters
diff_rasters <- function(r1, r2) {
  rout <- overlay(r1, r2,
                  fun = function(r1, r2) {return(r1 - r2)} )
  return(rout)
}

calc.zonal.stats <- function(r, polygon, stats = c('min', 'max', 'count', 'mean', 'median', 'stdev', 'quantile', 'coefficient_of_variation')) {
  # calculate zonal stats on a raster within given polygons
  # r: single layer raster
  # polygon: sf or SpatialPolygon object containing one or more polygons
  # stats: character vector of statistics to compute
  # outputs the polygon sf object with an additional column for each statistic
  if (grepl('quantile', paste(stats, collapse = ''), fixed = TRUE)) {
    stats_df <- cbind(
      exact_extract(r, polygon, fun = stats, force_df = TRUE, quantiles = c(0.25, 0.75), stack_apply = T),
      LandUse = polygon$LandUse) %>%
      merge(polygon, .)
  } else {
    stats_df <- cbind(
      exact_extract(r, polygon, fun = stats, force_df = TRUE, stack_apply = T),
      LandUse = polygon$LandUse) %>%
      merge(polygon, .)
  }
  return(stats_df)
}


calc.textures.from.evi.batch <- function (
  textures, n_grey = 256, evi_fp, aoi_fp, treemask_fp, swamp_fp, outdir) {
  # Function to calculate texture metrics and zonal statistics on evi files
  
  ## Inputs ##
  # textures: character vector of texture metrics desired from the glcm calculation
  # n_grey: the number of grey levels to use in the glcm calculation
  # evi_fp: path to directory containing evi files to process
  # aoi_fp: path to the shapefile of area of interest
  # treemask_fp: path to tif file to use as treemask
  # swamp_fp: path to shapefile to mask the swamp
  # outdir: output directory for texture metrics
  
  ## Outputs ##
  # A directory is created within 'outdir' for each texture metric in
  # 'textures'. The results of the glcm calculation are written to a tif file
  # in the corresponding directory.
  
  
  # get list of evi files
  evi_inputfilelist <- Sys.glob(file.path(evi_fp, "????????_evi.tif"))
  #evi_inputfilelist <- c(file.path(evi_fp, '20210602_evi.tif'))
  
  # extract dates
  filenames_dates <- evi_inputfilelist %>%
    basename() %>%
    str_remove_all(pattern = c('evi|.tif|_')) # remove file ending
  
  # create empty list for output filenames
  outfilenames_list <- list()
  
  # loop over textures to create lists of output filepaths and names
  for (i in 1:length(textures)) {
    # define output filepath for that texture
    outpath_texture <- file.path(outdir, glue('{textures[i]}'))
    # check if directory exists and create it if not
    if (!(dir.exists(outpath_texture))) {
      dir.create(outpath_texture)
      print(glue('Directory {outpath_texture} created.'))
    }
    # create list of output filenames for that texture
    outfilenames_list[[i]] <- filenames_dates %>%
      str_remove(pattern = ".tif") %>% # remove file ending
      purrr::map(~ file.path(outpath_texture, glue('evi_{.}_{textures[i]}_n{n_grey}.tif'))) %>%
      unlist()
    names(outfilenames_list)[i] <- textures[i]
  }
  
  
  ## Process Data
  
  # open files the same for all years
  treemask <- open.band(treemask_fp) # trees where treemask = 1
  swamp_shp <- readOGR(swamp_fp)
  study_area_shp <- st_read(aoi_fp)
  
  # Loop over years
  for (i in 1:length(evi_inputfilelist)) {
    
    # Read in evi data
    evi_image <- open.band(evi_inputfilelist[i])
    
    # Prepare Data for GLCM calculation
    evi_image_prepped <- evi_image %>%
      mask.raster.with.shp(swamp_shp) %>%
      mask.raster.with.raster(treemask) %>%
      rescale.outliers() %>%
      rescale.to.eightbit
    
    # testing
    #x <- raster(ncol=5, nrow=6)
    #values(x) <- 1:ncell(x)
    
    # Calculate texture metrics
    tic('glcm calculation, 3x3 window')
    evi_glcm <- calc.textures(evi_image_prepped, n_grey = n_grey)
    #evi_glcm <- calc.textures(x, n_grey = 256)
    toc()
    
    # create list of glcm outputs
    glcm_out_list <- names(evi_glcm)
    
    # write out rasters by looping over textures
    for (j in 1:length(names(evi_glcm))) {
      writeRaster(evi_glcm[[j]],
                  filename = outfilenames_list[[j]][i],
                  format = 'GTiff', overwrite = TRUE)
    }
  }
}

calc.zonal.stats.batch <- function (inputdir, outfile, aoi, rescale_outliers = c(0.02, 0.98)) {
  ## Inputs
  # inputdir: directory with input tifs. Will try and process all tifs
  # outfile: name of output file. Will be placed in same directory as input files
  # aoi: path to shapefile of aoi
  # rescale_outliers: vector with percentage intervals to rescale outliers.
  #                   Defaults to no rescaling
  
  # get list of input files
  inputfilelist <- Sys.glob(file.path(inputdir, "*.tif")) %>%
    str_sort(numeric = T)Nor
 
   # create year last
  year_list <- c('2021', '2020', '2019', '2018', '2017', '2016') %>%
    str_sort(numeric = T)
  
  # open shapefile
  study_area <- st_read(aoi) %>%
    dplyr::select('LandUse', 'geometry')
    
  # create empty sf object for outputs
  zs_sf <- st_sf(st_sfc())
  # set crs of new sf to that of study area sf
  st_crs(zs_sf) <- st_crs(study_area)
  
  # loop over years
  for (i in 1:length(inputfilelist)) {
    # open tif
    r <- open.band(inputfilelist[i])
    # rescale outliers
    zs <- rescale.outliers.probs(r, probs = rescale_outliers) %>% # rescale outliers to (0.02, 0.98)
      calc.zonal.stats(study_area) %>%  # calculate zonal stats per polygon
      dplyr::mutate('year' = year_list[i]) # add year as a column
    # add it to zs_df
    zs_sf <- rbind(zs_sf, zs)
  }
  
  # convert year column to numeric
  zs_sf$year <- as.numeric(zs_sf$year)
  
  # write out as geopackage
  st_write(zs_sf, outfile, append = F)
  
  return(zs_sf)
}
  
extract.values.by.polygon.batch <- function (
  inputdir, aoi, same_pixels = T, normalise = F) {
  ## Inputs
  # inputdir: path to input directory of raster files
  # aoi: path to shapefile containing polygons to extract values within
  # same_pixels: If true, selects only pixels that are not NA for every year
  # normalise: If true, normalises pixels to the mean value of Landuse over all years
  
  # read in shapefile
  polys <- st_read(aoi) %>%
    dplyr::select('LandUse', 'geometry')
  
  # obtain list of input files
  inputfilelist = Sys.glob(file.path(inputdir, "*.tif")) %>%
    str_sort(numeric = T)
  
  # create stack of files
  s <- raster::stack(inputfilelist)
  
  # extract raster values per polygon
  values_per_landuse <- raster::extract(s, polys, df = T)
  
  # tidy data and convert to longform dataframe
  values_per_landuse_tidy_cp <- values_per_landuse %>%
    setNames(c('id', 'x2016', 'x2017', 'x2018', 'x2019', 'x2020', 'x2021')) %>%
    # convert landuse ID column to names of landuse
    dplyr::mutate(LandUse =
                    dplyr::case_when(id == 1 ~ polys$LandUse[1],
                                     id == 2 ~ polys$LandUse[2])
    ) %>%
    dplyr::select(-id) %>%
    # add pixel id as column
    tibble::rowid_to_column("pixel_id")
  
  
  if (normalise == T) {
    # find mean value of homogeneity per landuse across all years
    mean_per_landuse <- values_per_landuse_tidy_cp %>%
      pivot_longer(cols = starts_with("x"), names_to = "year") %>%
      dplyr::group_by(LandUse) %>%
      dplyr::summarise(mean = mean(value, na.rm = T))
    
    # Normalise values of homogeneity to calculated mean values / landuse
    values_per_landuse_tidy_cp <- values_per_landuse_tidy_cp %>%
      pivot_longer(cols =  starts_with("x"), names_to = "year") %>%
      dplyr::mutate(value_norm =
                      dplyr::case_when(LandUse == mean_per_landuse$LandUse[1]
                                       ~ value / mean_per_landuse$mean[1],
                                       LandUse == mean_per_landuse$LandUse[2]
                                       ~ value / mean_per_landuse$mean[2])
                    )
  }
  
    # remove NA observations
    values_per_landuse_tidy_cp_na <- values_per_landuse_tidy_cp %>%
      dplyr::group_by(pixel_id) %>%
      dplyr::summarise(number_of_na = sum(is.na(value))) %>%
      dplyr::right_join(values_per_landuse_tidy_cp, by = 'pixel_id') %>%
      # remove pixels that have NA observations in any year
      {if (same_pixels == T) dplyr::filter(., number_of_na == 0) else . } %>%
      # remove pixels with NA observations for all years
      {if (same_pixels == F) dplyr::filter(., number_of_na < 6) else . } %>%
      dplyr::select(-number_of_na)
  
  return(values_per_landuse_tidy_cp_na)
}

sample.df <- function (df, n_sample, seed_n = 42) {
  # df: dataframe of texture metrics with landuse and year
  # n_sample: can be a single number or vector of two numbers.
  #           If a single number, no. of pixels sampled from each landuse per year
  #           If vector, they correspond to the number of pixels from each landuse area:
  #           c(n1, n2) <- (Livestock Rearing Area, Conservancy)
  # seed_n: seed number for reproducability
  
  set.seed(seed_n)
  
  df_nested <- df %>%
    dplyr::group_by(LandUse, pixel_id) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::group_by(LandUse) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n = n_sample)
  
  df_sample <- df_nested %>%
    dplyr::mutate(samp = map2(data, n, dplyr::sample_n)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = samp) %>%
    tidyr::unnest(cols = data)
  
  return(df_sample)
}

convert.dtpyes.for.plot <- function (df) {
  # Convert pixel_id and LandUse to factors
  df$pixel_id <- as.factor(df$pixel_id)
  df$LandUse <- as.factor(df$LandUse)
  
  # Convert year to numeric
  df <- df %>%
    tidyr::separate(col = year, into = c(NA, 'year'), sep = 1)
  df$year <- as.numeric(df$year)
  
  return(df)
}

#### Calc Texture Metrics and Zonal Stats from EVI ####

## Texture Metrics

# define filepaths
evi_fp <- file.path('veg_indices', 'evi') # location of evi files to process
treemask_file <- 'tree_mask/tree_mask_species_map_4_inverse_buffered.tif'
swamp_file <- 'swamp_shapefiles/swamp_clipped.shp'
study_area <- '../../other_data/study_area_shapefile/study_area_buffer_aggr_w_conservancy.shp'
evi_outdir <- file.path('texture_metrics', 'from_evi')

# define other inputs for glcm calculation
textures = c('hmg', 'var')
n_grey = 256

# call function to do calculations
# calc.textures.from.evi.batch(
#   textures = textures,
#   n_grey = 256,
#   evi_fp = evi_fp,
#   aoi_fp = study_area,
#   treemask_fp = treemask_file,
#   swamp_fp = swamp_file,
#   outdir = evi_outdir
#   )

## Zonal Stats

hmg_fp <- 'texture_metrics/from_evi/hmg/'
outfile <- paste0(inputdir, 'zs_from_evi.gpkg')

# call function to calculate zonal stats
# zs_sf <- calc.zonal.stats.batch(inputdir = hmg_fp,
#                        outfile = outfile,
#                        aoi = study_area
#                        )


## Extract all raster values per polygon

# call function
# values_per_landuse_norm <- extract.values.by.polygon.batch(
#   inputdir = hmg_fp,
#   aoi = study_area,
#   same_pixels = T,
#   normalise = T
#   )
# 
# # check number of pixels per landuse
# pixels_per_landuse <- values_per_landuse_norm %>%
#   dplyr::group_by(LandUse) %>%
#   dplyr::summarize(count = dplyr::n())

# Conservancy: 119976
# Livestock Rearing Area: 2674902
# Loads! 12 km2 of pixels still

## Sample dataframe for plotting

# number of sample pixels per landuse per year.
# total number of data points will be 6 times the sum
n = c(100, 100)

hmg_per_landuse_norm_sample1 <- 
  sample.df(values_per_landuse_norm, n_sample = n, seed_n = 42) %>%
  convert.dtpyes.for.plot()
hmg_per_landuse_norm_sample2 <- 
  sample.df(values_per_landuse_norm, n_sample = n, seed_n = 123) %>%
  convert.dtpyes.for.plot()
hmg_per_landuse_norm_sample3 <-
  sample.df(values_per_landuse_norm, n_sample = n, seed_n = 74) %>%
  convert.dtpyes.for.plot()
hmg_per_landuse_norm_sample4 <-
  sample.df(values_per_landuse_norm, n_sample = n, seed_n = 4834) %>%
  convert.dtpyes.for.plot()

# check the number of points per landuse is indeed equal
pixels_per_landuse_sample <- hmg_per_landuse_norm_sample %>%
   dplyr::group_by(LandUse) %>%
   dplyr::summarize(count = dplyr::n())
# it is indeed


## Plot data

# rescale year to start at 0
#hmg_per_landuse_norm_sample$year <-
#  hmg_per_landuse_norm_sample$year - min(hmg_per_landuse_norm_sample$year)

p1 <- ggplot(hmg_per_landuse_norm_sample1, 
       aes(x = year, y = value_norm,
           color = LandUse, shape = LandUse, group = pixel_id)) + 
  geom_point(alpha = 0.5) +
  geom_line(aes(color = LandUse), alpha = 0.2) +
  theme_bw()

p2 <- ggplot(hmg_per_landuse_norm_sample2, 
             aes(x = year, y = value_norm,
                 color = LandUse, shape = LandUse, group = pixel_id)) + 
  geom_point(alpha = 0.5) +
  geom_line(aes(color = LandUse), alpha = 0.2) +
  theme_bw()
p3 <- ggplot(hmg_per_landuse_norm_sample3, 
             aes(x = year, y = value_norm,
                 color = LandUse, shape = LandUse, group = pixel_id)) + 
  geom_point(alpha = 0.5) +
  geom_line(aes(color = LandUse), alpha = 0.2) +
  theme_bw()
p4 <- ggplot(hmg_per_landuse_norm_sample4, 
             aes(x = year, y = value_norm,
                 color = LandUse, shape = LandUse, group = pixel_id)) + 
  geom_point(alpha = 0.5) +
  geom_line(aes(color = LandUse), alpha = 0.2) +
  theme_bw()

combined <- p1 + p2 + p3 + p4 +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined

# save image
ggsave('test2.png', width = 8.82, height = 6.67, units = 'in', dpi = 300)


## Modelling

tic('Simple ANOVA run against LandUse')
hmg_anova <- aov(value ~ LandUse, data = hmg_values_per_landuse_tidy)
toc()
summary(hmg_anova)

tic('Linear Mixed Model')
hmg.lmer <- lmer(value ~ LandUse + (1|year), data = hmg_values_per_landuse_tidy)
toc()
summary(hmg.lmer)

## Linear Mixed Effects Model - NLME package
# incorporates temporal autocorrelation






## GLMM

# convert values column to natural numbers for compatability with package
values_per_landuse_glmm <- hmg_values_per_landuse_tidy %>%
  dplyr::mutate(values_nat = as.integer(value * 1000))

tic('GLMM')
hmg.glmm <- glmm(value_nat ~ LandUse, random = list(0 ~ year),
                 varcomps.names = c("year"),
                 data = values_per_landuse_glmm,
                 family.glmm = poisson.glmm,
                 m = 100,
                 debug = T)
toc()
summary(hmg.glmm)

glmm(Mate ~ 0 + Cross, random = list(~ 0 + Female,
                                     ~ 0 + Male), varcomps.names = c("F", "M"), data = salamander,
     family.glmm = bernoulli.glmm, m = 10^4, debug = TRUE)


## Analysing Zonal Stats for EVI

# filepath to folder where zonal stats are stored
zs_fp <- 'texture_metrics/from_evi/hmg/'
zs_filename <- 'zs_from_evi.gpkg'


## Divide stats by the mean of each year per LandUse

# read in zonal stats
zs_sf <- st_read(file.path(zs_fp, zs_filename))

# convert LandUse to factor for plotting
zs_sf$LandUse <- as.factor(zs_sf$LandUse)
# convert year to numeric
zs_sf$year <- as.character(zs_sf$year)

zs_sf_mean_norm <- zs_sf %>%
  dplyr::group_by(LandUse, year) %>%
  dplyr::mutate(min_norm = min / mean,
                max_norm = max / mean,
                mean_norm = mean / mean,
                median_norm = median / mean,
                stdev_norm = stdev / mean,
                q25_norm = q25 / mean,
                q75_norm = q75 / mean) %>%
  dplyr::ungroup() %>%
  dplyr::select(-min, -max, -mean, -median, -stdev, -q25, -q75) %>%
  dplyr::relocate(geom, .after = q75_norm)


# plot normalised stats
ggplot(zs_sf_mean_norm,
       aes(x = year,
           ymin = min_norm,
           lower = q25_norm, middle = median_norm, upper = q75_norm,
           ymax = max_norm,
           fill = LandUse)) +
  geom_boxplot(stat = 'identity')
  geom_point(aes(x = year, y = mean, fill = LandUse))
  


