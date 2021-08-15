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



#### Path definitions ####

# paths for all calculations
treemask_file <- 'tree_mask/tree_mask_species_map_4_inverse_buffered.tif'
swamp_file <- 'swamp_shapefiles/swamp_clipped.shp'
study_area <- '../../other_data/study_area_shapefile/study_area_shp_fixed_geoms.shp'

# variables to go in filenames
products = c('evi', 'b04', 'b08')
textures = c('hmg', 'var')
n_grey = 256
levels = n_grey - 1

### EVI

# define evi input file path
evi_fp <- file.path(wd, 'veg_indices', 'evi')
# set output directory for calucations from evi
evi_outdir <- file.path(wd, 'texture_metrics', 'from_evi')

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
  outpath_texture <- file.path(evi_outdir, glue('{textures[i]}'))
  # check if directory exists and create it if not
  if (!(dir.exists(outpath_texture))) {
    dir.create(outpath_texture)
    print(glue('Directory {outpath_texture} created.'))
  }
  # create list of output filenames for that texture
  outfilenames_list[[i]] <- filenames_dates %>%
    str_remove(pattern = ".tif") %>% # remove file ending
    purrr::map(~ file.path(outpath_texture, glue('evi_{.}_{textures[i]}_n{levels}.tif'))) %>%
    unlist()
  names(outfilenames_list)[i] <- textures[i]
}


#### Process Data ####

# open files the same for all years
treemask <- open.band(treemask_file) # trees where treemask = 1
swamp_shp <- readOGR(swamp_file)
study_area_shp <- st_read(study_area)

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
  
  # Zonal Stats
  # loop over evi_glcm rasterstack
  for (k in 1:length(names(evi_glcm))) {
    zonal_stats <- calc.zonal.stats(evi_glcm[[k]], study_area_shp)
    # create output filepath
    out_fp <- str_replace(outfilenames_list[[k]][i], 'tif', 'gpkg')
    # write geopackage
    st_write(zonal_stats, out_fp)
  }
  
}

## Analysing Zonal Stats for EVI

# filepath to folder where zonal stats are stored
evi_hmg_fp <- 'texture_metrics/from_evi/hmg/'

# create list of input files and list of years. Sort to make sure they're in the same order
inputfile_list_hmg_zs <- Sys.glob(file.path(wd, evi_hmg_fp, "*.gpkg")) %>%
  str_sort(numeric = T)

year_list <- c('2021', '2020', '2019', '2018', '2017', '2016') %>%
  str_sort(numeric = T)

# open study area as sf object to obtain crs
study_area_sf <- st_read(study_area_path)
# create empty sf object
merged_sf <- st_sf(st_sfc())
# set crs of merged_sf to that of study area sf
st_crs(merged_sf) <- st_crs(study_area_sf)

# loop over years and combine into single sf dataframe
for (i in 1:length(inputfile_list_hmg_zs)) {
  zonal_stats_sf <- st_read(inputfile_list_hmg_zs[i]) %>%
    dplyr::select('LandUse', 'min', 'max', 'count', 'mean', 'median', 'stdev', 'q25', 'q75', 'coefficient_of_variation', 'geom') %>%
    dplyr::mutate('year' = year_list[i])
  merged_sf <- rbind(merged_sf, zonal_stats_sf)
}

# convert LandUse to factor for plotting
merged_sf$LandUse <- as.factor(hmg_zonal_stats_sf$LandUse)

merged_sf_minus_buffer <- merged_sf %>%
  dplyr::filter(LandUse != c('Buffer Zone'))



