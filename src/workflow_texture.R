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


#### Working Dir ####

# set working directory
wd <- "/Users/taracunningham/projects/dissertation/sen2processing/processing/"
setwd(wd)


#### Functions ####

# open raster
open.band <- function(band) {
  band <- raster(paste0(wd, band))
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
  r <- as.integer(round(f * 255.0))
}

# calculate texture metrics function with defaults
calc.textures <- function(r, stats = c('homogeneity', 'variance'), n_grey = 101) {
  r_textures <- glcm(r[[1]], window = c(3,3), n_grey = n_grey,
                   shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                   statistics = stats,
                   min_x = 0, max_x = 100,
                   na_opt = 'any', na_val = NA)
  return(r_textures)
}

# function to difference two rasters
diff_rasters <- function(r1, r2) {
  rout <- overlay(r1, r2,
                  fun = function(r1, r2) {return(r1 - r2)} )
  return(rout)
}



#### Path definitions ####

# define input file paths
image_files_evi <- list('veg_indices/evi/20170509_evi.tif', 'veg_indices/evi/20180514_evi.tif')
treemask_file <- 'tree_mask/tree_mask_species_map_4_inverse_buffered.tif'
swamp_file <- 'swamp_shapefiles/swamp_clipped.shp'

# define output file paths
outpaths_hmg <- list('texture_metrics/from_evi/hmg_2017_rescaled_to_eightbit_n256_any.tif',
                     'texture_metrics/from_evi/hmg_2018_rescaled_to_eightbit_n256_any.tif')
outpaths_var <- list('texture_metrics/from_evi/var_2017_rescaled_to_eightbit_n256_any.tif',
                     'texture_metrics/from_evi/var_2018_rescaled_to_eightbit_n256_any.tif')

# open files the same for all years
treemask <- open.band(treemask_file) # trees where treemask = 1
swamp_shp <- readOGR(paste0(wd, swamp_file))


#### Process Data ####

# loop over years
for (i in 1:length(image_files_evi)) {
  # read in evi data
  evi_image <- open.band(image_files_evi[i])
  
  # repare Data for GLCM calculation
  evi_image_prepped <- evi_image %>%
    mask.raster.with.shp(swamp_shp) %>%
    mask.raster.with.raster(treemask) %>%
    rescale.outliers() %>%
    rescale.to.eightbit
  
  # calculate texture metrics
  tic('glcm calculation, 3x3 window')
  evi_glcm <- calc.textures(evi_image_prepped, n_grey = 256)
  toc()

  # write out rasters
  writeRaster(evi_glcm$glcm_homogeneity, filename = paste0(wd, outpaths_hmg[i]),
              format = 'GTiff', overwrite = TRUE)
  writeRaster(evi_glcm$glcm_variance, filename = paste0(wd, outpaths_var[i]),
              format = 'GTiff', overwrite = TRUE)
  
  print
}


#### Difference rasters ####

# do calculation for homogeneity and variance
diff_wet_dry_hmg <- diff_rasters(evi_2018$glcm_homogeneity, evi_2017$glcm_homogeneity)
diff_wet_dry_var <- diff_rasters(evi_2018$glcm_variance, evi_2017$glcm_variance)


#### Write out

# define output file names
outpath_hmg_diff <- 'texture_metrics/from_evi/diff_wet2018_minus_dry2017_hmg_rescaled_to_eightbit_n256_any.tif'
outpath_var_diff <- 'texture_metrics/from_evi/diff_wet2018_minus_dry2017_var_rescaled_to_eightbit_n256_any.tif'

# write out rasters
# NB overwrites files if they already exist
writeRaster(diff_wet_dry_hmg, filename = paste0(wd, outpath_hmg_diff),
            format = 'GTiff', overwrite = TRUE)
writeRaster(diff_wet_dry_var, filename = paste0(wd, outpath_var_diff),
            format = 'GTiff', overwrite = TRUE)