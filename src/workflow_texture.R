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
open_band <- function(band) {
  band <- raster(paste0(wd, band))
  # calc and save min and max values of raster to object
  band <- setMinMax(band)
  return(band)
}


##### 2017 #####

### Read in data

# define file paths
image_file_evi <- 'veg_indices/evi/20170509_evi.tif'
treemask_file <- 'tree_mask/tree_mask_species_map_4_inverse_buffered.tif'
swamp_file <- 'swamp_shapefiles/swamp_clipped.shp'

# read in data
evi_image <- open_band(image_file_evi)
treemask <- open_band(treemask_file) # trees where treemask = 1
swamp_shp <- readOGR(paste0(wd, swamp_file))



### Masking

# swamp
evi_image_masked <- mask(evi_image, mask = swamp_shp, inverse = TRUE)

# trees
evi_image_masked[!is.na(treemask)] <- NA


### Prepare data

# set values less than 0 to zero
# ref: Farwell et al. pg.3
evi_image_masked[evi_image_masked < 0] <- 0

# set values > 1 to NA because meaningless - cloud or water. Negligible pixel count.
# don't have a reference for this really.
evi_image_masked[evi_image_masked > 1] <- NA

# rescale data and change to integers
evi_image_masked <- as.integer(evi_image_masked*100)



#### GLCM calculation

tic('glcm calculation, 3x3 window')
evi_2017 <- glcm(evi_image_masked[[1]], window = c(3,3), n_grey = 101,
                 shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                 statistics = c('homogeneity', 'variance'),
                 min_x = 0, max_x = 100,
                 na_opt = 'center', na_val = NA)
toc()


##### 2018 #####

image_file_evi <- 'veg_indices/evi/20180514_evi.tif'
evi_image <- open_band(image_file_evi)

# mask
evi_image_masked <- mask(evi_image, mask = swamp_shp, inverse = TRUE)
evi_image_masked[!is.na(treemask)] <- NA

# prepare data
evi_image_masked[evi_image < 0] <- 0
evi_image_masked[evi_image > 1] <- NA

# rescale
evi_image_masked <- as.integer(evi_image_masked*100)

tic('glcm calculation, 3x3 window')
evi_2018 <- glcm(evi_image_masked[[1]], window = c(3,3), n_grey = 101,
                 shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                 statistics = c('homogeneity', 'variance'),
                 min_x = 0, max_x = 100,
                 na_opt = 'center', na_val = NA)
toc()


#### Subtract rasters

# define function
diff_rasters <- function(r1, r2) {
  rout <- overlay(r1, r2,
                  fun = function(r1, r2) {return(r1 - r2)} )
  return(rout)
}

# do calculation
diff_wet_dry_hmg <- diff_rasters(evi_2018$glcm_homogeneity, evi_2017$glcm_homogeneity)
diff_wet_dry_var <- diff_rasters(evi_2018$glcm_variance, evi_2017$glcm_variance)


#### Write out

# define output file names
outpath_hmg_2017 <- 'texture_metrics/from_evi/diff_2017_hmg_rescaled_n100.tif'
outpath_var_2017 <- 'texture_metrics/from_evi/diff_2017_var_rescaled_n100.tif'
outpath_hmg_2018 <- 'texture_metrics/from_evi/diff_2018_hmg_rescaled_n100.tif'
outpath_var_2018 <- 'texture_metrics/from_evi/diff_2018_var_rescaled_n100.tif'
outpath_hmg_diff <- 'texture_metrics/from_evi/diff_wet2018_minus_dry2017_hmg_rescaled_n100.tif'
outpath_var_diff <- 'texture_metrics/from_evi/diff_wet2018_minus_dry2017_var_rescaled_n100.tif'

# write out
writeRaster(evi_2017$glcm_homogeneity, filename = paste0(wd, outpath_hmg_2017),
            format = 'GTiff', overwrite = TRUE)
writeRaster(evi_2017$glcm_variance, filename = paste0(wd, outpath_var_2017),
            format = 'GTiff', overwrite = TRUE)
writeRaster(evi_2018$glcm_homogeneity, filename = paste0(wd, outpath_hmg_2018),
            format = 'GTiff', overwrite = TRUE)
writeRaster(evi_2018$glcm_variance, filename = paste0(wd, outpath_var_2018),
            format = 'GTiff', overwrite = TRUE)
writeRaster(diff_wet_dry_hmg, filename = paste0(wd, outpath_hmg_diff),
            format = 'GTiff', overwrite = TRUE)
writeRaster(diff_wet_dry_var, filename = paste0(wd, outpath_var_diff),
            format = 'GTiff', overwrite = TRUE)