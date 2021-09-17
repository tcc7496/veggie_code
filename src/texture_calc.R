#
# A script to calculate texture metrics
#

#### Libraries ----
library(glcm)
library(raster)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(sp)
library(tictoc)     # for timing operations
library(viridis)
library(tidyr)
library(sf)
library(exactextractr)
library(glue)
library(purrr)
library(stringr)
library(tibble)

## Working Dir ----

# set working directory
wd <- "/Users/taracunningham/projects/dissertation/sen2processing/processing/"
setwd(wd)

## Functions ----

source('../src/texture_functions.R')


# Calculate homogeneity texture metric from EVI ---

# define filepaths
evi_fp <- file.path('veg_indices', 'evi') # location of evi files to process
treemask_file <- 'tree_mask/tree_mask_species_map_4_inverse_buffered.tif'
swamp_file <- 'swamp_shapefiles/swamp_clipped.shp'
study_area <- '../../other_data/study_area_shapefile/study_area_buffer_aggr_w_conservancy.shp'
outdir <- file.path('texture_metrics', 'from_evi', 'hmg', 'abs')

# define other inputs for glcm calculation
textures = c('hmg')
n_grey = 256

# call function to calculate texture
calc.textures.from.evi.batch(
  textures = textures,
  n_grey = 256,
  evi_fp = evi_fp,
  aoi_fp = study_area,
  treemask_fp = treemask_file,
  swamp_fp = swamp_file,
  outdir = outdir
  )

## Rescale outliers in homogeneity tifs and write new rasters ----
# remove homogeneity data outside the 98th percentile

inputdir = 'texture_metrics/from_evi/hmg/abs/'
probs = c(0.02, 0.98)   # percentiles to remove

# function call
remove.hmg.outliers(inputdir = inputdir, probs = probs)


## Zonal Stats ----
# Calculate zonal statistics per land use for each year and write out geopackage

hmg_fp <- 'texture_metrics/from_evi/hmg/abs/'
outfile <- paste0(inputdir, 'zs_from_abs_evi.gpkg')

# call function to calculate zonal stats
calc.zonal.stats.batch(inputdir = hmg_fp,
                       outfile = outfile,
                       aoi = study_area
                       )


## Extract all raster values per polygon and save to csv ----

# Do 2 runs, one with the only the same pixels across all the years, one with all available pixels
# If normalise = T, output csv contains absolute and normalised values values
# Also outputs a csv file landuse_mean_values.csv which contains the mean value of homogeneity across all years
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


## Calculate normalised homogeneity rasters ----

# define filepaths
inputdir = 'texture_metrics/from_evi/hmg/abs/'
outputdir = 'texture_metrics/from_evi/hmg/norm/'
statsfile <- 'zs_from_evi_abs.gpkg'               # homogeneity zonal stats per land use for all years
meanfile <- '../landuse_mean_values.csv'    # mean values per landuse across all years

# function call
calc.normalised.hmg.raster(inputdir, outputdir, statsfile, meanfile)


## Perform Mann-Whitney U test ----
# equivalent to Wilcoxon Rank Sum test

# import homogeneity data
hmg_fp <- 'texture_metrics/from_evi/hmg/norm/'
inputfile <- paste0(hmg_fp, 'hmg_norm_same_pixels_F.csv')

hmg_df <- read.csv2(inputfile)

# Wilcoxon test
# formula: homogeneity ~ land use

# for all years
wil_test_all_yrs <- wilcox.test(
  hmg_df$value ~ hmg_df$LandUse,
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)

# within each year
wil_test_2016 <- wilcox.test(
  hmg_df$value[hmg_df$year == 2016] ~ hmg_df$LandUse[hmg_df$year == 2016],
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)

wil_test_2017 <- wilcox.test(
  hmg_df$value[hmg_df$year == 2017] ~ hmg_df$LandUse[hmg_df$year == 2017],
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)

wil_test_2018 <- wilcox.test(
  hmg_df$value[hmg_df$year == 2018] ~ hmg_df$LandUse[hmg_df$year == 2018],
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)

wil_test_2019 <- wilcox.test(
  hmg_df$value[hmg_df$year == 2019] ~ hmg_df$LandUse[hmg_df$year == 2019],
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)

wil_test_2020 <- wilcox.test(
  hmg_df$value[hmg_df$year == 2020] ~ hmg_df$LandUse[hmg_df$year == 2020],
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)

wil_test_2021 <- wilcox.test(
  hmg_df$value[hmg_df$year == 2021] ~ hmg_df$LandUse[hmg_df$year == 2021],
  alternative = 'two.sided',
  conf.int = T,
  digits.rank = 4)
