######################
# Script containing functions for processing texture metrics
######################

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
  r_masked <- raster::mask(r, mask = shp, inverse = TRUE)
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

calc.zonal.stats <- function(
  r, polygon,
  stats = c('min', 'max', 'count', 'mean', 'median', 'stdev', 'quantile', 'coefficient_of_variation')) {
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

remove.hmg.outliers <- function (inputdir, probs) {
  # A function to rescale the percentiles specified by probs in hmg rasters
  # and write out the new rescaled file.
  # output file will be the same as input file with _wo_outliers added to the end
  inputfilelist <- Sys.glob(file.path(inputdir, '*.tif' )) %>%
    str_sort(numeric = T)
  
  # create list of output files
  outfilelist <- inputfilelist %>%
    basename() %>%
    str_remove_all(pattern = c('.tif')) %>%
    paste0(inputdir, ., '_wo_outliers.tif')
  
  # loop over input files
  for (i in 1:length(inputfilelist)) {
    # open tif
    r <- open.band(inputfilelist[i])
    # rescale outliers
    r_rescaled <- rescale.outliers.probs(r, probs = probs)
    # write out tif
    writeRaster(r_rescaled,
                filename = outfilelist[i],
                format = 'GTiff', overwrite = TRUE)
    
  }
}

# if function above is run, don't need to have the rescale_outliers bit in this. Just read in the
# intermediate rasters
calc.zonal.stats.batch <- function (inputdir, outfile, aoi, rescale_outliers = c(0.02, 0.98)) {
  ## Inputs
  # inputdir: directory with input tifs. Will try and process all tifs
  # outfile: name of output file. Will be placed in same directory as input files
  # aoi: path to shapefile of aoi
  # rescale_outliers: vector with percentage intervals to rescale outliers.
  #                   Defaults to no rescaling
  
  # get list of input files
  inputfilelist <- Sys.glob(file.path(inputdir, "*.tif")) %>%
    str_sort(numeric = T)
  
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
  inputdir, outfile, aoi, same_pixels = T, normalise = F) {
  ## Inputs
  # inputdir: path to input directory of raster files
  # outfile: path to the dir output to be stored including name of file
  # aoi: path to shapefile containing polygons to extract values within
  # same_pixels: If true, selects only pixels that are not NA for every year
  # normalise: If true, normalises pixels to the mean value of Landuse over all years
  
  ## Outputs
  # writes a csv of the resulting dataframe to the output directory.
  # returns the dataframe also
  
  # read in shapefile
  polys <- st_read(aoi) %>%
    dplyr::select('LandUse', 'geometry')
  
  # obtain list of input files
  inputfilelist = Sys.glob(file.path(inputdir, "*n256.tif")) %>%
    str_sort(numeric = T)
  
  # create stack of files
  s <- raster::stack(inputfilelist)
  
  # rescale outliers
  # apply rescale.outliers function to one layer at a time and restack result
  s <- stack(lapply(1:nlayers(s), function(i){rescale.outliers.probs(s[[i]])}))
  
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
    
    # write out mean values per landuse
    #write.csv2(mean_per_landuse, file = 'landuse_normalised_values.csv')
    
    # Normalise values of homogeneity to calculated mean values / landuse
    values_per_landuse_tidy_cp <- values_per_landuse_tidy_cp %>%
      pivot_longer(cols =  starts_with("x"), names_to = "year") %>%
      dplyr::mutate(value_norm =
                      dplyr::case_when(LandUse == mean_per_landuse$LandUse[1]
                                       ~ value / mean_per_landuse$mean[1],
                                       LandUse == mean_per_landuse$LandUse[2]
                                       ~ value / mean_per_landuse$mean[2])
      )
  } else {
    # pivot df longer before completing rest of calculation if normalise = F
    values_per_landuse_tidy_cp <- values_per_landuse_tidy_cp %>%
      pivot_longer(cols =  starts_with("x"), names_to = "year")
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
  
  # remove x in front of years and change to numeric
  values_per_landuse_tidy_cp_na <- values_per_landuse_tidy_cp_na %>%
    tidyr::separate(col = year, into = c(NA, 'year'), sep = 1)
  values_per_landuse_tidy_cp_na$year <- as.numeric(values_per_landuse_tidy_cp_na$year)
  
  # write out dataframe
  write.csv2(values_per_landuse_tidy_cp_na, file = outfile)
  
  return(values_per_landuse_tidy_cp_na)
}

calc.normalised.hmg.raster <- function (inputdir, outputdir, statsfile, meanfile) {
  # A function that outputs a raster normalised to the mean of each polygon across all years
  # takes in rasters with hmg outliers (0.02, 0.98) removed
  
  # read in zonal stats
  stats_sf <- st_read(paste0(inputdir, statsfile))
  
  # read in normalised values per landuse
  mean_df <- read.csv2(paste0(inputdir, meanfile))
  
  # get list of evi files in year
  inputfilelist <- Sys.glob(file.path(inputdir, '*wo_outliers.tif' )) %>%
    str_sort(numeric = T)
  # create list of output files
  outfilelist <- inputfilelist %>%
    basename() %>%
    str_remove_all(pattern = c('.tif')) %>%
    paste0(outputdir, ., '_norm2.tif')
  
  # create year list
  year_list <- c('2021', '2020', '2019', '2018', '2017', '2016') %>%
    str_sort(numeric = T)
  
  # read in mean over years
  # [1,3]: Conservancy, [2,3]: Livestock Rearing Area
  mean_df <- read.csv2(paste0(inputdir, meanfile))
  
  # loop over inputfile list
  for (i in 1:length(inputfilelist)) {
    
    # open raster
    r <- open.band(inputfilelist[i])
    
    # create copy of raster for each landuse
    r_cp_1 <- r
    r_cp_2 <- r
    s <- raster::stack(r_cp_1, r_cp_2)
    
    ## Conservancy
    c_poly <- as_Spatial(stats$geom[1])
    # mask Livestock Rearing area
    s[[1]] <- mask(s[[1]], c_poly)
    # divide by mean to normalize
    s[[1]] <- s[[1]]/mean_df[1,3]
    # change NA values to zeros for addition
    s[[1]][is.na(s[[1]])] <- 0
    
    
    ## Livestock Rearing Area
    g_poly <- as_Spatial(stats$geom[2])
    # mask conservancy
    s[[2]] <- mask(s[[2]], g_poly)
    # divide by mean to normalize
    s[[2]] <- s[[2]]/mean_df[2,3]
    # change NA values to zeros for addition
    s[[2]][is.na(s[[2]])] <- 0
    
    # sum rasters
    r_new <- sum(s) # I think
    
    # set zero values back to NA
    r_new[r_new == 0] <- NA
    
    # write out raster
    writeRaster(r_new,
                filename = outfilelist[i],
                format = 'GTiff', overwrite = TRUE)
    
  }
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
  df$year <- as.numeric(df$year)
  
  return(df)
}

calc.rainfall <- function (inputdir, outfile) {
  
  year_list <- c('2021', '2020', '2019', '2018', '2017', '2016') %>%
    str_sort(numeric = T)
  
  # open shapefile to get crs
  study_area_shp <- st_read(study_area)
  
  # create empty df for results
  rainfall_df <- data.frame(LandUse = character(),
                            weighted_mean = double(),
                            year = character())
  
  # loop over years
  for (i in 1:length(year_list)) {
    
    # get files for one year
    inputfilelist <- Sys.glob(file.path(inputdir, glue('*{year_list[i]}*.tif') ) )
    
    # stack raster files
    s <- raster::stack(inputfilelist)
    
    # reproject rasterstack to study area crs
    s_reproj <- projectRaster(s, crs = crs(study_area_shp))
    
    # crop rasterstack to study area
    s_crop <- crop(s_reproj, study_area_shp, snap = "out")
    
    # sum raster stack
    rain_crop_sum <- calc(s_crop, fun = sum)
    
    # extract average rainfall per polygon
    rain_mean <- cbind(
      exact_extract(rain_crop_sum, study_area_shp,
                    fun = 'weighted_mean',
                    weights = 'area',
                    force_df = TRUE,
                    #include_cols = c('LandUse'),
                    progress = T),
      LandUse = study_area_shp$LandUse) %>%
      dplyr::mutate('year' = year_list[i]) # add year as a column
    
    # add it to zs_df
    rainfall_df <- rbind(rainfall_df, rain_mean)
  }
  
  # change name of weight_mean column
  names(rainfall_df)[names(rainfall_df) == "weighted_mean"] <- "rain"
    
  # write out data
  write.csv2(rainfall_df, file = outfile)
    
  return(rainfall_df)

}