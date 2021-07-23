'''
A script for the full workflow from original data to end products

'''

from cloudmask import *
from tree_mask import *
from band_process import *
from veg_indices import *


#######################################
'''
1 Atmsopheric Correction:
level-1C product to level-2A. Run script 'atm_cor.sh' from command line 

'''
#######################################
'''
2 Cloud Masking:

Fmask algorithm run on all files.
If fmask is unsatisfactory, cloud masks manually drawn in qgis to be used later

'''
inputdir = '/Users/taracunningham/projects/dissertation/sen2processing/original/sen2/'
outdir = '/Users/taracunningham/projects/dissertation/sen2processing/processing/fmask/'

batch_cloud_fmask(inputdir, outdir)

#######################################
'''
3 Tree Masking:

Species map provided - mask tree species 1, 3, 7

'''

file = '/Users/taracunningham/projects/dissertation/other_data/original/svmRadial_Multiple_season_time_series_raster_final.tif'

outfile = '/Users/taracunningham/projects/dissertation/sen2processing/processing/tree_mask/tree_mask_species_map_4.tif'

aoi = '/Users/taracunningham/projects/dissertation/other_data/study_area_shapefile/study_area.geojson'

tree_mask(file, outfile, aoi = aoi)

#######################################
'''
4 Calculate veg indices - NDVI, EVI:

NDVI calculation for masking brown vegetation. Threshold still tbc. Starting point 0.3 - 0.4.

'''

# construct output names for ndvi files



ndvi = ndvi_calc(filepath, cloudmask, outfile = outfile_ndvi, aoi = aoi)