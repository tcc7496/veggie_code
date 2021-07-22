'''
A script to mask three tree species in a plant species map of the study area. Outputs a tif where nodata = 0 and data = 255 in GDAL RFC 15 style mask
'''
# %%
import rasterio as rio
from rasterio.enums import Resampling
from matplotlib import pyplot as plt
import numpy as np
from band_process import *
import shutil
from rasterio.io import MemoryFile

# %%
def tree_mask(file, outfile, aoi = None):
    '''
    A function outputs a tif representing tree species masked
    '''
    # %%

    aoi = '/Users/taracunningham/projects/dissertation/other_data/study_area_shapefile/study_area.geojson'

    tmp = shutil.copy(file, f'{file[:-4]}_copy.tif')

    # resample data to 10m resolution
    upscale_factor = 2

    with rio.open(tmp, mode = 'r+') as src:
        # resample data to target shape
        species_map = src.read(
        out_shape=(
            src.count,
            int(src.height * upscale_factor),
            int(src.width * upscale_factor)
            ),
        resampling=Resampling.nearest
        )

        # copy profile
        profile = src.profile.copy()

        # scale image transform
        transform = src.transform * src.transform.scale(
            (src.width / species_map.shape[-1]),
            (src.height / species_map.shape[-2])
            )

        # update profile
        profile.update({
                 "height": species_map.shape[1],
                 "width": species_map.shape[2],
                 "transform": transform
                 })

        # write src to new dataset
        with MemoryFile() as memfile:
            with memfile.open(**profile) as mem:
                mem.write(species_map)

                if aoi is not None:
        
                    # open shapefile to crop raster to
                    with fiona.open(aoi, "r") as shapefile:
                        geoms = [feature["geometry"] for feature in shapefile]
        
                    species_map, transform = mask(mem, geoms, crop = True, filled = True)
                    species_map = species_map[0]

                    # update profile for new shape
                    profile.update({
                        "height": species_map.shape[0],
                        "width": species_map.shape[1],
                        "transform": transform
                        })

    # convert to masked array with species desired masked
    species_trees_masked = np.ma.array(species_map, mask = np.isin(species_map, [2, 4, 5, 6, 8, 9], invert = True))
    
    # set fill value to universal nodata value
    np.ma.set_fill_value(species_trees_masked, -9999)

    # update profile
    update_profile_nodata(profile)

    # write out tree mask as tif
    with rio.open(outfile, 'w', **profile) as dst:
        dst.write(species_trees_masked.filled(), 1)

    
#######################################

def tree_mask_bool(file, aoi = None):
    '''
    a function that outputs a boolean array indicating tree species to be masked.

    Parameters
    ----------
    file: tif file with different plant species indicated by integers.

    Returns
    ----------
    tree_mask_bool: boolean array of the treemask where species masked -> True
    '''

    # read in species_map
    if aoi == None:
        with rio.open(file) as src:
            species_map = src.read(1)
            tree_mask = src.read_masks(1)

    else:
        # open shapefile to crop raster to
        with fiona.open(aoi, "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]
        
        # open band
        with rio.open(file) as src:
            tree_mask, transform = mask(src, geoms, crop = True, filled = True)
            tree_mask = tree_mask[0]
            profile = src.profile.copy()
        # update profile for new shape
        profile.update({
                 "height": fmask.shape[0],
                 "width": fmask.shape[1],
                 "transform": transform
                 })

    # create starting mask with nodata values masked
    tree_mask_bool = np.isin(tree_mask, 0)

    # create list of species numbers to mask
    msk_species = [1, 3, 7]

    # mask each species in 
    for m in msk_species:
        msk = np.isin(species_map, m)
        tree_mask_bool = msk | tree_mask_bool

    

    return tree_mask_bool

    






# %%
if __name__ == "__main__":
    ''' Main block '''

    file = '/Users/taracunningham/projects/dissertation/other_data/original/svmRadial_Multiple_season_time_series_raster_final.tif'

    outfile = '/Users/taracunningham/projects/dissertation/sen2processing/processing/tree_mask/tree_mask_species_map_4.tif'

    aoi = '/Users/taracunningham/projects/dissertation/other_data/study_area_shapefile/study_area.geojson'

    tree_mask(file, outfile, aoi = aoi)

# %%
