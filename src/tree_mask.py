'''
A script to mask three tree species in a plant species map of the study area. Outputs a tif where nodata = 0 and data = 255 in GDAL RFC 15 style mask
'''
# %%
import rasterio as rio
from matplotlib import pyplot
import numpy as np

# %%
def tree_mask(file, outfile):
    '''
    A function outputs a tif representing tree species masked
    '''
    
    with rio.open(file) as src:
        # load species map to see where trees are
        species_map = src.read(1)

        # open with mask. nodata = 0
        tree_mask = src.read_masks(1)

        # copy profile for writing out
        profile = src.profile.copy()

    # set nodata value to 0 to be within int8 range
    profile['nodata'] = 0

    # set profile dtype to tree_mask dtype
    profile['dtype'] = tree_mask.dtype

    # mask where species_map is:
    # 1 Acacia | 3 Ficus Sur | 7 Prosopis
    tree_mask = np.where((species_map == 1) | (species_map == 3) | (species_map == 7), 0, tree_mask)

    # write out tree mask as tif
    with rio.open(outfile, 'w', **profile) as dst:
        dst.write(tree_mask, 1)
    
#######################################

def tree_mask_bool(file):
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
    with rio.open(file) as src:
        species_map = src.read(1)
        tree_mask = src.read_masks(1)

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
    outfile = '/Users/taracunningham/projects/dissertation/sen2processing/processing/tree_mask/tree_mask_species_map.tif'

    tree_mask(file, outfile)

# %%
