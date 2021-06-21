'''
A script to mask tree species in a species map. A masked tif if output
'''
# %%
import rasterio as rio
from matplotlib import pyplot
import numpy as np

# %%
def tree_mask(file, outfile):
    '''
    A function that masks trees in the input species map and outputs a new tif with trees masked.
    '''
    
    with rio.open(file) as src:
        # open with mask. Nodata values are True
        species_map_masked = src.read(1, masked=True)

        # read data into separate array to combine with new mask
        species_map = src.read(1)

        # create starter mask with nodata values
        tree_mask = species_map_masked.mask

        # copy profile for writing out
        profile = src.profile.copy()

    # set nodata value
    profile['nodata'] = -9999
    
    # mask where species_map is:
    # 1 Acacia | 3 Ficus Sur | 7 Prosopis
    tree_mask = np.where((species_map == 1) | (species_map == 3) | (species_map == 7), True, tree_mask)

    # create new masked array of species map and tree mask
    species_map_tree_mask = np.ma.masked_array(species_map, mask = tree_mask, fill_value=profile['nodata'])

    # write out species map with tree mask
    with rio.open(outfile, 'w', **profile) as dst:
        dst.write(species_map_tree_mask.filled(profile['nodata']), 1)
    

#######################################
# %%
if __name__ == "__main__":
    ''' Main block '''

    file = '/Users/taracunningham/projects/dissertation/other_data/original/svmRadial_Multiple_season_time_series_raster_final.tif'
    outfile = '/Users/taracunningham/projects/dissertation/sen2processing/processing/tree_mask/tree_mask_species_map.tif'

    tree_mask(file, outfile)

# %%
