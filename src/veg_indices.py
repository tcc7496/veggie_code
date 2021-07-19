'''
A script containing functions to calculate vegetation indices from sentinel 2 level-2A images
'''

from cloudmask import all
from tree_mask import all
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.mask import mask
import fiona
import glob
from band_process import open_band


def ndvi(filepath, aoi = None):
    '''
    A function to calculate ndvi = (B8 - B4) / (B8 + B4)

    Parameters
    ----------
    filepath: to SAFE directory for L2A
    aoi: optional parameter to crop raster to an area of interest by providing full path to a shapefile

    Returns
    ----------
    '''
    # %%
    # for testing
    filepath = '/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/S2A_MSIL2A_20180514T073611_N9999_R092_T36MZC_20210520T183726.SAFE/'
    # construct file path to individual bands
    bands_path = os.path.join(filepath, 'GRANULE/*/IMG_DATA/R10m/')
    # get list of full filepaths to bands 4 and 8
    bands_list = [glob.glob(b)[0] for b in [f'{bands_path}*B04*', f'{bands_path}*B08*']]


    # %%
    # create empty dictionary to store outputs
    bands_data = {}

    # create list of bands to loop over
    bands_to_process = ['b04', 'b08']

    # loop over bands list
    for band, f in zip(bands_to_process, bands_list):
        # add band to dictionary
        bands_data[band] = {}
        # add outputs to dictionary
        bands_data[band]['image'], bands_data[band]['mask'], bands_data[band]['profile'] = open_band(f, aoi = None)

        # convert masks to boolean masks
        bands_data[band]['image'].mask = np.isin(bands_data[band]['image'].data, 0)
        bands_data[band]['mask'] = bands_data[band]['image'].mask

        # change datatype so that band maths work
        bands_data[band]['image'] = bands_data[band]['image'].astype('float64')
        bands_data[band]['profile']['dtype'] = 'float64'    # update profile

    # %%

    # calculate ndvi
    ndvi = (bands_data['b08']['image'] - bands_data['b04']['image']) / (bands_data['b08']['image'] + bands_data['b04']['image'])

    # create profile for writing out
    out_profile = bands_data['b04']['profile']
    out_profile.update(
        dtype=ndvi.dtype,
        driver = 'Gtiff'
        )


    # write out ndvi test
    with rio.open('/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/ndvi_test3.tif', 'w', **kwargs) as dst:
        dst.write_band(1, ndvi)

    # %%
    '''
    rasters = [rasterio.open(p) for p in file_paths]
    stacked = np.dstack([r.read() for r in rasters])
    '''

#######################################

def check_crs():
    '''
    check two crs are the same. Return boolean
    '''

#######################################

def reproject_crs(file, dst_crs):
    '''
    '''
    


# %%
