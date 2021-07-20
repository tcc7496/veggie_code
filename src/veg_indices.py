'''
A script containing functions to calculate vegetation indices from sentinel 2 level-2A images
'''

from cloudmask import fmask_to_boolean_cloudmask
from tree_mask import tree_mask_bool
from band_process import open_band
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.mask import mask
import fiona
import glob
import os


def ndvi(filepath, cloudmask, outfile = None, aoi = None):
    '''
    A function to calculate ndvi = (B8 - B4) / (B8 + B4)

    Parameters
    ----------
    filepath: to SAFE directory for L2A
    cloudmask: full filepath to cloudmask to be used. tif or shapefile
    outfile: 
    aoi: optional parameter to crop raster to an area of interest by providing full path to a shapefile

    Outputs
    ----------
    Optional write out a tif of ndvi

    Returns
    ----------
    ndvi: masked numpy array of ndvi
    out_profile: profile for ndvi array

    '''

    # prepare boolean cloudmask from fmask
    clouds, clouds_profile = fmask_to_boolean_cloudmask(cloudmask, aoi = aoi)

    # construct file path to individual bands
    bands_path = os.path.join(filepath, 'GRANULE/*/IMG_DATA/R10m/')

    # get list of full filepaths to bands 4 and 8
    bands_list = [glob.glob(b)[0] for b in [f'{bands_path}*B04*', f'{bands_path}*B08*']]

    # create list of bands to loop over
    bands_to_process = ['b04', 'b08']

    # create empty dictionary to store outputs
    bands_data = {}

    # loop over bands list
    for band, f in zip(bands_to_process, bands_list):
        # add band to dictionary
        bands_data[band] = {}
        # add outputs to dictionary
        bands_data[band]['image'], bands_data[band]['mask'], bands_data[band]['profile'] = open_band(f, aoi = aoi)

        # convert mask to boolean mask
        bands_data[band]['image'].mask = np.isin(bands_data[band]['image'].data, 0)
        # add cloud mask with OR
        bands_data[band]['image'].mask = bands_data[band]['image'].mask | clouds
        # set mask variable
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

    if outfile is not None:
        # write out ndvi test
        with rio.open(outfile, 'w', **out_profile) as dst:
            dst.write_band(1, ndvi)

    return ndvi, out_profile

    # %%
    '''
    rasters = [rasterio.open(p) for p in file_paths]
    stacked = np.dstack([r.read() for r in rasters])
    '''

#######################################

def ndvi_mask(ndvi, threshold):
    '''
    function to produce a boolean mask below a threshold value of ndvi
    
    Parameters
    ----------
    ndvi: masked array of ndvi
    threshold: pixels below this value will be masked


    Returns
    ----------
    ndvi_mask: boolean mask with pixels below the threshold value masked 

    '''

    ndvi_mask = np.ma.masked_where(ndvi < threshold, ndvi, copy = True).mask
    
    return ndvi_mask




    

if __name__ == "__main__":
    ''' Main block '''

    # for testing
    filepath = '/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/S2A_MSIL2A_20180514T073611_N9999_R092_T36MZC_20210520T183726.SAFE/'

    cloudmask = '/Users/taracunningham/projects/dissertation/sen2processing/processing/fmask/S2A_MSIL1C_20180514T073611_N0206_R092_T36MZC_20180514T095515_cloud.tif'

    outfile_ndvi = '/Users/taracunningham/projects/dissertation/sen2processing/processing/ndvi/T36MZC_20180514_ndvi.tif'

    aoi = '/Users/taracunningham/projects/dissertation/other_data/study_area_shapefile/study_area.geojson'

    outfile_ndvi_mask = '/Users/taracunningham/projects/dissertation/sen2processing/processing/ndvi/T36MZC_20180514_ndvi_mask_0-3.tif'

    threshold = 0.3


    ndvi, ndvi_profile = ndvi(filepath, cloudmask, outfile = outfile_ndvi, aoi = aoi)

    ndvi_mask = ndvi_mask(ndvi, threshold)

    plt.imshow(ndvi_mask)
    plt.imshow(ndvi)







