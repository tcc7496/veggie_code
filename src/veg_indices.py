'''
A script containing functions to calculate vegetation indices from sentinel 2 level-2A images
'''
# %%
from cloudmask import fmask_to_boolean_cloudmask
from tree_mask import tree_mask_bool
from band_process import open_band, update_nodata_vals
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.mask import mask
import fiona
import glob
import os
from rasterio.plot import show


def ndvi_calc(filepath, cloudmask, outfile = None, aoi = None):
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
    outfile: optional write out a tif of ndvi

    Returns
    ----------
    ndvi: masked numpy array of ndvi
    out_profile: profile for ndvi array

    '''

    # prepare boolean cloudmask from fmask
    clouds, _ = fmask_to_boolean_cloudmask(cloudmask, aoi = aoi)

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
        bands_data[band]['image'], bands_data[band]['profile'] = open_band(f, aoi = aoi)

        # change datatype so that band maths work
        bands_data[band]['image'] = bands_data[band]['image'].astype(np.float32)

        bands_data[band]['profile']['dtype'] = 'float32'    # update profile


    # compute full mask from OR of both bands and cloud mask
    total_mask = bands_data['b04']['image'].mask | bands_data['b08']['image'].mask | clouds

    # assign new mask to both band masks
    bands_data['b04']['image'].mask = bands_data['b08']['image'].mask = total_mask

    # calculate ndvi
    ndvi = (bands_data['b08']['image'] - bands_data['b04']['image']) / (bands_data['b08']['image'] + bands_data['b04']['image'])

    # double check setting fill value to -9999
    np.ma.set_fill_value(ndvi -9999)

    # create profile for writing out
    out_profile = bands_data['b04']['profile']

    out_profile.update(
        dtype = ndvi.dtype,
        driver = 'Gtiff'
        )

    if outfile is not None:
        # write out ndvi test
        with rio.open(outfile, 'w', **out_profile) as dst:
            dst.write_band(1, ndvi.filled())

    return ndvi, out_profile

    '''
    rasters = [rasterio.open(p) for p in file_paths]
    stacked = np.dstack([r.read() for r in rasters])
    '''

#######################################

def mask_ndvi(ndvi, threshold):
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

    with rio.open(ndvi) as src:
        image = src.read(1, masked = True)

    ndvi_mask = np.ma.masked_where(image < threshold, image, copy = True).mask
    
    return ndvi_mask

#######################################

def evi_calc(filepath, cloudmask, treemask = None, ndvimask = None, outfile = None, aoi = None):
    '''
    A function to calculate EVI
    '''

    # prepare boolean cloudmask from fmask
    clouds, _ = fmask_to_boolean_cloudmask(cloudmask, aoi = aoi)

    # read in tree mask
    with rio.open(treemask) as src:
        trees = src.read(1, masked = True).mask

    # construct file path to individual bands
    bands_path = os.path.join(filepath, 'GRANULE/*/IMG_DATA/R10m/')

    # get list of full filepaths to bands 4 and 8
    bands_list = [glob.glob(b)[0] for b in [f'{bands_path}*B02*', f'{bands_path}*B04*', f'{bands_path}*B08*']]

    # create list of bands to loop over
    bands_to_process = ['b02', 'b04', 'b08']

    # create empty dictionary to store outputs
    bands_data = {}

    # loop over bands list
    for band, f in zip(bands_to_process, bands_list):
        # add band to dictionary
        bands_data[band] = {}
        # add outputs to dictionary
        bands_data[band]['image'], bands_data[band]['profile'] = open_band(f, aoi = aoi)

        # change datatype so that band maths work
        bands_data[band]['image'] = bands_data[band]['image'].astype(np.float32)

        bands_data[band]['profile']['dtype'] = 'float32'    # update profile

        # rescale band data for evi calculation
        bands_data[band]['image'] = bands_data[band]['image'] / 10000


    # compute full mask from OR of both bands and other masks
    if treemask is None:
        if ndvimask is None:
            total_mask = bands_data['b02']['image'].mask | bands_data['b04']['image'].mask | bands_data['b08']['image'].mask | clouds
        else:
            total_mask = bands_data['b02']['image'].mask | bands_data['b04']['image'].mask | bands_data['b08']['image'].mask | clouds | ndvimask
    else:
        if ndvimask is None:
            total_mask = bands_data['b02']['image'].mask | bands_data['b04']['image'].mask | bands_data['b08']['image'].mask | clouds | trees
        else:
            total_mask = bands_data['b02']['image'].mask | bands_data['b04']['image'].mask | bands_data['b08']['image'].mask | clouds | trees | ndvimask

    # assign new mask to both band masks
    bands_data['b02']['image'].mask = bands_data['b04']['image'].mask = bands_data['b08']['image'].mask = total_mask

    # calculate evi
    # 2.5(NIR - R) / (1 + NIR + 6R - 7.5B) = 2.5(b08 - b04) / (1 + b08 + 6*b04 - 7.5*b02)
    evi = 2.5*(bands_data['b08']['image'] - bands_data['b04']['image']) / (1 + bands_data['b08']['image'] + 6*bands_data['b04']['image'] - 7.5*bands_data['b02']['image'])
    plt.imshow(evi)

    # double check setting fill value to -9999
    np.ma.set_fill_value(evi, -9999)

    # create profile for writing out
    out_profile = bands_data['b02']['profile']

    out_profile.update(
        dtype = evi.dtype,
        driver = 'Gtiff'
        )

    if outfile is not None:
        # write out ndvi test
        with rio.open(outfile, 'w', **out_profile) as dst:
            dst.write_band(1, evi.filled())

    return evi, out_profile




if __name__ == "__main__":
    ''' Main block '''

    # %%
    # for testing
    filepath = '/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/S2A_MSIL2A_20180514T073611_N9999_R092_T36MZC_20210520T183726.SAFE/'

    cloudmask = '/Users/taracunningham/projects/dissertation/sen2processing/processing/fmask/S2A_MSIL1C_20180514T073611_N0206_R092_T36MZC_20180514T095515_cloud.tif'

    treemask = '/Users/taracunningham/projects/dissertation/sen2processing/processing/tree_mask/tree_mask_species_map_4.tif'

    outfile_ndvi = '/Users/taracunningham/projects/dissertation/sen2processing/processing/ndvi/T36MZC_20180514_ndvi.tif'

    aoi = '/Users/taracunningham/projects/dissertation/other_data/study_area_shapefile/study_area.geojson'

    ndvi = '/Users/taracunningham/projects/dissertation/sen2processing/processing/ndvi/T36MZC_20180514_ndvi.tif'

    outfile_ndvi_mask = '/Users/taracunningham/projects/dissertation/sen2processing/processing/ndvi/T36MZC_20180514_ndvi_mask_0-3.tif'

    outfile_evi = '/Users/taracunningham/projects/dissertation/sen2processing/processing/evi/T36MZC_20180514_evi_ndvi_mask_0-3.tif'

    threshold = 0.3


    #ndvi, ndvi_profile = ndvi_calc(filepath, cloudmask, outfile = outfile_ndvi, aoi = aoi)

    ndvi_mask = mask_ndvi(ndvi, threshold)

    evi, evi_profile = evi_calc(filepath, cloudmask, treemask, ndvimask = ndvi_mask, outfile = outfile_evi, aoi = aoi)



    #plt.imshow(ndvi_mask, cmap='gray')
    #plt.imshow(ndvi)








# %%
