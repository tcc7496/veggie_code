'''
A script containing functions to calculate vegetation indices from sentinel 2 level-2A images
'''

from cloudmask import all
from tree_mask import all
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio


def ndvi(filepath, aoi):
    '''
    A function to calculate ndvi = (B8 - B4) / (B8 + B4)

    Parameters
    ----------
    filepath to 

    Returns
    ----------
    '''

    # %%
    # for testing
    filepath = '/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/S2A_MSIL2A_20180514T073611_N9999_R092_T36MZC_20210520T183726.SAFE/GRANULE/L2A_T36MZC_A015104_20180514T075528/IMG_DATA/R10m'

    # read in red and NIR bands
    red = rio.open(f'{filepath}/T36MZC_20180514T073611_B04_10m.jp2', driver ='JP2OpenJPEG')
    nir = rio.open(f'{filepath}/T36MZC_20180514T073611_B08_10m.jp2', driver ='JP2OpenJPEG')

    # read metadata
    kwargs = red.meta

    # read values to do band maths
    b04 = red.read(1, masked = True).astype('float64')
    b08 = nir.read(1, masked = True)

    # close files
    red.close()
    nir.close()

    # calculate ndvi
    ndvi = (b08 - b04) / (b08 + b04)

    # update profile for writing out
    kwargs.update(
        dtype=ndvi.dtype,
        driver = 'Gtiff'
        )

# %%
    # write out ndvi test
    with rio.open('/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/ndvi_test.tif', 'w', **kwargs) as dst:
        dst.write_band(1, ndvi)

    # values not below zero anywhere. check stats. Divide by 10000? Shouldn't matter because of the ratio. Only matters in because formula adds numbers.

    '''
    rasters = [rasterio.open(p) for p in file_paths]
    stacked = np.dstack([r.read() for r in rasters])
    '''
    


# %%
