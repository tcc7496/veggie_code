'''
A script containing functions to calculate vegetation indices from sentinel 2 level-2A images
'''

from cloudmask import all
from tree_mask import all
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.tools.mask import mask
import fiona


def ndvi(filepath, aoi = None):
    '''
    A function to calculate ndvi = (B8 - B4) / (B8 + B4)

    Parameters
    ----------
    filepath: to folder containing individual bands
    aoi: optional parameter to crop raster to an area of interest by providing full path to a shapefile

    Returns
    ----------
    '''

    # %%
    # for testing
    filepath = '/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/S2A_MSIL2A_20180514T073611_N9999_R092_T36MZC_20210520T183726.SAFE/GRANULE/L2A_T36MZC_A015104_20180514T075528/IMG_DATA/R10m'

    if aoi is not None:
        with fiona.open("/Users/Cate/UK_Mainland.shp", "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]
            # can do it with the crop=True thing I think

        with rasterio.open("jan_clip.tif") as src:
            out_image, out_transform = mask(src, geoms, crop=True)
            out_meta = src.meta.copy()

        out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

        with rasterio.open("masked2.tif", "w", **out_meta) as dest:
            dest.write(out_image)
    


    # read in red and NIR bands
    red = rio.open(f'{filepath}/T36MZC_20180514T073611_B04_10m.jp2', driver ='JP2OpenJPEG')
    nir = rio.open(f'{filepath}/T36MZC_20180514T073611_B08_10m.jp2', driver ='JP2OpenJPEG')

    # read metadata
    kwargs = red.meta

    # read values to do band maths
    b04 = red.read(1, masked = True).astype('float64') # change dtype to float in prep for ndvi calc
    b08 = nir.read(1, masked = True).astype('float64')

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

    # write out ndvi test
    with rio.open('/Users/taracunningham/projects/dissertation/sen2processing/processing/l2a/ndvi_test2.tif', 'w', **kwargs) as dst:
        dst.write_band(1, ndvi)


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
