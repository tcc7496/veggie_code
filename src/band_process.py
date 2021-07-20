'''
A script to deal with tif files
'''

import rasterio as rio
from rasterio.mask import mask
import fiona
import pyproj
from shapely.ops import transform

def open_band(file, aoi = None):
    '''
    A function to open a sen2 raster file into memory
    
    Parameters
    ----------
    file: full filepath to raster file
    aoi: optional full filepath to .shp or geojson to crop raster to an area of interest

    Returns
    ----------
    masked numpy array where values to be masked are zeroes
    mask
    profile

    '''

    if aoi is None:
        with rio.open(file, driver = 'JP2OpenJPEG') as src:
            image = src.read(1, masked = True)
            msk = src.read_masks(1)
            profile = src.profile.copy()
           
    else:
        # open shapefile to crop raster to
        with fiona.open(aoi, "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]
        
        # open band
        with rio.open(file, driver = 'JP2OpenJPEG') as src:
            image, transform = mask(src, shapes = geoms, crop = True, filled = False)
            image = image[0]
            msk = src.read_masks(1)
            profile = src.profile.copy()
        # update profile for new shape
        profile.update({
                 "height": image.shape[0],
                 "width": image.shape[1],
                 "transform": transform
                 })
            
    return image, msk, profile

#######################################

def check_crs():
    '''
    check two crs are the same. Return boolean
    '''


#######################################

def reproject_crs(file, dst_crs):
    '''

    '''


if __name__ == "__main__":
    ''' Main block '''

    
