'''
A script to deal with tif files
'''

import rasterio as rio
from rasterio.mask import mask
import fiona
import numpy as np

#######################################

def open_band(file, aoi = None):
    '''
    A function to open a sen2 raster file into memory
    
    Parameters
    ----------
    file: full filepath to raster file
    aoi: optional full filepath to .shp or geojson to crop raster to an area of interest

    Returns
    ----------
    masked numpy array, nodata values masked = True and set to -9999 in data
    mask
    profile

    '''

    if aoi is None:
        with rio.open(file, driver = 'JP2OpenJPEG') as src:
            image = src.read(1, masked = True)
            #msk = src.read_masks(1)
            profile = src.profile.copy()
           
    else:
        # open shapefile to crop raster to
        with fiona.open(aoi, "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]
        
        # open band
        with rio.open(file, driver = 'JP2OpenJPEG') as src:
            image, transform = mask(src, shapes = geoms, crop = True, filled = True) # non-masked array, 0 outside polygons
            image = image[0]
            #msk = src.read_masks(1)
            profile = src.profile.copy()
        # update profile for new shape
        profile.update({
                 "height": image.shape[0],
                 "width": image.shape[1],
                 "transform": transform
                 })

    # update array and profile dtype to int16 to enable it to store -9999
    image = image.astype(np.int16)
    profile.update({
        "dtype": 'int16'
    })

    # convert to masked array
    image = np.ma.array(image, mask = np.isin(image, 0))

    # update masks to nodata masks
    #image.mask = np.isin(image.data == profile['nodata'], True)

    # update nodata value in profile
    update_profile_nodata(profile)

    # update fill value
    np.ma.set_fill_value(image, -9999)


    return image, profile
    # msk return removed

#######################################

def update_profile_nodata(profile):
    '''
    A function to update the nodatavalue to the same thing throughout processing
    '''

    profile.update({
        "nodata": -9999.0
    })

#######################################

def update_nodata_vals(arr, nodataval_in, arr_mask = None):
    '''
    A function to update masked array values to -9999 in the data
    '''
    if arr_mask == None:
        arr_mask = arr
    
    arr[arr_mask == nodataval_in] = -9999

#######################################

def update_profile_dtype(profile, dst_dtype):
    '''
    A function to update profile to universal dtype of float64
    '''

    profile.update({
        "dtype": dst_dtype
    })

#######################################


if __name__ == "__main__":
    ''' Main block '''

    

