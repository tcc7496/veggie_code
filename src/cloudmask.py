'''
A script to create a cloud mask of Sentinel 2 L1C images using python fmask
'''

import rasterio as rio
import numpy as np
import os
import subprocess
import glob
import fiona
from rasterio.mask import mask

#######################################

def cloud_fmask(file, outdir):
    '''
    A function to produce a cloud mask using fmask.
    Input file is TOA reflectance .SAFE directory
    Output is geotiff of cloud mask
    '''
    # obtain full file path to input file
    filename = os.path.basename(file)

    # check if output directory exists and create it if not
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # make temporary directory to store img files created by fmask commandline run
    # construct path to img directory and create directory if it doesn't already exist
    imgpath = os.path.join(outdir, 'img')
    if not os.path.exists(imgpath):
        os.mkdir(imgpath)

    # construct output file name
    outpath = f'{os.path.join(imgpath, filename[:-5])}_cloud.img'

    # run cloud mask from command line
    subprocess.run(f'fmask_sentinel2Stacked.py -o {outpath} --safedir {file} --pixsize 10 --parallaxtest', shell=True)

    # convert output .img file to geotiff
    raster_to_tif(outpath, outdir)

    # delete temperory img directory
    subprocess.run(f'rm -rf {imgpath}', shell=True)

#######################################

def batch_cloud_fmask(inputdir, outdir):
    '''
    A function to batch process a directory of images to create their cloud masks
    '''

    # create list of .SAFE directories in input directory
    inputfilelist = glob.glob(f'{inputdir}/*.SAFE')

    if len(inputfilelist) == 0:
        print(f'Input directory does not contain any .SAFE directories.')
    else:
        # loop over list of files
        for file in inputfilelist:
            print(f'Processing {file}...')
            cloud_fmask(file, outdir)


#######################################

def raster_to_tif(file, outdir):
    '''
    A function to convert a raster to geotiff format
    If output directory isn't specified, file is put in the same folder as input file
    '''
    # separate filename from filepath
    filename = os.path.basename(file)

    # find number of characters to remove at end of filename
    nchars = len(filename.split('.')[-1])
    # if nchars is not zero, +1 to remove the dot
    if nchars > 0:
        nchars = nchars + 1

    # think can also use os.path.splitext(filename).

    # run conversion
    # check if output directory specified
    if outdir == None:
        subprocess.run(f'gdal_translate -of Gtiff {file} {file[:-nchars]}.tif', shell=True)
    else:
        subprocess.run(f'gdal_translate -of Gtiff {file} {os.path.join(outdir, filename)[:-nchars]}.tif', shell=True)

    print(f'{filename} converted to geotiff')

#######################################

def fmask_to_boolean_cloudmask(file, aoi = None):
    '''
    A function to convert fmask output to boolean cloudmask array.

    Parameters
    ----------
    file: tif file output by fmask algorithm
    aoi: optional to crop file to an area of interest

    Returns
    ----------
    Boolean array of the cloudmask where only clear pixels are used. fmask = 1 -> mask = False. masked pixel -> True
    '''
    
    # read in result of fmask
    if aoi is None:
        with rio.open(file) as src:
            fmask = src.read(1)
            profile = src.profile.copy()

    else:
        # open shapefile to crop raster to
        with fiona.open(aoi, "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]
        
        # open band
        with rio.open(file) as src:
            fmask, transform = mask(src, geoms, crop = True, filled = True)
            fmask = fmask[0]
            profile = src.profile.copy()
        # update profile for new shape
        profile.update({
                 "height": fmask.shape[0],
                 "width": fmask.shape[1],
                 "transform": transform
                 })

    # convert to boolean by masking everywhere except where fmask = 1 = clear
    fmask_clear = np.isin(fmask, 1, invert = True)

    return fmask_clear, profile

#######################################

def rasterized_polygon_clouds_to_cloudmask(file, aoi = None):
    '''
    A function to convert rasterized cloud polygons to boolean cloudmask array.

    Parameters
    ----------
    file: tif file where the data indicates cloud.
    aoi: area of interest to crop data to

    Returns
    ----------
    Boolean array of the cloudmask where only clear pixels -> False and cloud pixel -> True
    '''

    if aoi is None:
        # read in mask. Reads in GDAL format => nodata = 0
        with rio.open(file) as src:
            msk = src.read_masks(1)
        
        # convert to boolean
        boolean_msk = np.isin(msk, 0, invert = True)

    else:
        with fiona.open(aoi, "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]
        
        # open band
        with rio.open(file) as src:
            msk, transform = mask(src, geoms, crop = True, filled = True)
            msk = msk[0]    # cloud = 1
            profile = src.profile.copy()
            # update profile for new shape
            profile.update({
                 "height": msk.shape[0],
                 "width": msk.shape[1],
                 "transform": transform
                 })
        # convert to boolean
        boolean_msk = np.isin(msk, 1)

    return boolean_msk

    
#######################################

if __name__ == "__main__":
    ''' Main block '''

    


