'''
A script to create a cloud mask and mask atmospherically corrected (2A) sentinel 2 images
'''

# %%
import rasterio as rio
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
from osgeo import gdal


#######################################

def cloud_fmask(filedir, file, outdir):
    '''
    function to produce a cloud mask using fmask.
    Input file is TOA reflectance .SAFE directory
    Output is geotiff of cloud mask
    '''
    # obtain full file path to input file
    filepath = os.path.join(filedir, file)

    # check if output directory exists and create it if not
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # construct output file name
    outpath = f'{os.path.join(outdir, file[:-5])}_cloud.img'

    # run cloud mask from command line
    subprocess.run(f'fmask_sentinel2Stacked.py -o {outpath} --safedir {filepath}', shell=True)

    # convert output .img file to geotiff
    raster_to_tif(outpath)

#######################################

def raster_to_tif(file):
    '''
    convert raster to geotiff format
    '''
    # find number of characters to remove at end of filename
    nchars = len(file.split('.')[-1])
    # if nchars is not zero, +1 to remove the dot
    if nchars > 0:
        nchars = nchars + 1
    # run conversion
    subprocess.run(f'gdal_translate -of Gtiff {file} {file[:-nchars]}.tif', shell=True)
    
    print(f'{file} converted to geotiff')


#######################################

if __name__ == "__main__":
    ''' Main block '''
# %%
    # hardcode file path for testing
    filepath = '/Users/taracunningham/projects/dissertation/sen2processing/original/sen2/'
    file = 'S2A_MSIL1C_20200523T073621_N0209_R092_T36MZC_20200523T095431.SAFE'
    outdir = '/Users/taracunningham/projects/dissertation/sen2processing/processing/fmask/'

    #cloud_fmask(filepath, file, outdir)
    outpath = f'{os.path.join(outdir, file[:-5])}_cloud.img'
    raster_to_tif(outpath)
    


# %%
