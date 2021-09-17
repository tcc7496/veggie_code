'''
A script to calcuate vegetation indices EVI and NDVI from atmospherically corrected Sentinel 2 images with cloud masking.

'''

#######################################

from veg_indices import *
import click

@click.command()
@click.argument('inputdir', type=click.Path(exists=True), help = 'directory containing Level 2A Sentinel 2 images to be processed')
@click.argument('path_to_cloudmasks', type=click.Path(exists=True), help ='full filepath to text file containing full filepaths to the corresponding cloud masks for each file')
@click.option('-a', '--aoi', default=None, type=click.Path(), help='optional area of interest shapefile')
def main(inputdir, path_to_cloudmasks, aoi):
    batch_veg_indices(inputdir, path_to_cloudmasks, aoi)
    

#######################################

if __name__ == "__main__":
    ''' Main block '''

    main()

#######################################