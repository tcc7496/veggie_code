'''
A script to extract the tree mask from the plant species map

'''

#######################################

from tree_mask import *
import click

@click.command()
@click.argument('file', type=click.Path(exists=True), help = 'full filepath to input plant species map')
@click.argument('outfile', type=click.Path(), help ='full filepath including filename of output geotiff')
@click.option('-a', '--aoi', default=None, type=click.Path(), help='optional area of interest shapefile')
def main(file, outfile, aoi):
    tree_mask(file, outfile, aoi)
    

#######################################

if __name__ == "__main__":
    ''' Main block '''

    main()

#######################################