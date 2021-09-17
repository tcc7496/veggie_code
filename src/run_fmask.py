'''
A script to batch run fmask on Sentinel 2 Level-1C images

'''

#######################################

from cloudmask import *
import click

@click.command()
@click.argument('inputdir', type=click.Path(exists=True), help = 'directory containing Sentinel 2 .SAFE directories to be processed')
@click.argument('outdir', type=click.Path(), help ='output directory to place output file')
def main(inputdir, outdir):
    batch_cloud_fmask(inputdir, outdir)
    

#######################################

if __name__ == "__main__":
    ''' Main block '''

    main()

#######################################

