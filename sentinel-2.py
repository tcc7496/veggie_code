from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
from datetime import date
import pandas as pd
import os

s1dir = '/Users/robwebster/Sync/msc_course/dissertation/data/sentinel-1/'
s2dir = '/Users/robwebster/Sync/msc_course/dissertation/data/sentinel-2/'
l8dir = '/Users/robwebster/Sync/msc_course/dissertation/data/landsat-8/'
aoidir = '/Users/robwebster/Sync/msc_course/dissertation/data/aoi/'

aoi_json = '/Users/robwebster/gsi/kenya_study_area.geojson'
#print(aoi_json)

# connect to the API
api = SentinelAPI(None, None)

# search by polygon, time, and SciHub query keywords
aoi = geojson_to_wkt(read_geojson(aoi_json))
from_date = '20210101'
to_date = 'NOW'


def s1(aoi, from_date, to_date, download):
    '''
    Finds online Sentinel-1 imagery
    '''

    products = api.query(aoi, date=(from_date, to_date), platformname='Sentinel-1', producttype='GRD', sensoroperationalmode='IW')

    # Creates GeoPandas GeoDataFrame with the metadata of the scenes and the footprints as geometries
    geodf = api.to_geodataframe(products)

    print(f'\nSentinel-1 Images from {from_date} to {to_date}:\n')
    
    if not geodf.empty:
        for id in geodf['uuid']:
            print(f"    Product: {api.get_product_odata(id)['title']}  (Online: {api.get_product_odata(id)['Online']})")

    else:
        print('     No valid data in this time period\n')

    if download:
        api.download_all(products, s1dir)

    return geodf

def s2(aoi, from_date, to_date, download):
    '''
    Finds online Sentinel-2 imagery
    '''

    products = api.query(aoi, date=(from_date, to_date), platformname='Sentinel-2', cloudcoverpercentage=(0, 30), producttype='S2MSI1C')

    # Creates GeoDataFrame
    df = api.to_dataframe(products)
    print(df.columns)

    print(f'\nSentinel-2 Images from {from_date} to {to_date}:\n')
    
    if not df.empty:
        for id in df['uuid']:
            print(f"    Product: {api.get_product_odata(id)['title']}  (Online: {api.get_product_odata(id)['Online']})")

    else:
        print('\tNo valid data in this time period\n')

    if download:
        api.download_all(products, s2dir)   

    return df


download = False
#s1(aoi, from_date, to_date, download)
s2(aoi, from_date, to_date, download)