# -*- coding: utf-8 -*-
"""==============================================================================
Title          :02_setting_dem_odc.py
Description    :Create the min and mean DEM maps needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 02_setting_dem_odc_neuro.py
Notes          :
python version :
=============================================================================="""
# =============================================================================
# Packages - Libraries
# =============================================================================
from pathlib import Path
import geopandas as gpd
import pandas as pd
import rasterstats
import rioxarray
import logging
import click

# ODC modules
import datacube

# Custom modules
import utils

# Modules for config.ini
from configparser import ConfigParser
config = ConfigParser()

#############################################################
# FUNCTIONS
#############################################################
def setup_logging(level: int = logging.INFO) -> logging.Logger:
    """Set up a simple logger -
    https://engineeringfordatascience.com/posts/python_logging/"""
    log = logging.getLogger(__name__)
    console = logging.StreamHandler()
    log.addHandler(console)
    log.setLevel(level)
    return log

def read_odc(dc_product:str, xmn:int, ymn:int, 
             xmx:int, ymx:int):
    
    log.info(f'\n --- \n Reading DEM from ODC\n --- \n')
    
    # Set ODC application
    dc = datacube.Datacube(app="dem_cop")

    # Get data
    # Load data from the datacube
    buffer = 0.125
    ds = dc.load(product=dc_product,
                lat=(ymn - buffer, ymx + buffer),
                lon=(xmn - buffer, xmx + buffer),
                )

    # Print output data
    return ds

def dem_catchm_values(wrk_path:str, dc_array, 
                      path_mask:str, df_geo,
                      nd_value, uid_field):
    
    log.info(f'\n --- \n Finding all elevation values for the catchment\n --- \n')
    
    # Create temp tif file with ODC data
    _tif = Path(wrk_path / 'temp_data/dem_cop.tif')
    dc_array['elevation'].rio.to_raster(_tif)
    
    # Get elevation mean, min, max values
    stats_lst = ['mean', 'min', 'max', 'percentile_10',
                 'percentile_50', 'percentile_90']
    
    # Find the mean elevation for each grid in the catchment
    # This produces a dict following the same order as the grid geodataframe
    catchm_stats = rasterstats.zonal_stats(str(path_mask), str(_tif), stats=stats_lst, all_touched=True, nodata=nd_value)

    # # Change the dictionary to a list to be added to the grid geodataframe
    df_stats = pd.DataFrame.from_dict(catchm_stats)
    
    # Add the elevation values to geodataframe
    df_geo = pd.concat([df_geo, df_stats], axis=1)
    
    # tidy the final dataframe
    # this is necessary to match neuralhydrology requirements for CAMELS GB
    col_names = stats_lst.insert(0, uid_field)
    df_geo_final = df_geo[stats_lst].copy()
    df_geo_final.rename(columns={uid_field:'gauge_id'}, inplace=True)

    # Delete temp .tif file
    utils.file_remove(_tif.parent, 'tif')


    return df_geo_final

    

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating topographic attributes files")
# @click.option("--ctchmnt_path", 
#               default=None, help="Enter the complete path to the catchment shapefile")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")
@click.option("--unique_value", 
              default=None, 
              help="Enter the unique identifier for the catchment")

def dem_attributes(unique_field, unique_value):
    
    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')
    
    # get iso code
    iso_code = config.get('iso_code', 'country')

    # Open Data Cube Product
    dc_data = config.get('dc_product', 'DEM')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    # Setting the path to the work environment
    dir_abs = Path(Path().resolve() / f'ODC_{iso_code}')
    
    # Read shp file to geopandas dataframe
    shp_path = Path(dir_abs / 
                f'active_data/{unique_field}_{unique_value}_WGS84.shp')
    
    _shp = gpd.read_file(shp_path)

    xmin, ymin, xmax, ymax = _shp.total_bounds
    
    # Getting data from ODC
    array_odc = read_odc(dc_data, xmin, ymin, xmax, ymax )
    
    # DEM values for SHETRAN mask
    df_DEM = dem_catchm_values(dir_abs, array_odc, 
                               shp_path, _shp, 
                               ND, unique_field)
    
    # save to csv
    df_DEM.to_csv(Path(dir_abs / 
                    f'temp_data/ODC_{iso_code}_'\
                    'topographic_attributes_'\
                    f'{unique_value}.csv'), 
                    index=False)
    
    log.info(f'\n --- \n Topographic attributes have been created\n --- \n')
        
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    dem_attributes()
    