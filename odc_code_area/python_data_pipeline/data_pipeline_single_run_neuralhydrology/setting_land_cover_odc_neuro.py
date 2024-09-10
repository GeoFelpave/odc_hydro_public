# -*- coding: utf-8 -*-
"""==============================================================================
Title          :05_climate_land_cover_odc.py
Description    :Create the land cover map needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 03_setting_land_cover_odc.py
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
    
    log.info(f'\n --- \n Reading Land Cover from ODC\n --- \n')
    
    # Set ODC application
    dc = datacube.Datacube(app="worldcover")

    # Get data
    # Load data from the datacube
    buffer = 0.125
    ds = dc.load(product=dc_product,
                lat=(ymn - buffer, ymx + buffer),
                lon=(xmn - buffer, xmx + buffer),
                time=('2021'),
                )

    # Print output data
    return ds

def lc_catchm_values(wrk_path:str, dc_array, 
                     path_mask:str, df_geo, 
                     uid_value, nd_value):
    
    log.info(f'\n --- \n Setting land cover attributes\n --- \n')
    
    # Create temp tif file with ODC data
    _tif = Path(wrk_path / 'temp_data/ea_worldcover.tif')
    dc_array['classification'].rio.to_raster(_tif)
    
    # Find the pixel with the most valus within each grid
    # This produces a dict following the same order as the grid geodataframe
    # with the pixel count
    catchm_stats_lc = rasterstats.zonal_stats(str(path_mask), str(_tif), 
                                              categorical=True,
                                              nodata=nd_value)
    
    # Change count to percentage
    stats_dict = catchm_stats_lc[0]
    pixel_total = sum(stats_dict.values())
    for k, v in stats_dict.items():
        stats_dict[k] = round(v * 100.0 / pixel_total,2)
        # print(round(v * 100.0 / sum(catchm_stats_lc[0].values()),2))
    
    
    # change land cover number to name
    lc_dict = {10:'tree_cover', 20:'shrubland', 30:'grassland',
               40:'cropland', 50:'built-up', 60:'sparse_veget',
               70:'snow_ice', 80:'perm_water', 90:'herb_wetlnd',
               95:'mangroves', 100:'moss_lichen'}
    
    # change land cover from number to categorical name
    stats_dict = dict([(lc_dict.get(k), v) for k, v in stats_dict.items()])
        
    # Set dictionary as dataframe
    df_lc = pd.DataFrame(stats_dict, index=[uid_value])
    
    # Add dominant land cover
    df_lc['dom_land_cover'] = df_lc.idxmax(axis=1)
    
    # Add suffix 
    df_lc = df_lc.rename(columns={col: f'{col}_perc' for col in 
                                  df_lc.columns if col not in 
                                  ['dom_land_cover']})
    
    # tidy the final dataframe
    # this is necessary to match neuralhydrology requirements for CAMELS GB
    df_lc.index.rename('gauge_id', inplace=True)
    df_lc = df_lc.reset_index()
    
    # Delete temp .tif file
    utils.file_remove(_tif.parent, 'tif')

    return df_lc

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("\n --- \n Creating land cover attributes files\n --- \n")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")
@click.option("--unique_value", 
              default=None, 
              help="Enter the unique identifier for the catchment")

def lc_map(unique_field, unique_value):
    
    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')
    
    # get iso code
    iso_code = config.get('iso_code', 'country')

    # Open Data Cube Product
    dc_data = config.get('dc_product', 'LC')

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
    df_final = lc_catchm_values(dir_abs, array_odc, shp_path, _shp, unique_value, ND)
    
    # save to csv
    df_final.to_csv(Path(dir_abs /
                        f'temp_data/ODC_{iso_code}_'\
                        'land_cover_attributes_'\
                        f'{unique_value}.csv'), 
                    index=False)
    
    log.info(f'\n --- \n Land cover attributes have been created\n --- \n')
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    lc_map()
    