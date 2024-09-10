# -*- coding: utf-8 -*-
"""==============================================================================
Title          :04_setting_lakes_odc.py
Description    :Create the lakes map needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 04_setting_lakes_odc.py
Notes          :
python version :
=============================================================================="""
# =============================================================================
# Packages - Libraries
# =============================================================================
from pathlib import Path
import geopandas as gpd
import pandas as pd
import numpy as np
import logging
import click

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

def grid_latlon(grid_df, crs_lcl):
    # Get the centroid for each grid
    temp_df = pd.DataFrame()
    temp_df = grid_df.to_crs(crs_lcl)
    temp_df['centroid'] = temp_df['geometry'].centroid

    # Create lat (Y) and lon (X) columns
    temp_df['lat'] = temp_df['centroid'].y.astype(int)
    temp_df['lon'] = temp_df['centroid'].x.astype(int)

    return temp_df


def lakes_catchm_values(df_grid, df_lakes,
                        crs_lcl, nd_value,
                        wrk_path):
    
    # Intersect grid with lakes to find which grid has lakes
    grid_lakes = gpd.sjoin(df_grid, df_lakes, how="left", predicate="intersects")

    # Replace NaN with empty sting
    grid_lakes = grid_lakes.replace(np.nan, '')


    # Assing LAKE ID if value GLWD_ID is empty
    grid_lakes["LAKE_ID"] = np.where(grid_lakes['Hylak_id'] == '', nd_value, 1)

    # Create lat and lot values using the function
    grid_centroids = grid_latlon(grid_lakes, crs_lcl)

    # Change geopandas to pandas ready to create csv file
    df_grid_lakes = pd.DataFrame(grid_centroids[['LAKE_ID', 'lat', 'lon']].copy())

    # Remove duplicates just in case there is an error
    df_grid_lakes = df_grid_lakes.drop_duplicates()

    # Pivoting dataframe to replicate SHETRAN format
    # Pivoting dataframe using lon as column and lat as row
    df_pivot_grid_lakes = df_grid_lakes.pivot(index='lat', columns='lon', values='LAKE_ID')
    df_pivot_grid_lakes = df_pivot_grid_lakes.sort_index(ascending=False)
    
    # Creating the text file for min elevation
    utils.shetran_csv_file(wrk_path, 'final_lake_map_SHETRAN', df_pivot_grid_lakes, 'd')
    log.info(f'\n --- \n SHETRAN Lake file created\n --- \n')

    

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating SHETRAN Lake file\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating lakes map")

def lakes_map():
    
    # Setting the path to the work environment
    dir_abs = Path().resolve()
    
    # Read shp file to geopandas dataframe
    grid_path = Path(dir_abs / 'shetran_data/active_data/final_mask_wgs84.shp')
    grid = gpd.read_file(grid_path)

    lakes_path = Path(dir_abs / 'shetran_data/active_data/HydroLAKES_polys_v10.shp')
    lakes = gpd.read_file(lakes_path)

    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'COL')

    # Open Data Cube Product
    dc_data = config.get('dc_product', 'LC')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    # Lakes SHETRAN mask
    lakes_catchm_values(grid, lakes,
                        crs_local, ND,
                        dir_abs)
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    lakes_map()
    