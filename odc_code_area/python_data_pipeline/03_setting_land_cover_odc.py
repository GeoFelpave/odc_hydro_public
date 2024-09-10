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
                     nd_value):
    
    log.info(f'\n --- \n Assigning Land Cover values to each grid\n --- \n')
    
    # Create temp tif file with ODC data
    _tif = Path(wrk_path / 'shetran_data/temp_data/ea_worldcover.tif')
    dc_array['classification'].rio.to_raster(_tif)
    
    # Find the pixel with the most valus within each grid
    # This produces a dict following the same order as the grid geodataframe
    catchm_stats_lc = rasterstats.zonal_stats(str(path_mask), str(_tif), categorical=True, nodata=nd_value)

    # Create a list containing the dictionary key (land cover type) with the Max value for each grid in the catchment
    lc_largest = [max(lc_grid, key=lc_grid.get) for lc_grid in catchm_stats_lc]

    # Add the largest lc to grid geodataframe
    df_geo['lc_largest'] = lc_largest
    
    # Replacing 0 with -9999
    df_geo.loc[df_geo['lc_largest'] == 0, 'lc_largest'] = nd_value
    
    # Delete temp .tif file
    utils.file_remove(_tif.parent, 'tif')

    return df_geo

def shetran_lc(df_geo, crs_lcl, wrk_path):
    
    log.info(f'\n --- \n Creating SHETRAN Land Cover map\n --- \n')
    
    # Get the centroid for each grid
    grid_bg = df_geo.to_crs(crs_lcl)
    grid_bg['centroid'] = grid_bg['geometry'].centroid

    # Create lat (Y) and lon (X) columns
    grid_bg['lat'] = grid_bg['centroid'].y.astype(int)
    grid_bg['lon'] = grid_bg['centroid'].x.astype(int)

    # Change geopandas to pandas ready to create csv file for the mean and min elevation
    df_grid_lc = pd.DataFrame(grid_bg[['lc_largest', 'lat', 'lon']].copy())

    # Pivoting dataframe to replicate SHETRAN format
    # Pivoting dataframe using lon as column and lat as row
    df_pivot_lc = df_grid_lc.pivot(index='lat', columns='lon', values='lc_largest')

    # Remove any NaN
    df_pivot_lc = df_pivot_lc.fillna(0)


    # Creating the text file for min elevation
    utils.shetran_csv_file(wrk_path, 'final_land_cover_SHETRAN', df_pivot_lc, 'd')
    log.info(f'\n --- \n SHETRAN Land Cover file created\n --- \n')

    

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating SHETRAN Land Cover file\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating land cover map")
# @click.option("--ctchmnt_path", 
#               default=None, 
#               help="Enter the complete path to the catchment shapefile")

def lc_map():
    
    # Setting the path to the work environment
    dir_abs = Path().resolve()

    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'COL')

    # Open Data Cube Product
    dc_data = config.get('dc_product', 'LC')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    # Read shp file to geopandas dataframe
    grid_path = Path(dir_abs / 'shetran_data/active_data/final_mask_wgs84.shp')
    grid = gpd.read_file(grid_path)

    xmin, ymin, xmax, ymax = grid.total_bounds
    
    # Getting data from ODC
    array_odc = read_odc(dc_data, xmin, ymin, xmax, ymax )
    
    # DEM values for SHETRAN mask
    grid_DEM = lc_catchm_values(dir_abs, array_odc, grid_path, grid, ND)
    
    # Creating frinal map
    shetran_lc(grid_DEM, crs_local, dir_abs)
    
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    lc_map()
    