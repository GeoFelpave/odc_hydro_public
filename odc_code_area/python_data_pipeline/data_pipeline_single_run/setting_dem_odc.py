# -*- coding: utf-8 -*-
"""==============================================================================
Title          :02_setting_dem_odc.py
Description    :Create the min and mean DEM maps needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 02_setting_dem_odc.py
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
                      nd_value):
    
    log.info(f'\n --- \n Assigning DEM values to each grid\n --- \n')
    
    # Create temp tif file with ODC data
    _tif = Path(wrk_path / 'temp_data/dem_cop.tif')
    dc_array['elevation'].rio.to_raster(_tif)
    
    # Find the mean elevation for each grid in the catchment
    # This produces a dict following the same order as the grid geodataframe
    catchm_stats_mean = rasterstats.zonal_stats(str(path_mask), str(_tif), stats="mean", all_touched=True, nodata=nd_value)
    catchm_stats_min = rasterstats.zonal_stats(str(path_mask), str(_tif), stats="min", all_touched=True, nodata=nd_value)

    # Change the dictionary to a list to be added to the grid geodataframe
    elevation_mean = [x['mean'] for x in catchm_stats_mean]
    elevation_min = [x['min'] for x in catchm_stats_min]

    # Add the elevation mean and min to grid geodataframe
    df_geo['elevation_mean'] = elevation_mean
    df_geo['elevation_min'] = elevation_min
    
    # Delete temp .tif file
    utils.file_remove(_tif.parent, 'tif')

    return df_geo

def shetran_dem(df_geo, crs_lcl, wrk_path,
                uid_field, uid_value):
    
    log.info(f'\n --- \n Creating SHETRAN DEM map\n --- \n')
    
    # Get the centroid for each grid
    grid_bg = df_geo.to_crs(crs_lcl)
    grid_bg['centroid'] = grid_bg['geometry'].centroid

    # Create lat (Y) and lon (X) columns
    grid_bg['lat'] = grid_bg['centroid'].y.astype(int)
    grid_bg['lon'] = grid_bg['centroid'].x.astype(int)

    # Change geopandas to pandas ready to create csv file for the mean and min elevation
    df_grid_mean = pd.DataFrame(grid_bg[['elevation_mean', 'lat', 'lon']].copy())
    df_grid_min = pd.DataFrame(grid_bg[['elevation_min', 'lat', 'lon']].copy())

    # Pivoting dataframe to replicate SHETRAN format
    # Pivoting dataframe using lon as column and lat as row
    df_pivot_mean = df_grid_mean.pivot(index='lat', columns='lon', values='elevation_mean')
    df_pivot_mean = df_pivot_mean.sort_index(ascending=False).round(0).astype(int)

    df_pivot_min = df_grid_min.pivot(index='lat', columns='lon', values='elevation_min')
    df_pivot_min = df_pivot_min.sort_index(ascending=False).round(0).astype(int)
    
    # Creating the text file for min elevation
    utils.shetran_csv_file(wrk_path, f'{uid_field}_{uid_value}_dem_min', df_pivot_min, 'd')
    log.info(f'\n --- \n {uid_field}_{uid_value} min elevation file created\n --- \n')

    utils.shetran_csv_file(wrk_path, f'{uid_field}_{uid_value}_dem_mean', df_pivot_mean, 'd')
    log.info(f'\n --- \n {uid_field}_{uid_value} mean elevation file created\n --- \n')
    
    
    

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating SHETRAN DEM file\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating catchment mask")
# @click.option("--ctchmnt_path", 
#               default=None, help="Enter the complete path to the catchment shapefile")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")
@click.option("--unique_value", 
              default=None, 
              help="Enter the unique identifier for the catchment")

def dem_map(unique_field, unique_value):
    
    # Setting the path to the work environment
    dir_abs = Path(Path().resolve() / f'{unique_field}_{unique_value}')

    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')

    # Open Data Cube Product
    dc_data = config.get('dc_product', 'DEM')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    # Read shp file to geopandas dataframe
    grid_path = Path(dir_abs / 
                     f'active_data/{unique_field}_{unique_value}_grid_WGS84.shp')
    grid = gpd.read_file(grid_path)

    xmin, ymin, xmax, ymax = grid.total_bounds
    
    # Getting data from ODC
    array_odc = read_odc(dc_data, xmin, ymin, xmax, ymax )
    
    # DEM values for SHETRAN mask
    grid_DEM = dem_catchm_values(dir_abs, array_odc, grid_path, grid, ND)
    
    # Creating frinal map
    shetran_dem(grid_DEM, crs_local, dir_abs, unique_field, unique_value)
    
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    dem_map()
    