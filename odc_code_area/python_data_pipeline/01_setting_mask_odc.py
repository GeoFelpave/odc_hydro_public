# -*- coding: utf-8 -*-
"""==============================================================================
Title          :01_setting_mask_odc.py
Description    :Create the catchment mask needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 01_setting_mask_odc.py --ctchmnt_path=<path>
Notes          :
python version :
=============================================================================="""
# =============================================================================
# Packages - Libraries
# =============================================================================
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
import numpy as np
import numbers
from pathlib import Path
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

def fishnet(path_shp:str,crs_gbl:int,crs_lcl:int,
            wdth:int, hght:int):
    
    # Read shp file to geopandas dataframe
    catchm = gpd.read_file(Path(path_shp))
    
    # If catchment in WGS need to project before doing anything - EPSG:3116
    # this code only transfroms from EPSG 4326 - for other CRS the code will need to be modified
    # catchm.crs.to_epsg(20) see https://gis.stackexchange.com/questions/326690/explaining-pyproj-to-epsg-min-confidence-parameter
    if catchm.crs.to_epsg(20) == crs_gbl:
        # user_crs = crs_lcl
        catchm = catchm.to_crs(f'{crs_lcl}')
        
        # Get the bbox for the polygon
        xmin,ymin,xmax,ymax =  catchm.total_bounds

        # Create the fishnet
        width = wdth
        height = hght
        rows = int(np.ceil((ymax-ymin) /  height))
        cols = int(np.ceil((xmax-xmin) / width))
        XleftOrigin = xmin
        XrightOrigin = xmin + width
        YtopOrigin = ymax
        YbottomOrigin = ymax- height
        polygons = []
        for i in range(cols):
            Ytop = YtopOrigin
            Ybottom =YbottomOrigin
            for j in range(rows):
                polygons.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
                Ytop = Ytop - height
                Ybottom = Ybottom - height
            XleftOrigin = XleftOrigin + width
            XrightOrigin = XrightOrigin + width
    else:
        user_crs = catchm.crs.to_epsg(20)
        # Get the bbox for the polygon
        xmin,ymin,xmax,ymax =  catchm.total_bounds

        # Create the fishnet
        width = wdth
        height = hght
        rows = int(np.ceil((ymax-ymin) /  height))
        cols = int(np.ceil((xmax-xmin) / width))
        XleftOrigin = xmin
        XrightOrigin = xmin + width
        YtopOrigin = ymax
        YbottomOrigin = ymax- height
        polygons = []
        for i in range(cols):
            Ytop = YtopOrigin
            Ybottom =YbottomOrigin
            for j in range(rows):
                polygons.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
                Ytop = Ytop - height
                Ybottom = Ybottom - height
            XleftOrigin = XleftOrigin + width
            XrightOrigin = XrightOrigin + width
        
    log.info(f'\n --- \n The grid mask has been generated as a list of values\n --- \n')
    return polygons

def shetran_fishnet(py_fishnet:str,crs_gbl:int, crs_lcl:int,
                    path_shp:str,nd_value:int, wrk_path:str):
    
    # Create grid as geopandas dataframe
    grid = gpd.GeoDataFrame({'geometry':py_fishnet})
    grid.set_crs(epsg=f'{crs_lcl}', inplace=True)
    
    # Read shp file to geopandas dataframe
    catchm = gpd.read_file(Path(path_shp))
    catchm = catchm.to_crs(f'{crs_lcl}')  #make sure is the projected CRS

    # Intersect grid with catchment to add SHETRAN ID
    # catchm is in the projected crs from the section above
    intersection = grid['geometry'].intersects(catchm.loc[0,'geometry'])

    # Add the intersection to the grid geodataframe
    grid['intersect'] = intersection

    # Assing SHETRAN ID if value is True
    grid["SHETRAN_ID"] = np.where(grid['intersect'] == True, 0, nd_value)

    # Save grid file as shp file - this will be an input on other scripts
    # add CRS before saving
    grid_wgs84 = grid.copy()
    # Change CRS to EPSG: 4326
    grid_wgs84 = grid.to_crs(crs_gbl)

    grid_wgs84.to_file(Path(wrk_path / 'shetran_data/active_data/final_mask_wgs84.shp'))


    # Get the centroid for each grid
    grid['centroid'] = grid['geometry'].centroid

    # Create lat (Y) and lon (X) columns
    grid['lat'] = grid['centroid'].y.astype(int)
    grid['lon'] = grid['centroid'].x.astype(int)

    # Change geopandas to pandas ready to create csv file
    df_grid = pd.DataFrame(grid[['SHETRAN_ID', 'lat', 'lon']].copy())
    df_grid

    # Pivoting dataframe to replicate SHETRAN format
    # Pivoting dataframe using lon as column and lat as row
    df_pivot = df_grid.pivot(index='lat', columns='lon', values='SHETRAN_ID')
    df_pivot = df_pivot.sort_index(ascending=False)
    
    # Creating the text file for min elevation
    utils.shetran_csv_file(wrk_path, 'final_catchment_mask_SHETRAN', df_pivot, 'd')
    
    log.info(f"\n --- \n The grid mask has been created check: \n"\
             f"{Path(wrk_path / 'shetran_data/active_data/final_mask_wgs84.shp')}\n --- \n")
    return grid
    

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating catchment mask\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating catchment mask")
@click.option("--ctchmnt_path", 
              default=None, help="Enter the complete path to the catchment shapefile")

def masking(ctchmnt_path:str):
    
    # Setting the path to the work environment
    dir_abs = Path().resolve()

    # Read config file values
    config.read('config.ini')
    
    # grid resolution setting
    WIDTH = config.getint('res_setting', 'WIDTH')
    HEIGHT = config.getint('res_setting', 'HEIGHT')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'COL')
    
    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    fishnet_py = fishnet(ctchmnt_path,crs_global,
                         crs_local,WIDTH, HEIGHT)
    
    shetran_fishnet(fishnet_py,crs_global, crs_local,
                    ctchmnt_path,ND, dir_abs)
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    masking()
    