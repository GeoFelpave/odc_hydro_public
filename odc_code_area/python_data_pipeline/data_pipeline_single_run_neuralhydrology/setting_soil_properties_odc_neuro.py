# -*- coding: utf-8 -*-
"""==============================================================================
Title          :06_setting_soil_properties_odc.py
Description    :Create the surface data needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 06_setting_soil_properties_odc.py
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
import numpy as np
import logging
import click

# XML libraries
from dicttoxml import dicttoxml
import xml.etree.ElementTree as ET
# from dict2xml import dict2xml
from xml.dom.minidom import parseString

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


def read_odc(dc_product:str, xmn:int, 
             ymn:int, xmx:int, ymx:int):
    
    
    log.info(f'\n --- \n Reading {dc_product} from ODC\n --- \n')
    
    # Set ODC application
    dc = datacube.Datacube(app="hihydrosoil")

    # Get data
    # Load data from the datacube
    buffer = 0.125
    ds = dc.load(product=dc_product,
                lat=(ymn - buffer, ymx + buffer),
                lon=(xmn - buffer, xmx + buffer),
                )

    # Print output data
    return ds

def soil_catchm_values(file_name, dc_array, 
                       df_mask, soil_var, 
                       wrk_path, soil_depth, 
                       nd_value, uid_field, uid_value):
    
    # Read shp file to geopandas dataframe
    mask_path = Path(wrk_path / 
                     f'active_data/{uid_field}_{uid_value}_WGS84.shp')
    
    # set dataframe
    temp_stats_df = pd.DataFrame()
    
    #create temp .tif file
    _tif = Path(wrk_path / f'temp_data/{file_name}.tif')
    dc_array[soil_var].rio.to_raster(_tif)
    
    # Find the mean elevation for each grid in the catchment
    # This produces a dict following the same order as the grid geodataframe
    # soilt_stats = rasterstats.zonal_stats(str(path_grid), str(_tif), categorical=True, nodata=-999)
    soil_stats =  rasterstats.zonal_stats(str(mask_path), str(_tif), 
                                          stats="mean", all_touched=True, nodata=nd_value)

    # Create a list containing the dictionary key (land cover type) with the Max value for each grid in the catchment
    prop_mean = [x['mean'] for x in soil_stats]

    # Add the largest lc to grid geodataframe
    df_mask['property_value'] = prop_mean
    
    if soil_depth != '200_2000cm':
        # change variables to expected number by multiplying by 0.0001
        df_mask['property_value'] = df_mask['property_value'] * 0.0001
    else:
        # df_mask['property_value'] = (10**(df_mask['property_value'] / 100))* 100000000000
        df_mask['property_value'] = (10**(df_mask['property_value'] / 100))* 10000000

    # add the soil depth for reference
    df_mask['soil_depth'] = soil_depth
    
    # add the soil property name for reference
    if soil_var == 'logKferr':
        soil_var = 'Ksat'
        
    df_mask['property_name'] = soil_var
    
    # Change geopandas to pandas ready to create csv file for the mean and min elevation
    df_properties_pd = pd.DataFrame()
    df_properties_pd = pd.DataFrame(df_mask[[uid_field, 'property_value', 'soil_depth', 'property_name']].copy())
    
    # delte temp .tif file
    utils.file_remove(Path(wrk_path / 'temp_data/'), 'tif')

    return df_properties_pd

# def soil_catchm_dataframe(soil_grid_df, crs_lcl, soil_var, soil_depth):
#     # Get the centroid for each grid
#     temp_df = pd.DataFrame()
#     temp_df = soil_grid_df.copy()
#     temp_df = soil_grid_df.to_crs(crs_lcl)
#     # temp_df['centroid'] = temp_df['geometry'].centroid

#     # # Create lat (Y) and lon (X) columns
#     # temp_df['lat'] = temp_df['centroid'].y.astype(int)
#     # temp_df['lon'] = temp_df['centroid'].x.astype(int)

#     # Change geopandas to pandas ready to create csv file for the mean and min elevation
#     # df_grid_lc = pd.DataFrame(temp_df[['lat', 'lon', 'property_value']].copy())
#     # df_grid_lc = pd.DataFrame(temp_df[['lat', 'lon', 'SHETRAN_ID', 'property_value']].copy())
#     df_grid_lc = pd.DataFrame(temp_df[['lat', 'lon', 'property_value']].copy())
    
#     # add the soil depth for reference
#     df_grid_lc['soil_depth'] = soil_depth
    
#     # add the soil property name for reference
#     df_grid_lc['property_name'] = soil_var
    
#     return df_grid_lc

def soil_properties(product_name,
                    xmn, ymn, xmx, ymx,
                    df_mask, wrk_path, 
                    nd_value, crs_lcl,
                    uid_field, uid_value):
    
    var_name = product_name.split("_")[1]
    depth = '_'.join(product_name.split("_")[2:])

    # ODC data
    array_odc = read_odc(product_name, xmn, ymn, xmx, ymx)

    
    # Hihydrosoil values for SHETRAN mask
    df_properties = soil_catchm_values(product_name, array_odc, 
                                       df_mask, var_name,
                                       wrk_path, depth, 
                                       nd_value,uid_field, uid_value) 
    
    
    
    # df_soil = soil_catchm_dataframe(df_properties, crs_lcl,
    #                                 var_name, depth)
    return df_properties
    
def geology(product_name, xmn, ymn, 
            xmx, ymx,
            df_mask, wrk_path, 
            nd_value, crs_lcl, 
            lst_dataframes,
            uid_field, uid_value):
    
    df_properties = pd.DataFrame()
    df_soil = pd.DataFrame()

    # Set variables
    value = product_name
    
    df_geology = soil_properties(value, xmn, 
                                 ymn, xmx, ymx,
                                 df_mask, wrk_path, 
                                 nd_value, crs_lcl,
                                 uid_field, uid_value)

    # Add data to the list of dataframes
    lst_dataframes.append(df_geology)

    # Create dataframe for other properties at the same depth
    prop_lst = ['ALFA','N','WCres','WCsat']
    prop_values = [0.01,5.0,0.2,0.3]

    for count, value in enumerate(prop_lst):
        # copy the Ksat dataframe to create the other
        # this only works because the soil properties are constant across grids
        df_temp = pd.DataFrame() # create clear df for each loop iteration
        df_temp = df_geology.copy()
        
        # Replace values in dataframe
        df_temp['property_value'] = prop_values[count]
        df_temp['property_name'] = value
        
        # add data to the list of dataframes 
        lst_dataframes.append(df_temp)
    
    log.info("\n --- \n Geology layer has been created\n --- \n")
    return lst_dataframes


#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("\n --- \n Creating soil attributes files\n --- \n")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")
@click.option("--unique_value", 
              default=None, 
              help="Enter the unique identifier for the catchment")

def surface(unique_field, unique_value):

    # # Make sure temp folder directory is empty
    # utils.file_remove(Path(dir_abs / 'temp_data/'), 'all')

    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')
    
    # get iso code
    iso_code = config.get('iso_code', 'country')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    # Setting the path to the work environment
    dir_abs = Path(Path().resolve() / f'ODC_{iso_code}')
    
    # Get the dates for the file name
    date_start = config.get('time_period', 'start_date')
    date_end = config.get('time_period', 'end_date')

    # Open Data Cube Product
    lst_dc_data = config.get('dc_product', 'lst_soil_products').split(',')
    
    # Read shp file to geopandas dataframe
    shp_path = Path(dir_abs / 
                     f'active_data/{unique_field}_{unique_value}_WGS84.shp')
    _shp = gpd.read_file(shp_path)

    xmin, ymin, xmax, ymax = _shp.total_bounds
        
    # Getting data from ODC and process the soil data from Hihydrosoil
    # Create column of values using hihydrosil data
    
    # list to stor dataframe
    lst_df = []
    
    for count, value in enumerate(lst_dc_data):
        
        df_soil = soil_properties(value, xmin, 
                                  ymin, xmax, ymax,
                                  _shp, dir_abs,
                                  ND, crs_local,
                                  unique_field, unique_value)
        # Append dataframe to list
        lst_df.append(df_soil)
            
    # Set geology
    # Adding geology dataframe to the list of dataframes
    geo_product = 'ghlymps_logKferr_200_2000cm'
    lst_soil_column = geology(geo_product, xmin, 
                           ymin, xmax, ymax,
                           _shp, dir_abs,
                           ND, crs_local, lst_df,
                           unique_field, unique_value)
    
    # Create csv file with all the properties of the soild column
    df_soil_all = pd.concat(lst_soil_column)
    
    # create column for new col names
    df_soil_all['new_col_name'] = df_soil_all[
        ['property_name','soil_depth']].apply(lambda row: '_'
                                              .join(row.values.astype(str)), 
                                              axis=1)
    
    # only work with the right
    df_soil_all = df_soil_all[[unique_field,
                               'property_value',
                               'new_col_name']]
    
    df_final = df_soil_all.pivot(index=unique_field, 
                                 columns='new_col_name', 
                                 values='property_value').reset_index()
    
    # this is necessary to match neuralhydrology requirements for CAMELS GB
    df_final.rename(columns={unique_field:'gauge_id'}, inplace=True)
    
    # save to csv
    df_final.to_csv(Path(dir_abs /
                        f'temp_data/ODC_{iso_code}_'\
                        'soil_attributes_'\
                        f'{unique_value}.csv'), 
                    index=False)
    
    log.info(f'\n --- \n Soil attributes have been created\n --- \n')
    
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    surface()
    