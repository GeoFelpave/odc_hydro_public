# -*- coding: utf-8 -*-
"""==============================================================================
Title          :setting_model_inputs_neuro.py
Description    :
Author         :Luis FP Velasquez
Date           :2024-01
Version        :1.0
Usage          :python setting_model_inputs_neuro.py --ctchmnt_path=<path> --unique_field=<column id>
Notes          :
python version :
=============================================================================="""
# =============================================================================
# Packages - Libraries
# =============================================================================
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
import subprocess
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

def catchm_ids(path_shp:str, 
               uid_field:str
               ):
    
    # Read shp file to geopandas dataframe
    catchm = gpd.read_file(Path(path_shp))
    
    # Return list of unique catchment ids
    uid_lst = catchm[uid_field].to_list()
    
    return uid_lst

def create_catchm_shp(path_shp,
                      active_path, 
                      uid_field,
                      uid_value,
                      crs_gbl
                      ):
    
    # Read shp file to geopandas dataframe
    catchm = gpd.read_file(Path(path_shp))
    
    # query geopandas with all data
    if isinstance(uid_value, str):
        exp = f'{uid_field} == "{uid_value}"'
    else:
        exp = f'{uid_field} == {uid_value}'
        
    final_gpd = catchm.query(exp)
    
    # create individual files
    final_gpd.to_file(Path(active_path / f'{uid_field}_{uid_value}.shp'))
    
    # create individua files in wgs84
    # Change geodataframe to WGS84
    _shp_wgs84 = final_gpd.to_crs(crs_gbl)
    _shp_wgs84.to_file(Path(active_path / f'{uid_field}_{uid_value}_WGS84.shp'))
    

def folder_struct(lst_uid:list, 
                  wrk_path:str,
                  path_shp:str, 
                  uid_field:str,
                  crs_gbl:str,
                  code_iso:str,
                  ):
    
    # create top layer
    uid_path = Path(wrk_path / f'ODC_{code_iso}')
    Path(uid_path).mkdir(parents=True, exist_ok=True)
    
    # create subfolders
    Path(Path(uid_path / 'active_data')).mkdir(parents=True, exist_ok=True)
    Path(Path(uid_path / 'timeseries')).mkdir(parents=True, exist_ok=True)
    Path(Path(uid_path / 'temp_data')).mkdir(parents=True, exist_ok=True)
    
    # Go through uid creating the folders needed
    for uid in lst_uid:
                
        # save individual catchments
        create_catchm_shp(path_shp, Path(uid_path / 'active_data'),uid_field,uid,crs_gbl)

    

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating model input files\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating input files")
@click.option("--ctchmnt_path", 
              default=None, 
              help="Enter the complete path to the catchment shapefile")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")

def shetran_inputs(ctchmnt_path,
                   unique_field):
    
    # Setting the path to the work environment
    dir_abs = Path().resolve()
    
    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')
    
    # get iso code
    iso_code = config.get('iso_code', 'country')
    
    
    ################################
    # PROCESS
    ################################

    # Get uid
    uid_lst = catchm_ids(ctchmnt_path, unique_field)
    
    # # Set up folder structure
    log.info("\n --- \n Setting folder structure\n --- \n")
    folder_struct(uid_lst, dir_abs, 
                  ctchmnt_path, unique_field, 
                  crs_global, iso_code)
    
    # # ################################
    # # # CREATING INPUTS
    # # ################################
    

    for uid in uid_lst:
        
        log.info("\n------------------------------------- \n")
        log.info(f"------Working on catchment = {uid}--------")
        log.info("\n------------------------------------- \n")
        
        # Creating dem input file
        subprocess.run(['python','setting_dem_odc_neuro.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
        # Creating Land cover input file
        subprocess.run(['python','setting_land_cover_odc_neuro.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
        # Creating subsurface and project input file
        subprocess.run(['python','setting_soil_properties_odc_neuro.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
    #     # # Creating Lake input file
    #     # subprocess.run(['python','setting_lakes_odc.py',
    #     #                 f'--unique_field={unique_field}',
    #     #                 f'--unique_value={uid}'])
        
    #     # # Creating metereological input file
    #     # subprocess.run(['python','setting_climate_variables_odc_neuro.py',
    #     #                 f'--unique_field={unique_field}',
    #     #                 f'--unique_value={uid}'])
        
    #     break
        
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    shetran_inputs()
