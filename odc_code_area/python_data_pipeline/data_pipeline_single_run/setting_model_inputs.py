# -*- coding: utf-8 -*-
"""==============================================================================
Title          :00_setting_model_inputs.py
Description    :Create the input files for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-09
Version        :1.0
Usage          :python setting_model_inputs.py 
                --ctchmnt_path='/home/geofelpave/Documents/1_PhD_SharePoint/
                LFPV - PhD - Documents/00_PhD_main/006_GitHub/0065_odc_hydro/odc_code_area/
                python_data_pipeline/data_pipeline_single_run/catchments/Test catchments/UK_irina_luis_analysis.shp'
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
                      uid_value
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

def folder_struct(lst_uid:list, 
                  wrk_path:str,
                  path_shp:str, 
                  uid_field:str
                  ):
    
    # Go through uid creating the folders needed
    for uid in lst_uid:
        # create top folder
        uid_path = Path(wrk_path / f'{uid_field}_{uid}')
        Path(uid_path).mkdir(parents=True, exist_ok=True)
        
        # create subfolders
        Path(Path(uid_path / 'active_data')).mkdir(parents=True, exist_ok=True)
        Path(Path(uid_path / 'model_input')).mkdir(parents=True, exist_ok=True)
        Path(Path(uid_path / 'temp_data')).mkdir(parents=True, exist_ok=True)
        
        # save individual catchments
        create_catchm_shp(path_shp, Path(uid_path / 'active_data'),uid_field,uid)

    

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
    
    ################################
    # PROCESS
    ################################

    # Get uid
    uid_lst = catchm_ids(ctchmnt_path, unique_field)
    
    # # Set up folder structure
    log.info("\n --- \n Setting folder structure\n --- \n")
    folder_struct(uid_lst, dir_abs, ctchmnt_path, unique_field)
    
    # ################################
    # # CREATING INPUTS
    # ################################

    for uid in uid_lst:
        
        log.info("\n------------------------------------- \n")
        log.info(f"------Working on ID = {uid}-----------")
        log.info("\n------------------------------------- \n")
        
        # Creating mask input file
        basin_path = Path(dir_abs / f'{unique_field}_{uid}/active_data/{unique_field}_{uid}.shp')
        subprocess.run(['python','setting_mask_odc.py',
                        f'--ctchmnt_path={basin_path}',
                        f'--unique_field={unique_field}'])
    
        # Creating dem input file
        subprocess.run(['python','setting_dem_odc.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
        # Creating Land cover input file
        subprocess.run(['python','setting_land_cover_odc.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
        # Creating Lake input file
        subprocess.run(['python','setting_lakes_odc.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
        # Creating metereological input file
        subprocess.run(['python','setting_climate_variables_odc.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
        
        # Creating subsurface and project input file
        subprocess.run(['python','setting_soil_properties_odc.py',
                        f'--unique_field={unique_field}',
                        f'--unique_value={uid}'])
                
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    shetran_inputs()
