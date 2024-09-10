# -*- coding: utf-8 -*-
"""==============================================================================
Title          :setting_climate_variables_odc_neuro.py
Description    :Create the precipitation, evaporation, temperature data
Author         :Luis FP Velasquez
Date           :2024-01
Version        :1.0
Usage          :python setting_climate_variables_odc_neuro.py
Notes          :
python version :
=============================================================================="""
# =============================================================================
# Packages - Libraries
# =============================================================================
import os, sys
from pathlib import Path
import xarray as xr
import pandas as pd
import geopandas as gpd
import rioxarray
import rasterstats
from datetime import datetime
import time
import logging
import click
import datetime
from tqdm.dask import TqdmCallback
import warnings
import rasterio
warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)

# config_handler.set_global(length=40, bar='smooth', spinner='fish2')


# libraries for working with year and months
import calendar # use to check for leap years
from calendar import monthrange

# ODC modules
import datacube
import numpy as np

# Custom modules
import utils

# Modules for config.ini
from configparser import ConfigParser
config = ConfigParser()


# Multiprocessing
import multiprocessing
from multiprocessing import get_start_method
from multiprocessing import Pool

# Set multiprocessing context
if get_start_method() is not None:
    print(f'Context is already set to: {get_start_method()}')
else:
    multiprocessing.set_start_method('fork') # this is necessary for ipython
    print(f'Context set to: {get_start_method()}')

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

def grid_latlon(soil_grid_df, crs_num):
    """_summary_

    Args:
        soil_grid_df (_type_): _description_
        crs_num (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Get the centroid for each grid
    temp_df = pd.DataFrame()
    temp_df = soil_grid_df.to_crs(crs_num)
    temp_df['centroid'] = temp_df['geometry'].centroid

    # Create lat (Y) and lon (X) columns
    temp_df['lat'] = temp_df['centroid'].y.astype(int)
    temp_df['lon'] = temp_df['centroid'].x.astype(int)

    # Return the columns needed
    df_grid = pd.DataFrame(temp_df[['lat', 'lon', 'geometry']].copy())

    
    return df_grid

def mask_uid(wrk_path, crs_lcl,
             crs_gbl, uid_field,
             uid_value):
    """_summary_

    Args:
        wrk_path (_type_): _description_
        crs_lcl (_type_): _description_
        crs_gbl (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    # Read shp file to geopandas dataframe
    shp_path = Path(wrk_path / 
                     f'active_data/{uid_field}_{uid_value}_WGS84.shp')
    _shp = gpd.read_file(shp_path)

    xmin, ymin, xmax, ymax = _shp.total_bounds
    
    return xmin, ymin, xmax, ymax, _shp

def read_odc(start_date, end_date,
             dc_product:str, xmn:int, 
             ymn:int, xmx:int, ymx:int):
    
    """_summary_

    Args:
        start_date (str): _description_
        end_date (str): _description_
        dc_product (str): _description_
        xmn (int): _description_
        ymn (int): _description_
        xmx (int): _description_
        ymx (int): _description_

    Returns:
        _type_: _description_
    """
    
    log.info(f'\n --- \n Reading {dc_product} from ODC\n --- \n')
    
    # Set ODC application
    dc = datacube.Datacube(app="era_five")

    # Get data
    # Load data from the datacube
    buffer = 0.125
    START_DATE = start_date
    END_DATE = end_date
    ds = dc.load(product=dc_product,
                lat=(ymn - buffer, ymx + buffer),
                lon=(xmn - buffer, xmx + buffer),
                time=(START_DATE, END_DATE),
                dask_chunks={'time': 1, 'longitude': 200, 'latitude': 200}
                )

    # Print output data
    return ds

# Define function to use in the pool
def meteo_values(climate_var, dc_array,
                 df_mask, time_split, 
                 path_mask, wrk_path,
                 nd_value):
    """_summary_

    Args:
        climate_var (_type_): _description_
        dc_array (_type_): _description_
        df_mask (_type_): _description_
        time_split (_type_): _description_
        path_mask (_type_): _description_
        wrk_path (_type_): _description_
        nd_value (_type_): _description_
    """
    

    # Create list of timestamps
    lst_time = dc_array.coords['time'].values
    
    # Create file name using the date
    date_str = np.datetime_as_string(lst_time[time_split], unit='D')
    # date_str = file_name.replace('-', "")
    
    # Read nc file as xarray and select one timestamp at the time
    # prep_array = lst_cds[count_cds].to_xarray()
    clim_tif = dc_array.isel(time=slice(time_split, time_split + 1))

    # Create tif file need for rasterstats
    _tif = Path(wrk_path / f'temp_data/{climate_var}_{date_str}.tif')
    clim_tif[climate_var].rio.to_raster(_tif)

    
    ##################################################################################
    # RUNNING THE STATS FOR THE CATCHMENT
    ##################################################################################
    # Find the mean precipitation for each grid in the catchment
    # This produces a dict following the same order as the grid geodataframe
    catchm_stats_mean = rasterstats.zonal_stats(str(path_mask), 
                                                str(_tif), stats="mean", 
                                                all_touched=True, nodata=nd_value)

    # Change the dictionary to a list to be added to the mask
    stats_value = [x['mean'] for x in catchm_stats_mean]

    # Temp dataframe to store rasterstats dataset
    catch_stats = pd.DataFrame() # need to make sure there are not errors
    catch_stats = pd.DataFrame([[date_str, stats_value[0]]], columns=['date', climate_var])    
    
    _csv = Path(wrk_path / f'temp_data/{climate_var}_{date_str}.csv')
    catch_stats.to_csv(_csv, index=False)
    
    print(f'{climate_var}_{date_str} done')

def meteo_processing(meteo_var, dc_array, df_catch,
                     wrk_path, nd_value, idx_lst,
                     uid_field, uid_value):
    """_summary_

    Args:
        meteo_var (_type_): _description_
        dc_array (_type_): _description_
        df_catch (_type_): _description_
        wrk_path (_type_): _description_
        nd_value (_type_): _description_
        idx_lst (_type_): _description_
    """
    # Make sure Temp folder is empty to avoid errors
    utils.file_remove(Path(Path(wrk_path / 'temp_data/')), 'all')
    
    # Read shp file to geopandas dataframe
    mask_path = Path(wrk_path / 
                     f'active_data/{uid_field}_{uid_value}_WGS84.shp')
    
    # Pool use async to wait for all task to finish
    results = []
    with Pool(8) as pool:
        start=time.time()
        result = pool.starmap_async(meteo_values, [(meteo_var, dc_array, 
                                              df_catch, t, mask_path, 
                                              wrk_path, nd_value) for t in idx_lst])
        # wait for task to finish
        result.wait()
        # pool.terminate()
        print("Time Taken: ", str(time.time()-start), flush=True) 


def meteo_csv(wrk_path, meteo_var,
              start_date, end_date,
              uid_field, uid_value):
    """_summary_

    Args:
        wrk_path (_type_): _description_
        meteo_var (_type_): _description_
        start_date (_type_): _description_
        end_date (_type_): _description_
    """
    # Delete tif file to avoid any error
    utils.file_remove(Path(Path(wrk_path / 'temp_data/')), 'tif')
    log.info("\n --- \n Temp folder is ready for the process")
    
    log.info(f" Creating {meteo_var} csv file ready for SHETRAN\n --- \n")
    
    # Add file to a list 
    '''this is needed as when reding the files they might no be read in the right order'''
    lst_files = []

    for path in Path(Path(wrk_path / f'temp_data/')).glob(f'*.csv'):
            lst_files.append(path)
    
    
    # Make sure the list is order by file name
    lst_files.sort()

    # Using the list create dataframe
    df_prep = pd.concat((pd.read_csv(f) for f in lst_files), ignore_index=True)

    # # Make sure we only using the right columns
    # df_prep = df_prep[['uid',f'{meteo_var}']].copy()
    # df_prep.columns = ['uid',meteo_var]
    
    # Make sure to have positive values for evapotranspiration      
    if 'evaporation' in meteo_var:
        df_prep[meteo_var] = df_prep[meteo_var]\
            .map(lambda a: float(a) * -1 if float(a) < 0 else float(a))
    
    # Make sure to have positive values for evapotranspiration      
    if 'evapotranspiration' in meteo_var:
        df_prep[meteo_var] = df_prep[meteo_var]\
            .map(lambda a: float(a) * -1 if float(a) < 0 else float(a))
    
    # Pivot the table using the uid
    # df_prep = df_prep.pivot(columns='uid')[meteo_var]
    df_prep_final = df_prep.copy()

    # delete all csv file before creating the final csv for climate variable
    utils.file_remove(Path(Path(wrk_path / 'temp_data/')), 'csv')
    
    # Saving file to csv
    file_name = f'{uid_field}_{uid_value}_{meteo_var}_{start_date.split("-")[0]}'\
        f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
        
    _csvFinal = Path(wrk_path / f'timeseries/{file_name}.csv')
    df_prep_final.to_csv(_csvFinal, index=False)
        
    # log.info(f" {uid_field}_{uid_value} {meteo_var} file created")

def obs_disch(path_csv, uid_value,
              start_date, end_date):
    
    # Put subdirectories in a list 
    lst_subdir = os.listdir(path_csv)
    
    # Work with the loop
    for _subdir in lst_subdir:
        
        # only use the subdir catchmet id
        flag = _subdir.split('_')[-1]
        
        if flag == uid_value:
            file_path = Path(f'{path_csv}/{_subdir}')
            _csv = list(file_path.iterdir())
            
            # Check if there is more than one file
            # set error check
            if len(_csv) == 1: 
                
                # Open discharge data
                # The skiprows have been set for the data from NRFA
                df_disc = pd.read_csv(_csv[0], skiprows=6)
                
                # set col names
                col_name = ['date', 'QObs(mm/d)']
                df_disc.columns = col_name
                
                # select data between the config dates
                df_qObs = df_disc.query('date >= @start_date and date <= @end_date')
                
                return df_qObs
            else:
                log.info('\n --- \n************ ERROR ************\n'\
                    ' Too many csv files in folder, please check \n' \
                        ' Only one file should exist in the catchment folder\n'\
                            '*******************************\n --- \n')
                sys.exit()
    

def timeseries_csv(wrk_path, start_date,
                   end_date, uid_field, 
                   uid_value, df_dischargeObs):
    
    log.info(f" Creating timeseries csv file ready for Neurohydrology\n --- \n")
    
    # Add file to a list 
    '''this is needed as when reding the files they might no be read in the right order'''
    lst_files = []

    for path in Path(Path(wrk_path / f'timeseries/')).glob(f'{uid_field}_{uid_value}*.csv'):
            lst_files.append(path)
    
    # ID_STRING_94001_total_precipitation_201101_201112.csv

    # Make sure the list is order by file name
    lst_files.sort()

    # Using the list create dataframe
    df_timeseries = pd.concat((pd.read_csv(f) for f in lst_files), ignore_index=False, axis=1)
    
    # remove the duplicate columns
    df_timeseries = df_timeseries.loc[:,~df_timeseries.columns.duplicated()].copy()
    
    # add the observed discharge dataset
    df_final = df_timeseries.merge(df_dischargeObs, on='date', how='inner')

    # Saving file to csv
    file_name = f'CAMELS_GB_hydromet_timeseries_{uid_value}_'\
        f'{start_date.split("-")[0]}'\
        f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
        
    _csvFinal = Path(wrk_path / f'timeseries/{file_name}.csv')
    df_final.to_csv(_csvFinal, index=False)
        
    # Remove unnecessary files
    # files in timeseries
    for f in lst_files:
        f.unlink()
        
#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating SHETRAN meteorological files\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating meteorological data map and csv files")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")
@click.option("--unique_value", 
              default=None, 
              help="Enter the unique identifier for the catchment")
@click.option("--csv_path", 
              default=None, 
              help="Enter the path to the folder with observational discharge data")

def meteo(unique_field,
          unique_value,
          csv_path):
    
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

    
    # Get the dates for the file name
    date_start = config.get("time_period", "start_date")
    date_end = config.get("time_period", "end_date")
    

    # Set data cube product to read
    ''' [era5_reanalysis_tp_daily, era5_reanalysis_pev_daily, era5_reanalysis_tmean_daily,
    era5_reanalysis_tmax_daily,era5_reanalysis_tmin_daily]'''
    DC_PRODUCT = config.get('dc_product', 'lst_era_products').split(',')

    # Set era5 var
    ''' [total_precipitation, potential_evaporation, 2m_temperature_mean,
    2m_temperature_max,2m_temperature_min]'''
    ERA_VAR = config.get('era_land', 'lst_era_var').split(',')
    
    # Adding UID to catchment mask
    '''Python also comes with an unpacking operator, 
    which is denoted by *. Say that we only cared 
    about the first item returned. We still need 
    to assign the remaining values to another variable, 
    but we can easily group them into a single variable, 
    using the unpacking operator. The last element denoted
    by * are returned packed into a list'''
    xmin, ymin, xmax, ymax, *df_catchm = mask_uid(dir_abs, crs_local, 
                                                    crs_global, unique_field,
                                                    unique_value)
    
    # # Working with each meteorological dataset
    arr_lst = []
        
    for count, value in enumerate(ERA_VAR):
                
        # Getting data from ODC        
        array_odc = read_odc(date_start, date_end,
                            DC_PRODUCT[count], xmin, 
                            ymin, xmax, ymax)
        
        # Processing meteo data in a parallel process
        # Create list of timestamps to bconfig.get('time_period', 'end_date')e used for the pool
        log.info(f"\n --- \n Processing {value} data in a pool process")
        lst_time = array_odc.coords['time'].values # list of date values
        lst_index = list(range(len(lst_time))) # list index values
        
        meteo_processing(value, array_odc,
                         df_catchm[0], dir_abs,
                         ND, lst_index, unique_field,
                         unique_value)
        
        # Create csv file ready for SHETRAN
        meteo_csv(dir_abs, value,
                  date_start, date_end, 
                  unique_field, unique_value)
    
    # # Create timeseries file ready for neurohydrology
    # utils.file_remove(Path(Path(dir_abs / 'temp_data/')), 'csv')
    
    # Get obs discharge for each value as dataframe
    qObs_data = obs_disch(csv_path, unique_value,
                          date_start, date_end)
    
    timeseries_csv(dir_abs, date_start, 
                   date_end, unique_field, 
                   unique_value, qObs_data)

    
    log.info(f"\n --- \n Files for {unique_field}_{unique_value} {ERA_VAR} have been created \n"\
             f"between the periods of {date_start} and {date_end} \n --- \n ")
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    meteo()

    
    