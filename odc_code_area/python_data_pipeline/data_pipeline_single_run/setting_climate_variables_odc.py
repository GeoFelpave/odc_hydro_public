# -*- coding: utf-8 -*-
"""==============================================================================
Title          :05_setting_climate_variables_odc.py
Description    :Create the precipitation, evaporation, temperature data needed for SHETRAN
Author         :Luis FP Velasquez
Date           :2023-07
Version        :1.0
Usage          :python 05_setting_climate_variables_odc.py
Notes          :
python version :
=============================================================================="""
# =============================================================================
# Packages - Libraries
# =============================================================================
import os
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

def justify(a, invalid_val=0, axis=1, side='left'):    
    """
    Justifies a 2D array

    Parameters
    ----------
    A : ndarray
        Input array to be justified
    axis : int
        Axis along which justification is to be made
    side : str
        Direction of justification. It could be 'left', 'right', 'up', 'down'
        It should be 'left' or 'right' for axis=1 and 'up' or 'down' for axis=0.

    """

    if invalid_val is np.nan:
        #change to notnull
        mask = pd.notnull(a)
    else:
        mask = a!=invalid_val
    justified_mask = np.sort(mask,axis=axis)
    if (side=='up') | (side=='left'):
        justified_mask = np.flip(justified_mask,axis=axis)
    #change dtype to object
    out = np.full(a.shape, invalid_val, dtype=object)  
    if axis==1:
        out[justified_mask] = a[mask]
    else:
        out.T[justified_mask.T] = a.T[mask.T]
    return out

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
        
    # # # Loop inside each nc file using the timestamp
    # # for count, value in enumerate(lst_time):
        
    # Temp dataframe to store rasterstats dataset
    grid_stats = pd.DataFrame() # need to make sure there are not errors
    grid_stats = df_mask[['uid']].copy()
    
    # Create file name using the date
    file_name = np.datetime_as_string(lst_time[time_split], unit='D')
    file_name = file_name.replace('-', "")
    
    # Read nc file as xarray and select one timestamp at the time
    # prep_array = lst_cds[count_cds].to_xarray()
    clim_tif = dc_array.isel(time=slice(time_split, time_split + 1))

    # Create tif file need for rasterstats
    _tif = Path(wrk_path / f'temp_data/{climate_var}_{file_name}.tif')
    clim_tif[climate_var].rio.to_raster(_tif)
    # prep_tif[VAR_NAME].rio.to_raster(f'Data/temp/{ERA_VAR}_{file_name}.tif')
    
    # Store file name
    # lst_fileNames.append(f'{climate_var}_{file_name}.tif')
    
    ##################################################################################
    # RUNNING THE STATS FOR EACH GRID IN THE CATCHMENT
    ##################################################################################
    # Find the mean precipitation for each grid in the catchment
    # This produces a dict following the same order as the grid geodataframe
    catchm_stats_mean = rasterstats.zonal_stats(str(path_mask), 
                                                str(_tif), stats="mean", 
                                                all_touched=True, nodata=nd_value)

    # Change the dictionary to a list to be added to the mask
    stats_value = [x['mean'] for x in catchm_stats_mean]

    # Add the climate information to grid geodataframe
    grid_stats[climate_var] = stats_value
    
    # Save to csv file ready for next stage
    # Path(Path(wrk_path / f'temp_data/{climate_var}_csv')).mkdir(parents=True, exist_ok=True)
    # _csv = Path(wrk_path / f'temp_data/{climate_var}_csv/{climate_var}_{file_name}.csv')
    # grid_stats.to_csv(_csv, index=False)
    
    _csv = Path(wrk_path / f'temp_data/{climate_var}_{file_name}.csv')
    grid_stats.to_csv(_csv, index=False)
    
    print(f'{climate_var}_{file_name} done')

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
    grid_path = Path(wrk_path / 
                     f'active_data/{uid_field}_{uid_value}_grid_WGS84.shp')
    grid = gpd.read_file(grid_path)

    # Create lat and lot values using the function
    '''3116 is the EPSG in Colombia'''
    grid_centroids = grid_latlon(grid, crs_lcl)

    # Get list of unique values of lon
    '''This will be use to create the right order for the grid id'''
    lst_lon = grid_centroids.lon.values
    lst_lon = np.unique(lst_lon, axis=0)

    # Create list to store the dataframes
    lst_df = []

    # Loop through list of longitude to create the unique values
    for count, value in enumerate(lst_lon):
        temp = grid_centroids.copy()
        temp = temp.loc[grid_centroids['lon'] == value]

        # create list with uid values
        # set range values
        start = count
        end = len(temp.lat.to_list())*len(set(grid_centroids.lon.to_list()))
        increment = len(set(grid_centroids.lon.to_list()))
        
        # create list of numbers'
        '''this uses the dimensions of the grid'''
        lst_number = np.arange(start, end, increment).tolist()
        
        # # # Add prefix to numbers to make the final uid
        # lst_uid = list(map(lambda x: f'grid_{str(x).zfill(5)}', lst_number))
        
        # Add uid as a column in dataframe
        temp['uid'] = lst_number
        
        # add temp to lst_df 
        lst_df.append(temp)


    # combine all dataframe
    grid_final = pd.concat(lst_df)
    grid_final = grid_final.reset_index()

    # create geodataframe
    grid_final = gpd.GeoDataFrame(grid_final)

    # Change geodataframe to WGS84
    grid_wgs84 = grid_final.to_crs(crs_gbl)
    
    # Increase value in UID for config file
    grid_wgs84['uid'] = grid_wgs84['uid'] + 1

    # # Get bbox
    xmin, ymin, xmax, ymax = grid_wgs84.total_bounds

    return xmin, ymin, xmax, ymax, grid_wgs84

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

def meteo_processing(meteo_var, dc_array, df_uid,
                     wrk_path, nd_value, idx_lst,
                     uid_field, uid_value):
    """_summary_

    Args:
        meteo_var (_type_): _description_
        dc_array (_type_): _description_
        df_uid (_type_): _description_
        wrk_path (_type_): _description_
        nd_value (_type_): _description_
        idx_lst (_type_): _description_
    """
    # Make sure Temp folder is empty to avoid errors
    utils.file_remove(Path(Path(wrk_path / 'temp_data/')), 'all')
    
    # Read shp file to geopandas dataframe
    mask_path = Path(wrk_path / 
                     f'active_data/{uid_field}_{uid_value}_grid_WGS84.shp')

    # Pool use async to wait for all task to finish
    results = []
    with Pool(8) as pool:
        start=time.time()
        result = pool.starmap_async(meteo_values, [(meteo_var, dc_array, 
                                              df_uid, t, mask_path, 
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

    # for path in Path(Path(wrk_path / f'temp_data/')).glob(f'*.csv'):
    #         lst_files.append(path)
    
    # Work with each csv file individually to avoid memory issues
    _df_lst = []
    _lst_csv = list(Path(wrk_path / f'temp_data/').iterdir())
    _lst_csv.sort()

    
    for csv in _lst_csv:
        log.info(f" Working with {csv.name}")
        
        df_prep = pd.read_csv(csv)
        df_prep.reset_index(drop=True, inplace=True)
        
        # Make sure we only using the right columns
        df_prep = df_prep[['uid',f'{meteo_var}']].copy()
        # df_prep.columns = ['uid',meteo_var]
        
        # Make sure to have positive values for evapotranspiration      
        if 'evaporation' in meteo_var:
            df_prep[meteo_var] = df_prep[meteo_var]\
                .map(lambda a: float(a) * -1 if float(a) < 0 else float(a))
        
        # Make sure to have positive values for evapotranspiration      
        if 'evapotranspiration' in meteo_var:
            df_prep[meteo_var] = df_prep[meteo_var]\
                .map(lambda a: float(a) * -1 if float(a) < 0 else float(a))
        
        # Make sure to have positive values for precipitation     
        if 'precipitation' in meteo_var:
            df_prep[meteo_var] = df_prep[meteo_var]\
                .map(lambda a: 0 if float(a) < 0 else float(a))
        
        # Pivot the table using the uid
        df_prep_pvt = df_prep.pivot(columns='uid')[meteo_var]
        df_prep_pvt.reset_index()
        
        # Create boolean mask - using nan values
        mask = df_prep_pvt.isnull().values
        
        # sort bolean mask - using argsort
        mask_sort = df_prep_pvt.values[np.argsort(mask, axis=0, kind='mergesort'), np.arange(mask.shape[1])]
        
        # create final dataframe
        df_prep_final = pd.DataFrame(mask_sort, columns=df_prep_pvt.columns)
        

        # Remove all rows where all values are NaN
        df_prep_final = df_prep_final.dropna(axis=0, how='all')
        
        # Ensure all values are float
        df_prep_final = df_prep_final.astype(float)
        
        # Append dataframe to list to work with it later
        _df_lst.append(df_prep_final)
    
    # Concatenate the list of dataframes to produce the final dataframe
    df_meteo = pd.concat(_df_lst, axis=0, ignore_index=False)

    # Saving file to csv
    file_name = f'{uid_field}_{uid_value}_{meteo_var}_{start_date.split("-")[0]}'\
        f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
        
    _csvFinal = Path(wrk_path / f'model_input/{file_name}.csv')
    df_meteo.to_csv(_csvFinal, index=False)
        
    log.info(f" {uid_field}_{uid_value} {meteo_var} file created")

def meteo_map(df_uid, wrk_path,
              meteo_var, uid_field,
              uid_value):
    """_summary_

    Args:
        df_uid (_type_): _description_
        wrk_path (_type_): _description_
        meteo_var (_type_): _description_
    """
    log.info(f"\n --- \n Creating {meteo_var} map file ready for SHETRAN")
    
    pivoted = df_uid.pivot(index='lat', columns='lon', values='uid')
    pivoted_order = pivoted.sort_values('lat', ascending=False)
    
    # Creating the map file
    utils.shetran_csv_file(wrk_path, f'{uid_field}_{uid_value}_{meteo_var}_map', 
                           pivoted_order, 's')
    log.info(f'{uid_field}_{uid_value} {meteo_var} map file created\n --- \n')

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

def meteo(unique_field, unique_value):
    
    # Setting the path to the work environment
    dir_abs = Path(Path().resolve() / f'{unique_field}_{unique_value}')

    # # # Make sure temp folder directory is empty
    utils.file_remove(Path(dir_abs / 'temp_data/'), 'all')    

    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')

    # Open Data Cube Product
    dc_data = config.get('dc_product', 'LC')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
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
    xmin, ymin, xmax, ymax, *df_grid_uid= mask_uid(dir_abs, crs_local, 
                                                   crs_global, unique_field,
                                                   unique_value)
    
    # Working with each meteorological dataset
    for count, value in enumerate(ERA_VAR):
        print(value)
                
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
                         df_grid_uid[0], dir_abs,
                         ND, lst_index, unique_field,
                         unique_value)
        
        # Create csv file ready for SHETRAN
        meteo_csv(dir_abs, value,
                  date_start, date_end, 
                  unique_field, unique_value)
        
        # Create map file ready for SHETRAN
        meteo_map(df_grid_uid[0], dir_abs, value,
                  unique_field, unique_value)        
    
    # Empty temp folder once more to avoid errors
    utils.file_remove(Path(Path(dir_abs / 'temp_data/')), 'all')
        
    log.info(f"\n --- \n Files for {unique_field}_{unique_value} {ERA_VAR} have been created \n"\
             f"between the periods of {date_start} and {date_end} \n --- \n ")
    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    meteo()

    
    