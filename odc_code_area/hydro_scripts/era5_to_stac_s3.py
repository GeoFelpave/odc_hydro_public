# INITIAL LIBRARIES
import calendar
import json
import sys
from datetime import datetime
from typing import Optional, Tuple

import logging
from typing import Sequence
from uuid import UUID, uuid5

# libraries for working with year and months
import calendar # use to check for leap years
from calendar import monthrange

import click
import pystac
import requests
from odc.aws import s3_dump, s3_head_object
from pystac.utils import datetime_to_str
from rasterio.io import MemoryFile
from rio_cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles
from rio_stac import create_stac_item

# ERA5 LIBRARIES

# this package enables connection to cds
import climetlab as cml
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
import pprint

#############################################################
# SOURCE
"""https://github.com/digitalearthafrica/deafrica-scripts/
blob/d3c5dbbcdbe5fa2c3f7ba369096cadaa421ac9fd/deafrica/data/chirps.py#L11"""
############################################################# 

''' 
ERA 5 VARIABLES
total_precipitation - tp - era5_reanalysis_tp_daily
total_evaporation - e - era5_reanalysis_te_daily
potential_evaporation - pev - era5_reanalysis_pev_daily


'''
ERA_VAR = 'total_precipitation'
VAR_NAME = 'tp'
PRODUCT_NAME = "era5_reanalysis_tp_daily"


#############################################################
# FUNCTIONS TO BE ADDED TO A UTILS.PY FILE
############################################################# 

def setup_logging(level: int = logging.INFO) -> logging.Logger:
    """Set up a simple logger -
    https://engineeringfordatascience.com/posts/python_logging/"""
    log = logging.getLogger(__name__)
    console = logging.StreamHandler()
    log.addHandler(console)
    log.setLevel(level)
    return log

def odc_uuid(
    algorithm: str,
    algorithm_version: str,
    sources: Sequence[UUID],
    deployment_id: str = "",
    **other_tags,
) -> UUID:
    """
    Generate deterministic UUID for a derived Dataset.
    :param algorithm: Name of the algorithm
    :param algorithm_version: Version string of the algorithm
    :param sources: Sequence of input Dataset UUIDs
    :param deployment_id: Some sort of identifier for installation that performs
                          the run, for example Docker image hash, or dea module version on NCI.
    :param **other_tags: Any other identifiers necessary to uniquely identify dataset
    """
    tags = [f"{k}={str(v)}" for k, v in other_tags.items()]

    stringified_sources = (
        [str(algorithm), str(algorithm_version), str(deployment_id)]
        + sorted(tags)
        + [str(u) for u in sorted(sources)]
    )

    srcs_hashes = "\n".join(s.lower() for s in stringified_sources)
    return uuid5(UUID("6f34c6f4-13d6-43c0-8e4e-42b6c13203af"), srcs_hashes)

def file_remove(file_check_path):
    dir_abs = Path().resolve()
        
    file_name = file_check_path.stem
    file_path = file_check_path.parent
    lst_extensions = ['cpg', 'dbf', 'prj', 'shp', 'shx', 'tif']

    for extension in lst_extensions:
        file_other_name= f'{file_name}.{extension}'
        file_to_check = Path(dir_abs / Path(file_path / file_other_name))
        
        if file_to_check.is_file():
            file_to_check.unlink()
            print(f'{file_other_name} DELETED!')

# Set log level to info
log = setup_logging()

log.info("Starting ERA5 downloader")

#############################################################
# FUNCTIONS
#############################################################

def check_values(
    year: str, month: str, day: Optional[str]
) -> Tuple[str, str, Optional[str]]:
    # Assert that the year should be 4 characters long
    assert len(year) == 4, "Year should be 4 characters long"

    # Ensure the month is zero padded   
    if len(month) == 1:
        month = f"0{month}"

    # As is the day, if it is provided
    if day is not None and len(day) == 1:
        day = f"0{day}"

    return year, month, day

def check_for_url_existence(href):
    response = requests.head(href)
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError:
        log.error(f"{href} returned {response.status_code}")
        return False
    return True


def cds_response(variable_name, year_data, month_data, day_data, lst_horas):
    ## Getting the data from cds
    print(variable_name, year_data, month_data, day_data, lst_horas)
    source_prep = cml.load_source(
        "cds",
        "reanalysis-era5-land",
        variable=variable_name,
        product_type="reanalysis",
        area=[90, -180, -90, 180], #[latmax, longmin, latmin, lonmax] - this is for the World
        # area=[13, -79, -4.5, -66], #[latmax, longmin, latmin, lonmax] - this is for Colombia
        # area=[61, -10.80, 49.80, 2], #[latmax, longmin, latmin, lonmax] - this is for UK
        year = year_data,
        month = month_data,
        day = day_data,
        time = lst_horas,
        format="netcdf",
    )
    return source_prep

def cds_netcdf(
    years: list,
    months: list,
    hours: list,
    era_var: str,
):
    # Set list to store response
    lst_cds = []
    
    # Calling the data from the START_DATE to END_DATE7
    # Each individual response is stored in a list to be called later
    for year in years:
        for month in months:
            # Create list of days based on the year and month
            # this will deal with leap yeards too
            days_of_month = monthrange(year, month)
            
            '''the list - range need to start at one, because
            calendar.monthrange(year, month) returns weekday of 
            first day of the month and number of days in month, 
            for the specified year and month. '''
            lst_days = list(range(1, days_of_month[-1] + 1, 1))
            
            # Get data from CDS and store to be used later on
            # this will work for temperature if needed
            lst_cds.append(cds_response(era_var, year, month, lst_days, hours))

            # for dia in lst_days:
            #     print(era_var, year, month, dia, hours)
                # lst_cds.append(cds_response(era_var, year, month, dia, hours))
    return lst_cds

def download_and_cog_era5(
    years: list,
    months: list,
    hours: list,
    era_var: str,
    s3_dst: str,
):
    era_data = cds_netcdf(years, months, hours, era_var)
        
    # Cleaning and sanity checks
    s3_dst = s3_dst.rstrip("/")
    
    # Loop through all the CDS API responses 
    for count_cds, value_cds in enumerate(era_data):
        
        # Open nc file using 'with' for statbility
        # resample from hourly to daily using 
        with xr.open_dataset(value_cds) as ds:
            array_data = ds
                
        # Create list of timestamps
        ## lst_time = xr.open_dataset(daily_precip).coords['time'].values
        lst_time = array_data.coords['time'].values
                    
        # Loop inside each nc file using the timestamp
        for count, value in enumerate(lst_time):

            # for precipitation only
            date_value = value - np.timedelta64(1, 'D') # this is needed only for var with accumulation values - see ERA5 Land Accumulation
            
            # Read nc file as xarray and select one timestamp at the time
            array_tif = array_data.isel(time=slice(count, count + 1))
            
            # Make sure to change the time in xarray so it is reflected in .tif
            array_tif.coords['time'] = array_tif.coords['time'] - np.timedelta64(1, 'D')
            
            # Set CRS
            array_tif = array_tif.rio.write_crs("epsg:4326", inplace=True)        

            # Change values from meters to mm
            array_tif = array_tif * 1000
            
            # Create tif file
            # Create file name using the date
            file_name = np.datetime_as_string(date_value, unit='h')
            file_name = file_name.replace('-', "")
            array_tif[VAR_NAME].rio.to_raster(f'temp/{ERA_VAR}_{file_name}.tif')
            
            # Set variables for file names
            y = pd.to_datetime(date_value).year
            m = pd.to_datetime(date_value).month
            d = pd.to_datetime(date_value).day
            h = pd.to_datetime(date_value).hour
            
            file_base = f"{s3_dst}/{y}/{m}/{PRODUCT_NAME}.{y}."\
                f"{format(m,'02')}.{format(d,'02')}"
            out_data = f"{file_base}.tif"
            out_stac = f"{file_base}.stac-item.json"

            start_datetime = f"{y}-{format(m,'02')}-{format(d,'02')}T{format(h,'02')}:00:00Z"
            end_datetime = f"{y}-{format(m,'02')}-{format(d,'02')}T{format(h,'02')}:59:59Z"
            product_name = PRODUCT_NAME
            in_data = f'temp/{ERA_VAR}_{file_name}.tif'
            
            try:
                # Check if file already exists
                log.info(f"Working on {value}")
                # if not overwrite and s3_head_object(out_stac) is not None:
                #     log.warning(f"File {out_stac} already exists. Skipping.")
                #     return

                # COG and STAC
                with MemoryFile() as mem_dst:
                    
                    # Creating the COG, with a memory cache and no download. Shiny.
                    cog_translate(
                        in_data,
                        mem_dst.name,
                        cog_profiles.get("deflate"),
                        in_memory=True, 
                        nodata=-9999,
                    )
                    
                    # Creating the STAC document with appropriate date range
                    _, end = calendar.monthrange(int(y), int(m))
                    item = create_stac_item(
                        mem_dst,
                        id=str(odc_uuid("era5-land", "5.0", [f'{ERA_VAR}_{file_name}'])),
                        with_proj=True,
                        input_datetime=pd.to_datetime(value - np.timedelta64(1, 'D')),
                        properties={
                            "odc:processing_datetime": datetime_to_str(datetime.now()),
                            "odc:product": product_name,
                            "start_datetime": start_datetime,
                            "end_datetime": end_datetime,
                        },
                    )
                    
                    item.set_self_href(out_stac)
                    # Manually redo the asset
                    del item.assets["asset"]
                    item.assets[ERA_VAR] = pystac.Asset(
                        href=out_data,
                        title="ERA5-reanalysis-daily",
                        media_type=pystac.MediaType.COG,
                        roles=["data"],
                    )
                    # Let's add a link to the source
                    item.add_links(
                        [
                            pystac.Link(
                                target='https://cds.climate.copernicus.eu/cdsapp#!/home',
                                title="Source file",
                                rel=pystac.RelType.DERIVED_FROM,
                                media_type="application/netcdf",
                            )
                        ]
                    )
                    # Dump the data to S3
                    mem_dst.seek(0)
                    log.info(f"Writing DATA to: {out_data}")
                    s3_dump(mem_dst, out_data, ACL="bucket-owner-full-control")
                    # Write STAC to S3
                    log.info(f"Writing STAC to: {out_stac}")
                    s3_dump(
                        json.dumps(item.to_dict(), indent=2),
                        out_stac,
                        ContentType="application/json",
                        ACL="bucket-owner-full-control",
                    )
                    # All done!
                    # Delete tif file to avoid overloading
                    file_remove(Path(f'temp/{ERA_VAR}_{file_name}.tif'))
                    log.info(f"Completed work on {value} \n ####################################")

            except Exception as e:
                message = f"Failed to handle {value} with error {e}"



#############################################################
# COMMAND LINE UTILITIES

#############################################################
@click.command("download-era5")
@click.option("--start_date", default="1980-01-01")
@click.option("--end_date", default="1980-01-01",
              help="the end date needs to be n+1 due to way ERA5 Land data is retrieved")
@click.option("--start_hour", default="00:00:00")
@click.option("--end_hour", default="00:00:00")
@click.option("--s3_dst", default=None, 
              help="Enter the url to the S3 bucket")
# @click.option("--overwrite", is_flag=True, default=False)


def reanalysis_era5_land(
    start_date, 
    end_date,
    start_hour, 
    end_hour,
    s3_dst
):
    start_date = datetime.strptime(start_date, '%Y-%m-%d')
    end_date = datetime.strptime(end_date, '%Y-%m-%d')
    
    # set months in the right order
    mth_strt = start_date.month
    mth_end = end_date.month
    if mth_strt > mth_end:
        mth_strt = end_date.month
        mth_end = start_date.month
    
    # Make sure the months are done correctly
    # This will do a whole year for different periods in different years
    # For periods with the same year, it will do the months asked for
    if (start_date.year != end_date.year):
        lst_months = list(range(mth_strt, 13, 1))
        # Create list of years, month, day and hours
        if mth_end != 12:
            lst_years = list(range(start_date.year, end_date.year + 1, 1))
        else:
            lst_years = list(range(start_date.year, end_date.year, 1))
    else:
        lst_months = list(range(mth_strt, mth_end + 1, 1))
        # Create list of years, month, day and hours
        lst_years = list(range(start_date.year, end_date.year + 1, 1))
        
    lst_hours = pd.date_range(start_hour,end_hour,freq="60min").time
    lst_hours = [t.strftime('%H:%M:%S') for t in lst_hours]

    download_and_cog_era5(lst_years, lst_months, lst_hours, ERA_VAR, s3_dst)
    

# #############################################################
# # FUNCTIONS
# #############################################################

if __name__ == "__main__":
    reanalysis_era5_land()
    
    # # Set the start and end datas
    # START_DATE = datetime.strptime('2020-01-01', '%Y-%m-%d')
    # END_DATE = datetime.strptime('2020-12-31', '%Y-%m-%d')
    # START_HOUR = '00:00:00'
    # END_HOUR = '23:00:00'
    
    # reanalysis_era5_land(START_DATE, END_DATE, START_HOUR, END_HOUR)
    
