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

PRODUCT_NAME = 'era5_reanalysis_pet_daily'
ERA_VAR = 'potential_evapotranspiration'

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


# Set log level to info
log = setup_logging()

log.info("Starting hihydrosoils uploader")

#############################################################
# FUNCTIONS
#############################################################


def hydro_local(
    pth_file:str,
):
    # Set list to store response
    lst_cds = []
    
    # Read files locally
    p = Path(pth_file).glob('**/*')
    _files = [x for x in p if x.is_file()]

    return _files

def download_and_cog_soils(
    pth_file: str,
    s3_dst: str,
):
    soil_data = hydro_local(pth_file)
    soil_data = sorted(soil_data)
    
    # # Loop through all the soil files
    for count_cds, value_pet in enumerate(soil_data):
        file_name = Path(value_pet).name
        log.info(f"Working on {file_name}")
        
        # Put file name into a list - helps to do file naming
        lst_name = Path(value_pet).name.split("_")
        
        # create stac variables
        start_datetime = f'{lst_name[2]}-{format(lst_name[3],"02")}-{format(lst_name[4].split(".")[0],"02")}T00:00:00Z'
        end_datetime = f'{lst_name[2]}-{format(lst_name[3],"02")}-{format(lst_name[4].split(".")[0],"02")}T11:59:59Z'
        product_name = PRODUCT_NAME
        in_data = value_pet
        
        file_base = f'{s3_dst}/{lst_name[2]}/{lst_name[3]}/{PRODUCT_NAME}.{lst_name[2]}.'\
                f'{format(lst_name[3],"02")}.{format(lst_name[4].split(".")[0],"02")}'
        out_data = f"{file_base}.tif"
        out_stac = f"{file_base}.stac-item.json"

                    
        try:

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
                _, end = calendar.monthrange(int(2020), int(12))
                item = create_stac_item(
                    mem_dst,
                    id=str(odc_uuid("era5-land", "5.0", [f'{ERA_VAR}_{file_name}'])),
                    with_proj=True,
                    input_datetime=pd.to_datetime(f'{lst_name[2]}-{format(lst_name[3],"02")}-{format(lst_name[4].split(".")[0],"02")}T11:59:59Z'),
                    properties={
                        "odc:processing_datetime": datetime_to_str(datetime.now()),
                        "odc:product": product_name,
                        # "odc:region_code": lst_name[-1].split(".")[0],
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
                            target='https://doi.org/10.1038/s41597-021-01003-9',
                            title="Source file",
                            rel=pystac.RelType.DERIVED_FROM,
                            media_type="application/.nc",
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

        except Exception as e:
            message = f"Failed to handle {file_name} with error {e}"
        
#############################################################
# COMMAND LINE UTILITIES

#############################################################
@click.command("uploading HiHydrosoil to S3")
@click.option("--files_path", default="/home/geofelpave/Documents/6_Test_Area/data/test_s3_data/")
@click.option("--s3_dst", default=None,
              help="Enter the url to the S3 bucket")



def hihydrosoils(
    files_path,
    s3_dst,
):
    main_path = Path(files_path).glob('*')
    _dir = [x for x in main_path if x.is_dir()]
    lst_srt = sorted(_dir)
    for subdir in lst_srt:
        download_and_cog_soils(subdir, s3_dst)
    

# #############################################################
# # FUNCTIONS
# #############################################################

if __name__ == "__main__":
    hihydrosoils()
    
