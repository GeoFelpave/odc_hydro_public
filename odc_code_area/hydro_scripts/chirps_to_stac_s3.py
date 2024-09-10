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

#############################################################
# SOURCE
"""https://github.com/digitalearthafrica/deafrica-scripts/
blob/d3c5dbbcdbe5fa2c3f7ba369096cadaa421ac9fd/deafrica/data/chirps.py#L11"""
############################################################# 



DAILY_URL_TEMPLATE = "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/tifs/p05//{year}/{in_file}"


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

log.info("Starting CHIRPS downloader")

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

def download_and_cog_chirps(
    year: str,
    month: str,
    s3_dst: str,
    day: str = None,
):
    # Cleaning and sanity checks
    s3_dst = s3_dst.rstrip("/")
    # Set up file strings
    if day is not None:
        # Set up a daily process
        in_file = f"chirps-v2.0.{year}.{month}.{day}.tif.gz"
        in_href = DAILY_URL_TEMPLATE.format(year=year, in_file=in_file)
        
        """/vsigzip/ is a file handler that allows on-the-fly reading of GZip
        (.gz) files without decompressing them in advance -GDAL"""
        
        in_data = f"/vsigzip//vsicurl/{in_href}"

        if not check_for_url_existence(in_href):
            log.warning("Couldn't find the gzipped file, trying the .tif")
            in_file = f"chirps-v2.0.{year}.{month}.{day}.tif"
            in_href = DAILY_URL_TEMPLATE.format(year=year, in_file=in_file)
            in_data = f"/vsicurl/{in_href}"

            if not check_for_url_existence(in_href):
                log.error("Couldn't find the .tif file either, aborting")
                sys.exit(1)

        file_base = f"{s3_dst}/{year}/{month}/chirps-v2.0_{year}.{month}.{day}"
        out_data = f"{file_base}.tif"
        out_stac = f"{file_base}.stac-item.json"

        start_datetime = f"{year}-{month}-{day}T00:00:00Z"
        end_datetime = f"{year}-{month}-{day}T23:59:59Z"
        product_name = "rainfall_chirps_daily"
    
    try:
        # Check if file already exists
        log.info(f"Working on {in_file}")
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
            _, end = calendar.monthrange(int(year), int(month))
            item = create_stac_item(
                mem_dst,
                id=str(odc_uuid("chirps", "2.0", [in_file])),
                with_proj=True,
                input_datetime=datetime(int(year), int(month), int(day)),
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
            item.assets["rainfall"] = pystac.Asset(
                href=out_data,
                title="CHIRPS-v2.0",
                media_type=pystac.MediaType.COG,
                roles=["data"],
            )
            # Let's add a link to the source
            item.add_links(
                [
                    pystac.Link(
                        target=in_href,
                        title="Source file",
                        rel=pystac.RelType.DERIVED_FROM,
                        media_type="application/gzip",
                    )
                ]
            )
            # Use statment below to check what the STAC json looks like
            # print(json.dumps(item.to_dict(), indent=4))

            # Dump the data to S3
            # 
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
            log.info(f"Completed work on {in_file}")

    except Exception as e:
        message = f"Failed to handle {in_file} with error {e}"


#############################################################
# COMMAN LINE UTILITIES

#############################################################
# @click.command("download-chirps")
# @click.option("--year", default="2020")
# @click.option("--month", default="01")
# @click.option("--day", default="01")
# @click.option("--s3_dst", default="s3://hydro-odc/chirps/")
# @click.option("--overwrite", is_flag=True, default=False)

def cli_daily(year, month, day, s3_dst):
    """
    Download CHIRPS Africa daily tifs, COG, copy to
    S3 bucket.
    GeoTIFFs are copied from here:
        https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/tifs/p05/
    Example:
    download-chirps-daily
        --s3_dst s3://deafrica-data-dev-af/rainfall_chirps_monthy/
        --year 1983
        --month 01
        --day 01
    Available years are 1981-2021.
    """

    yr, mt, dt = check_values(year, month, day)

    download_and_cog_chirps(
        year=yr,
        month=mt,
        day=dt,
        s3_dst=s3_dst,
    )

#############################################################
# FUNCTIONS
#############################################################

if __name__ == "__main__":
    START_DATE = datetime.strptime('2020-01', '%Y-%m')
    END_DATE = datetime.strptime('2020-01', '%Y-%m')
    
    # Create list of years, month, day and hours
    lst_years = list(range(START_DATE.year, END_DATE.year + 1, 1))
    lst_months = list(range(START_DATE.month, END_DATE.month + 1, 1))
    
    for year in lst_years:
        for month in lst_months:
            # Create list of days based on the year and month
            # this will deal with leap yeards too
            days_of_month = monthrange(year, month)
            
            '''the list - range need to start at one, because
            calendar.monthrange(year, month) returns weekday of 
            first day of the month and number of days in month, 
            for the specified year and month. '''
            lst_days = list(range(1, days_of_month[-1] + 1, 1))
            
            for dia in lst_days:
                y, m, d = (str(year), str(month), str(dia))
                cli_daily(y,m,d,"s3://odc-hydro-data/chirps-global/")

    """Use the lines below to test
    - cli_daily()
    - setup_logging()
    - odc_uuid("chirps", "2.0", ["chirps-v2.0.2020.01.01.tif.gz"])
    
    The code can be run in terminal using this line to pass parameters --year or --month
    when using the command line utilities code
    
    python chirps_to_stac.py --day "01" - """