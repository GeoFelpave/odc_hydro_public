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
                     f'active_data/{uid_field}_{uid_value}_grid_WGS84.shp')
    
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

    # # Replacing 0 with -9999
    # grid.loc[grid[variable_name] == 0, variable_name] = -9999
    
    # delte temp .tif file
    utils.file_remove(Path(wrk_path / 'temp_data/'), 'tif')

    return df_mask

def soil_catchm_dataframe(soil_grid_df, crs_lcl, soil_var, soil_depth):
    # Get the centroid for each grid
    temp_df = pd.DataFrame()
    temp_df = soil_grid_df.to_crs(crs_lcl)
    temp_df['centroid'] = temp_df['geometry'].centroid

    # Create lat (Y) and lon (X) columns
    temp_df['lat'] = temp_df['centroid'].y.astype(int)
    temp_df['lon'] = temp_df['centroid'].x.astype(int)

    # Change geopandas to pandas ready to create csv file for the mean and min elevation
    # df_grid_lc = pd.DataFrame(temp_df[['lat', 'lon', 'property_value']].copy())
    df_grid_lc = pd.DataFrame(temp_df[['lat', 'lon', 'SHETRAN_ID', 'property_value']].copy())
    
    # add the soil depth for reference
    df_grid_lc['soil_depth'] = soil_depth
    
    # add the soil property name for reference
    df_grid_lc['property_name'] = soil_var
    
    return df_grid_lc

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
    
    # Prepare soil property dataframe
    # set var name if doing geology
    if product_name == 'ghlymps_logKferr_200_2000cm':
        # Set variables
        var_name = 'Ksat'
    
        
    df_soil = soil_catchm_dataframe(df_properties, crs_lcl,
                                    var_name, depth)
    
    return df_soil
    
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

def soil_column(lst_dataframes):
    log.info("\n --- \n Creating data for the soil column\n --- \n")
    
    # combine all dataframe
    df_soil_all = pd.concat(lst_dataframes)
    df_soil_all = df_soil_all.reset_index()
    
    # Round values before contiuing with process
    df_soil_all['property_value'] = df_soil_all['property_value'].round(2)
    
    # Pivot table - soil properties as columns
    pivoted = df_soil_all.pivot(index=['lat','lon', 'soil_depth', 'SHETRAN_ID'], 
                                columns='property_name', 
                                values='property_value').reset_index()
    
    # Order pivot table
    pivoted_order = pivoted.copy()
    pivoted_order = pivoted_order.sort_values(by=['lat', 'lon', 'soil_depth'])


    ## Create the depth column
    layer_order = [1, 2, 3, 4, 5, 6, 7]
    depth_conditions = [
        pivoted_order['soil_depth'] == '0_5cm', 
        pivoted_order['soil_depth'] == '5_15cm',
        pivoted_order['soil_depth'] == '15_30cm',
        pivoted_order['soil_depth'] == '30_60cm',
        pivoted_order['soil_depth'] == '60_100cm',
        pivoted_order['soil_depth'] == '100_200cm',
        pivoted_order['soil_depth'] == '200_2000cm'
        ]

    # Adding the soil depth column
    pivoted_order['layer_order'] = np.select(depth_conditions, layer_order)

    # Sort by category and soil layer
    pivoted_order = pivoted_order.sort_values(by=['lat', 'lon', 'layer_order'])
    pivoted_order.reset_index(drop=True, inplace=True)

    # Drop unnecessary columns
    pivoted_order = pivoted_order.drop('layer_order', axis=1)
    
    # Set of columns needed
    lst_cols =['lat', 'lon', 'ALFA', 'Ksat', 'N', 'WCres', 'WCsat', 'soil_depth', 'SHETRAN_ID'] 

    # Copy dataframe to avoid any error or crashes
    df_soils_types = pivoted_order[lst_cols].copy().reset_index(drop=True)

    # Create list of soil types ID
    # This create the list without comparing. 
    # The comparison between soil properties will happen at latter stage
    lst_soil_type = []
    for i in range(0, len(df_soils_types)):
        # lst_soil_type.append(f'st_{str(i).zfill(5)}')
        lst_soil_type.append(i+1)
        
    # Sort by category and soil layer
    # df_soils_types = df_soils_types.sort_values(by=['soil_depth'])

    # Add list to dataframe
    df_soils_types['soil_type_ID'] = lst_soil_type

    # Add soil category column
    df_soils_types['soil_category_ID'] = np.nan

    # df_soils_types.sort_index()

    df_soils_types

    return df_soils_types

def soil_type_categ(df_soilColumn):
    log.info("\n --- \n Defining soil type and category using the soil column information\n --- \n")
    
    df_soils_types_final = df_soilColumn.copy()

    # List of soild characteristics names
    lst_cols_charact =['ALFA', 'Ksat', 'N', 'WCres', 'WCsat']

    # Create dataframe with soil characteristics only
    df_check = df_soilColumn[lst_cols_charact].copy().reset_index(drop=True)
    df_check
    
    # List of soild type ID and category ID
    # This will be used to avoid repetition in a further loop
    lst_soil_type_ID = []
    # lst_soil_ctgr_ID = []
    cnt_cat = 1

    # Creating the process of comparison through the whole dataframe
    # for i in range(20):
    for i in range(len(df_soilColumn)):
        # empty variables
        tmp_lst = []
        tmo_df = []
        
        # Populate temp dataframe and remove row for dataframe
        tmp_df = pd.DataFrame()
        tmp_df = df_check.copy()
        
        #######################################################
        ### WORKING WITH THE SOIL TYPE ###
        #######################################################
        
        # Create soil type ID variable
        soil_type_ID = df_soils_types_final.loc[i, ['soil_type_ID']].values[0]
        
        # Get row as a list and remove from dataframe
        tmp_lst = tmp_df.values.tolist()[i]
        
        # Remove the row (tmp_lst) from dataframe before comparison
        tmp_df = tmp_df.drop(index=i)
        
        # Check if the set of soil characteristics (tmp_lst) appears again in the dataframe
        '''This takes the list of soil characteristics and compares it to each row of the dataframe individually. It returns
        a dataframe with the rows that are equal to the list'''
        duplicate_df = df_check.loc[df_check.isin(tmp_lst).astype(int).sum(axis=1) == len(tmp_lst), :]
        # duplicate_df = tmp_df.loc[tmp_df.isin(tmp_lst).astype(int).sum(axis=1) == len(tmp_lst), :]
        
        # check if duplicate df is not empty - returns true of empty
        flag = duplicate_df.empty
        
        # if the duplicate_df is not empty then make sure the rows have the same soil type
        if flag is False:
            
            # Check the soil type with the lst of soil types already checked
            if (soil_type_ID not in lst_soil_type_ID):
                '''For example if we are in row n+1 and the soild ID has already been changed for the one in 
                row n, then we don't want to redo that row again'''
                
                # Store soil type ID in list to avoid repetition in future loops
                lst_soil_type_ID.append(soil_type_ID)
                
                # Add the soil type ID to the rows with equal soil characteristics
                '''Here the code uses de index inside the soil duplicate dataframe to change the soil_type_ID
                in the final dataframe using the index and the soil_type_ID of the value that is currently being read'''
                for value in list(duplicate_df.index.values):
                    df_soils_types_final.loc[value, 'soil_type_ID'] = soil_type_ID
        
        #######################################################
        ### WORKING WITH THE SOIL CATEGORY###
        ####################################################### 
        
        # Create soil type ID variable
        soil_category_ID = df_soils_types_final.loc[i, ['soil_category_ID']].values[0]
        
        # Create coordinates (soild category) as a list - ready for comparing
        lst_soil_categor = list(df_soils_types_final.loc[i, ['lat', 'lon']])
        
        # Find the grids where the set of coordinates are duplicate 
        '''This takes the list of soil characteristics and compares it to each row of the dataframe individually. It returns
        a dataframe with the rows that are equal to the list based on the lat and lon i.e the different soil depths'''
        
        grid_depths = df_soils_types_final.loc[df_soils_types_final.isin(lst_soil_categor).sum(axis=1) == len(lst_soil_categor), :]
        
        # check if duplicate df is not empty - returns true of empty
        flag_cat = grid_depths.empty
        
        # if the duplicate_df is not empty then make sure the rows have the same soil type
        
        if flag_cat is False:
            # Check the category is not nan
            if (pd.isna(soil_category_ID) is True):
                # print(grid_depths)
                            
                # Add the soil category ID to the rows with equal coordinates / grid
                '''Here the code uses de index inside the soil duplicate dataframe to change the soil_category_ID
                in the final dataframe using the index and the soil_category_ID of the value that is currently being read'''
                
                for value in list(grid_depths.index.values):
                    # print(value)
                    df_soils_types_final.loc[value, 'soil_category_ID'] = cnt_cat
                    
                # Increase the soil category ID
                cnt_cat += 1
    
    # Make sure category ID is Integer
    df_soils_types_final['soil_category_ID'] = df_soils_types_final['soil_category_ID'].astype(int)

    log.info("\n --- \n Soil type and category attributes created\n --- \n")
    return df_soils_types_final

def project_file(df_soil_types_category, wrk_path,
                 uid_field, uid_value,
                 start_date, end_date):
    
    log.info("\n --- \n Creating final project file\n")
    
    df_xml_v1 = df_soil_types_category.copy()

    ## 1. Create the depth column
    depth_values = [0.05, 0.15, 0.30, 0.60, 1.00, 2.00, 20.00]
    layer_values = [1, 2, 3, 4, 5, 6, 7]
    depth_conditions = [
        df_xml_v1['soil_depth'] == '0_5cm', 
        df_xml_v1['soil_depth'] == '5_15cm',
        df_xml_v1['soil_depth'] == '15_30cm',
        df_xml_v1['soil_depth'] == '30_60cm',
        df_xml_v1['soil_depth'] == '60_100cm',
        df_xml_v1['soil_depth'] == '100_200cm',
        df_xml_v1['soil_depth'] == '200_2000cm'
        ]

    # 2. Adding the soil depth column
    df_xml_v1['soil_depth_meters'] = np.select(depth_conditions, depth_values)

    # 3. Adding the layer value column
    df_xml_v1['soil_layer'] = np.select(depth_conditions, layer_values)

    # 4. Add soil number column based on the index
    df_xml_v1['soil_number'] = df_xml_v1.index.values + 1


    # 5. Create soil properties column ready for xml
    lst_soil_properties = ['soil_number', 'soil_type_ID',
                           'WCsat', 'WCres', 'Ksat',
                           'ALFA','N']

    # 6. Create soil properties column ready for xml
    lst_soil_details = ['soil_category_ID', 'soil_layer',
                        'soil_type_ID', 'soil_depth_meters']

    # 7. Create Soil Property and Detail columns
    df_xml_v1['SoilProperty'] = df_xml_v1[lst_soil_properties].values.tolist()
    df_xml_v1['SoilDetail'] = df_xml_v1[lst_soil_details].values.tolist()

    # 8. Sort by category and soil layer
    df_xml_v1 = df_xml_v1.sort_values(by=['soil_category_ID', 'soil_layer'])
    
    # 9. Change soil property and soil detail to list
    '''Changing the column to list generates a list of lists. Before converting to XML
    it is necessary to change each list to a single string element separated by a comma.
    The final output should be a list where each element is the list of soil properties
    as a string. It should go from a list of lists to a simple list of elements'''

    lst_soilProp_header = config.get('xml', 'soilProp_header').split(',')
    lst_soilDet_header = config.get('xml', 'soilDet_header').split(',')

    # 10. Soil Property to list of list
    lst_soilPrprty_xml = list(df_xml_v1.SoilProperty)

    # 11. Changing to list of list ready for dicitonary
    for i, lst_item in enumerate(lst_soilPrprty_xml):
        lst_soilPrprty_xml[i] = ", ".join([str(item) for item in lst_item])

    # 12. Soil Detail to list of list
    lst_soilDetail_xml = list(df_xml_v1.SoilDetail)

    # 13. Changing to list of list ready for dicitonary
    for i, lst_item in enumerate(lst_soilDetail_xml):
        lst_soilDetail_xml[i] = ", ".join([str(item) for item in lst_item])

    # 14. Add header to the list of soil properties
    lst_soilPrprty_xml = [", ".join([str(item) for item in lst_soilProp_header])] + lst_soilPrprty_xml

    # 15. Add header to the list of soil details
    lst_soilDetail_xml = [", ".join([str(item) for item in lst_soilDet_header])] + lst_soilDetail_xml
    
    # 16. Create dicts ready for XML
    # Dict flie tags
    dict_tags = {
        'ProjectFile':f'{uid_field}_{uid_value}_LibraryFile',
        'CatchmentName':f'{uid_field}_{uid_value}',
        'DEMMeanFileName':f'{uid_field}_{uid_value}_dem_mean.txt',
        'DEMminFileName':f'{uid_field}_{uid_value}_dem_min.txt',
        'MaskFileName':f'{uid_field}_{uid_value}_mask.txt',
        'VegMap':f'{uid_field}_{uid_value}_lc.txt',
        'SoilMap':f'{uid_field}_{uid_value}_subsurface.txt',
        'LakeMap':f'{uid_field}_{uid_value}_lakes.txt',
        'PrecipMap':f'{uid_field}_{uid_value}_total_precipitation_map.txt',
        'PeMap':f'{uid_field}_{uid_value}_potential_evapotranspiration_map.txt'
    }

    # Vegetation details
    '''This will need to be reviewed and coded rather than having it hard coded into the script'''
    # lst_vegdetails = [
    #     'VegetationTypeNumber, VegetationType,'\
    #     'CanopyStorageCapacity, LeafAreaIndex,'\
    #     'MaximumRootingDepth, AE/PEatFieldCapacity,'
    #     'StricklerCoefficient',
    #     '1, Arable, 1.5, 1, 0.8, 0.6, 2.5',
    #     '2, BareGround, 0, 0, 0.1, 0.4, 3.0',
    #     '3, Grass, 1.5, 1, 1.0, 0.6, 0.5',
    #     '4, Forest, 5, 1, 1.8, 1.0, 1.0',
    #     '5, FloodedVegetation, 0, 0, 0.1, 0.4, 3.0',
    #     '6, Urban, 0.3, 0.3, 0.5, 0.4, 5.0'
    # ]
    lst_vegdetails = [
        'VegetationTypeNumber, VegetationType,'\
        'CanopyStorageCapacity, LeafAreaIndex,'\
        'MaximumRootingDepth, AE/PEatFieldCapacity,'
        'StricklerCoefficient',
        '1, Arable, 1.0, 0.8, 0.8, 0.6, 0.6',
        '2, BareGround, 0, 0, 0.1, 0.4, 3.0',
        '3, Grass, 1.5, 1.0, 1.0, 0.6, 0.5',
        '4, Forest, 5.0, 1.0, 2.0, 1.0, 0.25',
        '5, FloodedVegetation, 0, 0, 0.1, 0.4, 3.0',
        '6, Urban, 0.3, 0.3, 0.5, 0.4, 5.0'
    ]

    tp_file_name = f'{uid_field}_{uid_value}_total_precipitation_{start_date.split("-")[0]}'\
        f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
        
    pe_file_name = f'{uid_field}_{uid_value}_potential_evapotranspiration_{start_date.split("-")[0]}'\
    f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
    
    ## Settign name for temp files in XML 
    # tmax_file_name = f'{uid_field}_{uid_value}_2m_temperature_max_{start_date.split("-")[0]}'\
    # f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
    
    # tmin_file_name = f'{uid_field}_{uid_value}_2m_temperature_min_{start_date.split("-")[0]}'\
    # f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
    
    tmean_file_name = f'{uid_field}_{uid_value}_2m_temperature_mean_{start_date.split("-")[0]}'\
    f'{start_date.split("-")[1]}_{end_date.split("-")[0]}{end_date.split("-")[1]}'
    
    # ensuring the date is correct in the xml file
    """Shetran needs the n day + 1 for the end of the simulation"""
    
    
    #################################
    # This step below needs reviewing
    #################################
    end_d = end_date.split("-")[2]
    end_m = end_date.split("-")[1]
    end_y = end_date.split("-")[0]
    
    # if end_d == '31' and end_m == '12':
    #     end_d = 1
    #     end_m = 1
    #     end_y = int(end_y) + 1
    # else:
    #     end_d = int(end_d) + 1
    #     end_m = int(end_m) + 1
    #################################
    ################################# 
     
    end_dict_tags = {
        'InitialConditions':'0',
        'PrecipitationTimeSeriesData':f'{tp_file_name}.csv',
        'PrecipitationTimeStep':24,
        'EvaporationTimeSeriesData':f'{pe_file_name}.csv',
        'EvaporationTimeStep':24,
        # 'MaxTempTimeSeriesData':f'{tmax_file_name}.csv',
        # 'MinTempTimeSeriesData':f'{tmin_file_name}.csv',
        'MaxTempTimeSeriesData':f'{tmean_file_name}.csv',
        'MinTempTimeSeriesData':f'{tmean_file_name}.csv',
        'StartDay':f'{start_date.split("-")[2]}',
        'StartMonth':f'{start_date.split("-")[1]}',
        'StartYear':f'{start_date.split("-")[0]}',
        'EndDay':f'{str(end_d).zfill(2)}',
        'EndMonth':f'{str(end_m).zfill(2)}',
        'EndYear':f'{str(end_y).zfill(2)}',
        # 'EndDay':f'{end_date.split("-")[2]}',
        # 'EndMonth':f'{end_date.split("-")[1]}',
        # 'EndYear':f'{end_date.split("-")[0]}',
        'RiverGridSquaresAccumulated':'2',
        'DropFromGridToChannelDepth':'2',
        'MinimumDropBetweenChannels':'0.5',
        'RegularTimestep':'1.0',
        'IncreasingTimestep':'0.05',
        'SimulatedDischargeTimestep':'24.0',
        'SnowmeltDegreeDayFactor':'0.0002'    
    }

    # 17. Add vegetation detail tag to dict
    dict_tags.update(VegetationDetails=lst_vegdetails)

    # 18. Add soilt property tag to dict
    dict_tags.update(SoilProperties=lst_soilPrprty_xml)

    # 19. Add soilt detail tag to dict
    dict_tags.update(SoilDetails=lst_soilDetail_xml)

    # 20. Add end tag to dict
    dict_tags.update(end_dict_tags)
    
    # 21. Form XML file
    # Conver dict to XML
    my_item_func = lambda x: 'SoilProperty' if x == 'SoilProperties' else ('SoilDetail' if x == 'SoilDetails' else 'VegetationDetail')
    xml = dicttoxml(dict_tags, custom_root='ShetranInput', item_func = my_item_func, attr_type = False)

    # 22. Saving to a file
    # make sure to remove any indent - needed for the executable
    root = ET.fromstring(xml)
    tree = ET.ElementTree(root)
    ET.indent(tree, space="")

    treestring = "<?xml version='1.0'?>" + ET.tostring(root, encoding='unicode', method='xml')
    with open(Path(wrk_path / f'model_input/{uid_field}_{uid_value}_LibraryFile.xml'), "w") as f:
        f.write(treestring)
    
    log.info("Project file created\n --- \n")

def subsurface_map(df_soil_types_category, wrk_path,
                   nd_value,uid_field,uid_value):
    
    log.info("\n --- \n Creating subsurface map\n")
    
    # 1. Copy df soil types to avoid any issues
    df_subsurface = df_soil_types_category.copy()

    # 2. Change the soil category for SHETRAN_ID equal to -9999
    '''This is needed as the subsurface should be the same as the mask'''
    df_subsurface.loc[df_subsurface.SHETRAN_ID == nd_value, 'soil_category_ID'] = nd_value

    # 3. Drop any duplicates
    '''This is needed as each grid has 7 rows'''
    df_subsurface = df_subsurface.drop_duplicates(subset=['lat', 'lon', 'soil_category_ID'])

    # 4. Only work with the columns needed
    df_subsurface = pd.DataFrame(df_subsurface[['soil_category_ID', 'lat', 'lon']].copy())


    # Pivoting dataframe to replicate SHETRAN format
    # Pivoting dataframe using lon as column and lat as row
    df_subsurface_pivot = df_subsurface.pivot(index='lat', columns='lon', values='soil_category_ID')
    df_subsurface_pivot = df_subsurface_pivot.sort_index(ascending=False).round(0).astype(int)

    df_subsurface_pivot
    
    # Creating the subsurface map file
    utils.shetran_csv_file(wrk_path, f'{uid_field}_{uid_value}_subsurface', df_subsurface_pivot, 'd')
    log.info("SHETRAN Subsurface file created\n --- \n")

#############################################################
# Logging to console
#############################################################
# Set log level to info
log = setup_logging()

log.info("\n --- \n Creating SHETRAN subsurface files\n --- \n")


#############################################################
# COMMAND LINE UTILITIES
#############################################################
@click.command("creating soil data map and config file")
@click.option("--unique_field", 
              default='HYBAS_ID', 
              help="Enter the field to be used as unique ids - if using "
                  "HydroSHEDS then HYBAS_ID will be used as default")
@click.option("--unique_value", 
              default=None, 
              help="Enter the unique identifier for the catchment")

def surface(unique_field, unique_value):
    
    # Setting the path to the work environment
    dir_abs = Path(Path().resolve() / f'{unique_field}_{unique_value}')

    # Make sure temp folder directory is empty
    utils.file_remove(Path(dir_abs / 'temp_data/'), 'all')

    # Read config file values
    config.read('config.ini')

    # Setting CRS
    crs_global = config.getint('crs_setting', 'GLB')
    crs_local = config.getint('crs_setting', 'LCL')

    # No data value
    ND = config.getint('res_setting', 'NO_DATA')
    
    # Get the dates for the file name
    date_start = config.get('time_period', 'start_date')
    date_end = config.get('time_period', 'end_date')

    # Open Data Cube Product
    lst_dc_data = config.get('dc_product', 'lst_soil_products').split(',')
    
    # Read shp file to geopandas dataframe
    grid_path = Path(dir_abs / 
                     f'active_data/{unique_field}_{unique_value}_grid_WGS84.shp')
    grid_soil = gpd.read_file(grid_path)

    xmin, ymin, xmax, ymax = grid_soil.total_bounds
    
    # Getting data from ODC and process the soil data from Hihydrosoil
    # Create column of values using hihydrosil data
    
    # list to stor dataframe
    lst_df = []
    
    for count, value in enumerate(lst_dc_data):
        
        df_soil = soil_properties(value, xmin, 
                                  ymin, xmax, ymax,
                                  grid_soil, dir_abs,
                                  ND, crs_local,
                                  unique_field, unique_value)
        
        # Append dataframe to list
        lst_df.append(df_soil)
            
    # Set geology
    # Adding geology dataframe to the list of dataframes
    geo_product = 'ghlymps_logKferr_200_2000cm'
    lst_soil_column = geology(geo_product, xmin, 
                           ymin, xmax, ymax,
                           grid_soil, dir_abs,
                           ND, crs_local, lst_df,
                           unique_field, unique_value)
    
    # Prepare dataframe for soil type and categories
    df_soil_column = soil_column(lst_soil_column)
    # print(df_soil_column.shape)
    # print(df_soil_column[(df_soil_column.lat == 848266) & (df_soil_column.lon == 709018)])
    
    # Create soil type and soil categories
    df_soil_type_categ = soil_type_categ(df_soil_column)
    
    # Create project file
    project_file(df_soil_type_categ, dir_abs,
                 unique_field, unique_value,
                 date_start, date_end)
    
    # Create subsurface map
    subsurface_map(df_soil_type_categ, dir_abs, ND,
                   unique_field, unique_value)

    
#############################################################
# PROCESS
#############################################################
if __name__ == "__main__":
    surface()
    