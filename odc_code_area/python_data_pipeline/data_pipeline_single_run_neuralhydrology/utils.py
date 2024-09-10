from pathlib import Path
import glob
import numpy as np
import os

# Modules for config.ini
from configparser import ConfigParser
config = ConfigParser()

# Read config file values
config.read('config.ini')

def folder_exists(folder_path):
    return Path(folder_path).is_dir()

def create_folder(folder_path):
    try:
        if not folder_exists(folder_path):
            Path(folder_path).mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"Error: {e.strerror}")

def delete_folder(folder_path):
    try:
        if folder_exists(folder_path):
            input(f'This will delete files from {folder_path} are you sure?')
            for file_path in Path(folder_path).glob(f'*.csv'):
                try:
                    file_path.unlink()
                except OSError as e:
                    print("Error: %s : %s" % (file_path, e.strerror))
            # Path(folder_path).rmdir()
            
    except OSError as e:
        print(f"Error: {e.strerror}")

def file_remove(dir_path, single_extension):
    # dir_abs = Path().resolve()
        
    # file_name = file_check_path.stem
    # file_path = file_check_path.parent
    lst_extensions = ['cpg', 'dbf', 'prj', 'shp', 'shx', 'tif', 'csv']
    ext_flag = single_extension
    
    # delete all files
    if ext_flag == 'all':
        for extension in lst_extensions:
            temp_lst = []
            temp_lst = list(dir_path.glob(f"*.{extension}"))
            
            if len(temp_lst) != 0:
                for f in temp_lst:
                    if f.is_file():
                        f.unlink()
                        print(f'{f.name} DELETED!')
    # delete selected extension
    else:
        temp_lst = list(dir_path.glob(f"*.{ext_flag}"))
        if len(temp_lst) != 0:
            for f in temp_lst:
                if f.is_file():
                    f.unlink()
                    print(f'{f.name} DELETED!')
    

    # for extension in lst_extensions:
    #     file_other_name= f'{file_name}.{extension}'
    #     file_to_check = Path(dir_abs / Path(file_path / file_other_name))
        
    #     if file_to_check.is_file():
    #         file_to_check.unlink()
    #         print(f'{file_other_name} DELETED!')

def sync_folder(folder_path):
    delete_folder(folder_path)
    create_folder(folder_path)
    
def shetran_csv_file(main_dir, file_name, pivot_dataframe, specifier):
    # Saving dataframe as text file
    filename = Path(main_dir / f'model_input/{file_name}.txt')
    np.savetxt(filename, pivot_dataframe.values, fmt=f'%{specifier}')

    # reating header needed for SHETRAN
    ncols = pivot_dataframe.shape[1]
    nrows = pivot_dataframe.shape[0]
    xllcorner = int(list(pivot_dataframe.columns)[0])
    yllcorner = int(pivot_dataframe.index[-1])
    cellsize = config.getint('res_setting', 'width')
    no_data_val = config.getint('res_setting', 'no_data')

    # Adding header to text file
    # Copying current information in text file
    append_copy = open(filename, "r")
    original_text = append_copy.read()
    append_copy.close()
    # Adding header information - this deletes any information in the text file
    append_copy = open(filename, "w")
    append_copy.write(
        "ncols         " + str(ncols) + "\n" + 
        "nrows         " + str(nrows) +  "\n" +
        "xllcorner     " + str(xllcorner) +  "\n" +
        "yllcorner     " + str(yllcorner) + "\n" +
        "cellsize      " + str(cellsize) + "\n" +
        "NODATA_value  " + str(no_data_val) + "\n")
    # Pasting content that was in text file before adding header
    append_copy.write(original_text)
    # Saving text file
    append_copy.close()

def get_local_folder():
    return os.path.dirname(os.path.realpath(__file__))
