# -*- coding: utf-8 -*-
"""
A few functions for loading data from files to be plotted.
    load_data(): coordinates and map data from formatted *.txt or *.csv file
    load_ibexfile(): load data from a file from IBEX data release
    get_ibex_filepath(): determine the directory path to specific IBEX data file
    load_json(): load *.json configuration file

Version: v1.3 (2025-10-07)
Created on Tue Mar 12 10:28:27 2024

@author: jgazer (Jonathan Gasser, SwRI)
"""

import numpy as np

def load_data(filename, path=None):
    '''
    Load map data from a *.txt or *.csv file that is formatted as follows:
        * Comments: All lines in the file that start with '#' are disregarded.
        * Lon: the first non-commented line should contain n+1 comma-separated float values in ascending order,
            representing the longitude pixel edges in [deg].
        * Lat: the 2nd non-commented line should contain m+1 comma-separated float values in ascending order,
            representing the latitude pixel edges in [deg].
        * Data: the file should contain at least additional m lines with n comma-separated float values each,
            representing the mxn pixel values (in any unit).
    
    Args:
        filename (str): the filename incl file extension of the data file.
        path (str): the directory path of the folder containing the data file.
    
    Returns:
        lon (numpy.array): n+1 array containing the longitude coordinates
        lat (numpy.array): m+1 array containing the latitude coordinates
        data (numpy.array): m x n array containing the map pixel values
    '''
    if path is None:
        path = load_json()['path_data']
    if not path.endswith('/'):
        path +='/'
    if filename.endswith('.txt'):
        filename = filename.replace('.txt','.csv')
    elif not filename.endswith('.csv'):
        filename += '.csv'
    fid = open(path +filename, 'r')
    
    lon, lat, buf = [],[],[]
    for line in fid:
        if line[0] == '#':
            continue
        elif lon ==[]:
            lon= [ float(s) for s in line.replace('\n','').split(',')]
        elif lat ==[]:
            lat= [ float(s) for s in line.replace('\n','').split(',')]
        else:
            nums= [ float(s) for s in line.replace('\n','').split(',') if s]
            buf.append(nums)
    fid.close()
    
    data = np.array(buf)
    return np.array(lon), np.array(lat), data
# END OF FUNCTION

def load_ibexfile(filename, path=None):
    """
    Load data file formatted as in IBEX data releases, based on filename and path.
    IBEX Data Releases can be downloaded from the IBEX website at https://ibex.princeton.edu/DataRelease
    This function does NOT modify data to correct for
        * zero values that should be interpreted as NAN (in flux, fvar, ..)
        * inconsistent NAN entries between different quantities
    To load data files from a locally stored Data Release directory, 
        based on desired data parameters, execute:
        > load_ibexfile( *get_ibex_filepath )
    
    Args:
        filename (str): the subfolder and filename of datafile.
        path (str): the data containing directory path.
            If None, the correct path for Lo or Hi data is taken from JSON file.
    Returns:
        lon (numpy.array): Ecl. Longitude 1d array
        lat (numpy.array): Ecl. Latitude 1d array
        data (numpy.array): 2d array with data from file.
    """
    if path is None:
        if '-hi-' in filename:
            path = load_json()['path_ibex_hi']
        elif '-lo-' in filename:
            path = load_json()['path_ibex_lo']
        else:
            path ='' #runtime path
    
    if path and path[-1] != '/':
        path +='/'
    if filename.endswith('.csv'):
        filename = filename.replace('.csv','.txt')
    elif not filename.endswith('.txt'):
        filename += '.txt'
    fid = open(path +filename, 'r')
    
    hdline= fid.readline()
    lat_num= int(hdline[5:7])
    lon_num= int(hdline[8:10])
    
    buf =[]
    for line in fid:
        if line[0] == '#':
            continue
        nums= [ float(s) for s in line.split(' ') if s]
        buf.append(nums)
    fid.close()
    
    data = np.array(buf) # flux data set from file
    lon= np.linspace(0, 360, lon_num+1) # ecl.longitude bin edges
    lat= np.linspace(-90, 90, lat_num+1) # ecl.latitude bin edges

    return lon, lat, data
# END OF FUNCTION

def load_ibexfile_dr(**kwargs):
    """ Wrapper for  'load_ibexfile( *get_ibex_filepath()' ).
    See docstring of 'get_ibex_filepath()' and 'load_ibexfile()'.
    """
    return load_ibexfile( get_ibex_filepath(**kwargs) )
                         
def get_ibex_filepath( path0=None, is_hi=True, dataprod='noSP_ram', subf_ext='',
                     ebin=3, year=2009, qty='flux'):
    """
    This function constructs the path of data files within a local IBEX data release directory,
        based on desired data parameters (data product, quantity, energy bin, ..)

    Args:
        path0 (str): path of the data folder containing the different dataprod subfolders.
            for DR17 (IBEX-Lo): '~/DR_17/Lo-Hydrogen/'
            for DR18 (IBEX-Hi): '~/DR_18/Hi/'
        is_hi (bool): if True, gets the path for IBEX-Hi data, else IBEX-Lo data.
        dataprod (str): specifies the type of data (e.g., CG? SP? ram/antiram?)
        subf_ext (str): extended subfolders below the dataprod level (e.g., 'prod/map_h/')
        ebin (str|int): energy bin number
        year (str|int): year of the data. For some dataproducts of Hi, must specify 'A' or 'B'.
        qty (str): data quantity to be loaded (e.g., 'flux', 'fvar', 'fexp',..)
        
    Returns:
        filename (str): data file name.
        filepath (str): directory path to the folder containing the data file.
    """
    if is_hi:
        if path0 is None:
            path0 = load_json()['path_ibex_hi']
        subpath = 'hvset_'+dataprod +'_'+str(year) +'/' +subf_ext
        filename = 'hv60.hide-trp-flux100-hi-'+str(ebin) +'-'+qty+'.txt'
    else:
        if path0 is None:
            path0 = load_json()['path_ibex_lo']
        subpath = 'lvset_h_'+dataprod +'_hb_'+str(year) +'/' +subf_ext
        filename = 'lv60.lohb-trp-flux100-lo-'+str(ebin) +'-'+qty+'.txt'
    if not path0.endswith('/'):
        path0 +='/'
    filepath = path0 +subpath
    return filename, filepath
# END OF FUNCTION

def load_json(filename='config'):
    """
    Load JSON file containing the default config variables (path,...).
    
    Args:
        filename (str): relative path and file name of json file.
    Returns:
        JSON object of the loaded .json file
    """
    import json
    rout = __file__.replace('\\','/')
    this_path = rout[0: 1+rout.rindex('/')]
    if filename.endswith('.json'):
        filename = filename[:-5]
    with open(this_path +filename +'.json', 'r') as fid:
        jsn = json.load(fid)
    return jsn
