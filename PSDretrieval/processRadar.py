#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
reads and processes radar data
'''

def loadSpectra(loadSample=True,dataPath=None):
    '''
    load Doppler spectra data

    loadSample: 
        if True: load processed data from sample_data directory
        if False: load data from path
    '''
    import xarray as xr
    from os import path
    
    this_dir = path.dirname(path.realpath(__file__)) #path of this script

    if loadSample:
        xrSpec = xr.open_dataset(this_dir + "/sample_data/20190130_1430.nc") #this is window with several times and heights
    else:
        if dataPath is None:
            print("error: provide path to radar data if loadSample=False")
            sys.exit(0)

    return xrSpec

def 
