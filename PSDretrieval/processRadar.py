#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
reads and processes radar data
'''
import sys

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
        xrSpec = xr.open_dataset(this_dir + "/sample_data/20190122_1455_1000m.nc") #this is window with several times and heights
    else:
        if dataPath is None:
            print("error: provide path to radar data if loadSample=False")
            sys.exit(0)

    return xrSpec

def selectSingleTimeHeight(xrSpec,centered=True,time=None,height=None):
    '''
    select data from a single time-height pixel

    centered:
            True: get data from center in time and range space
            False: provide time and height (TODO: what format??)
    '''


    if centered:
        centeredTime    = xrSpec.time[0] + (xrSpec.time[-1] - xrSpec.time[0])/2.
        centeredHeight  = xrSpec.range[0] + (xrSpec.range[-1] - xrSpec.range[0])/2.
        xrSpecSingle    = xrSpec.sel(time=centeredTime,range=centeredHeight,method="nearest")
    else:
        if (time is None) or (height is None):
            print("time or height cant be None when centered=False")
        else:
            print("TODO: implement time and height selection")

    return xrSpecSingle

def plotObsSpectra(xrSpec,ax):
    '''
    plot observed spectra
    '''

    #plot simple spectrum
    colors      = ["b","g","r"]
    labels      = ["X","Ka","W"]
    
    for i_var,var in enumerate(["XSpecH","KaSpecH","WSpecH"]):
        ax.plot(-xrSpec[var]["doppler"],xrSpec[var],color=colors[i_var],label=labels[i_var],lw=1)
        #if not var=="WSpecH":
        #    ax.axhline(y=xrSpec[var + "specNoise"]+delNoiseLevel,lw=0.1,ls="--",color=colors[i_var])

        #ax0.text(0.95,0.05, "h={}m\nt={}".format(height,dateSpec),horizontalalignment='right',verticalalignment='bottom',transform=ax0.transAxes)
        ax.legend()
        ax.set_xlabel("v [m/s]")
        ax.set_ylabel("z [dBz]")
        ax.set_xlim([-0.5,2])

    return ax
