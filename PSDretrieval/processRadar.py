#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
reads and processes radar data
'''
import sys
import numpy as np
import pandas as pd
import xarray as xr
import glob
from IPython.core.debugger import set_trace


def getPandasTime(date=None,time=None):
    '''
    get pandas time format from strings of date and time
    Arguments:
        INPUT:    
            date: date in format yyyymmdd
            time: time in format HH:MM
        OUTPUT:
            tStep: time step in pandas format
    '''

    tStep = pd.to_datetime(date+' ' + time)

    return tStep

def regridSpec(data,newVelRes=0.01,windowWidth=10,vRange=[-5,1]):
    '''
    regrid the spectra to a common doppler velocity grid
    The regridded spectra has typically a finer, but smoothed grid compared to the input spectra.
    Arguments:
        INPUT:    
            data: data from all frequencies
        optional
            newVelRes [m/s]: velocity resolution of the regridded spectra
            windowWidth [bins]: width of moving average window
            vRange [m/s]: velocity range of regridded spectra
        OUTPUT:
            dataOut: regridded spectra
    '''
    #xarray doesnt allow reindexing and rename the doppler velocity coordinate inside of the original dataset -> create new dataset
    dataOut = xr.Dataset()

    #new resolution and doppler velocity vector
    newVel = np.arange(vRange[0], vRange[1], newVelRes) 
    #normalize spectrum
    dvX = np.abs(np.diff(data['dopplerX'].values)[0])
    X = data['XSpecH']/dvX
    # interpolate to new velocity
    XSpecInt = X.interp(dopplerX=newVel) 
    XSpecInt = XSpecInt*newVelRes
    # rename doppler velocity grid
    XSpecInt = XSpecInt.rename({'dopplerX':'doppler'})
    print('interp X done')

    #normalization, interpolation and renaming for Ka-Band same as for X-Band above (TODO: wouldnt that be better in a for-loop?)
    dvKa = np.abs(np.diff(data['dopplerKa'].values)[0])
    Ka = data['KaSpecH']/dvKa
    KaSpecInt = Ka.interp(dopplerKa=newVel)
    KaSpecInt = KaSpecInt*newVelRes
    KaSpecInt = KaSpecInt.rename({'dopplerKa':'doppler'})
    print('interp Ka done')

    #normalization, interpolation and renaming for W-Band same as for X-Band above (TODO: wouldnt that be better in a for-loop?)
    dvW = np.abs(np.diff(data['dopplerW'].values)[0])
    W = data['WSpecH']/dvW
    WSpecInt = W.interp(dopplerW=newVel)
    WSpecInt = WSpecInt*newVelRes
    WSpecInt = WSpecInt.rename({'dopplerW':'doppler'})
    print('interp W done')A

    ##correct noise: make true linear (saved are linear units to which the log-lin transformation was applied once to much)
    data["KaSpecNoiseH"] = 10.*np.log10(10.*np.log10(data["KaSpecNoiseH"]))
    data["XSpecNoiseH"] = 10.*np.log10(10.*np.log10(data["XSpecNoiseH"]))

    ##overwrite dataframe
    dataOut = xr.merge([dataOut,WSpecInt,KaSpecInt,XSpecInt,data["XSpecNoiseH"],data["KaSpecNoiseH"]],compat='override')
    print('merging datasets done')
    dataOut = dataOut.rolling(doppler=windowWidth,min_periods=1,center=True).mean() #smooth dataset

    return dataOut

def addOffsets(data,data2,test_interp=False):
    '''
    add offsets and atmospheric attenuation stored in the level 2 data to the level 0 data
    Arguments:
        INPUT:    
            data2: level 2 data containing the offset information
        INPUT and OUTPUT:
            data: level 0 data containing the spectral data
        optional
            test_interp: test interpolation against level 2 data
    '''
    data['XSpecH'] = 10*np.log10(data['XSpecH']) + data2.rain_offset_X + data2.offset_x + data2.pia_x
    print('X offsets added')
    data['KaSpecH'] = 10*np.log10(data['KaSpecH']) + data2.rain_offset_Ka +  data2.pia_ka
    data['WSpecH'] = 10*np.log10(data['WSpecH']) + data2.rain_offset_W + data2.offset_w + data2.pia_w

    if test_interp==True:
        data['linKaSpec'] = 10**(data['KaSpecH']/10)
        data['ZeKa'] = data['linKaSpec'].sum(dim='doppler')
        data['Ka_DBZ'] = 10*np.log10(data['ZeKa'])

        data['linWSpec'] = 10**(data['WSpecH']/10)
        data['ZeW'] = data['linWSpec'].sum(dim='doppler')
        data['W_DBZ'] = 10*np.log10(data['ZeW'])
        #dataDWR['W_DBZ_LV0'] = 10*np.log10(data['W_Z_H'])+offsetW+data2.pia_w
        #data1['W_DBZ_H'] = data1['W_DBZ_H']+offsetW+data2.pia_w
        data['linXSpec'] = 10**(data['XSpecH']/10)
        data['ZeX'] = data['linXSpec'].sum(dim='doppler')
        data['X_DBZ'] = 10*np.log10(data['ZeX'])

    data['DWR_X_Ka'] = data['XSpecH'] - data['KaSpecH']
    data['DWR_Ka_W'] = data['KaSpecH'] - data['WSpecH']

    return data

def loadTripexPol(dataPath="/data/obs/campaigns/tripex-pol/processed/",date="20190113",time="06:18:04",tRange=0,hRange=180,hcenter=1000.):
    '''
    load spectra (+ regrid, offset-correction, ...) from level 0 data of the Tripex-pol dataset 
    Arguments:
        optional INPUT:    
            dataPath:  path to the level 0 and level 2 data
            date [yyyymmdd]:    date
            time [HH:MM:SS]:    time
            tRange [Delta min]: time range
            hRange [Delta m]:   height range
            hcenter [m]:        center of height range
        OUTPUT:
            xrSpec: spectra which has been read in, regridded, corrected for offsets + information on pressure and noise
    '''
  
    #convert strings to pandas time format 
    tStep = getPandasTime(date=date,time=time)

    ####load LEVEL1 data
    #get list of files
    dataLV0List = glob.glob(dataPath + "tripex_pol_level_0/" +tStep.strftime('%Y')+'/'+tStep.strftime('%m')+'/'+tStep.strftime('%d')+'/'
                                    +tStep.strftime('%Y')+    tStep.strftime('%m')+    tStep.strftime('%d')+'_*' 
                            '_tripex_pol_3fr_spec_filtered_regridded.nc')
    #load data
    dataLV0 = xr.open_mfdataset(dataLV0List)
    #change time attribute
    dataLV0.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    #decode according to cf-conventions
    dataLV0 = xr.decode_cf(dataLV0) 
    #select time step
    xrSpec = dataLV0.sel(time=slice(tStep-pd.Timedelta(minutes=tRange/2.),tStep+pd.Timedelta(minutes=tRange/2.)),range=slice(hcenter-hRange/2,hcenter+hRange/2))
    xrSpec = regridSpec(xrSpec,windowWidth=10)

    ####load LEVEL2 data
    dataLV2 = xr.open_mfdataset(dataPath + "tripex_pol_level_2/" + date + '_tripex_pol_3fr_L2_mom.nc')
    #select time step
    dataLV2 = dataLV2.sel(time=slice(tStep-pd.Timedelta(minutes=tRange/2.),tStep+pd.Timedelta(minutes=tRange/2.)),range=slice(hcenter-hRange/2,hcenter+hRange/2))
        
    #get offsets from LV2 data
    xrSpec       = addOffsets(xrSpec,dataLV2)
    #add pressure to the file (needed for density correction of fall speed)
    xrSpec["pa"] = dataLV2["pa"] 
    
    return xrSpec

def loadSpectra(loadSample=True,dataPath=None,createSample=False,date="20190113",time="06:18:04",tRange=1,hRange=180,hcenter=1600):
    '''
    load Doppler spectra data from 1) any dataset (requires some additional implementions-  maybe in a git-branch) 2) the Tripex-pol dataset  3) the sample_data folder
    Arguments:
        INPUT:    
            loadSample [boolean]: 
                if True: load processed data from sample_data directory
                if False: load data from path
            dataPath: path to the spectra files
            createSample:  Create a sample from the data, so only the subset of the data must be read in later (saves time, especially for developing the retrieval)
            tRange [Delta min]:     time range to read
            hRange [Delta m]:       height range to read
            hcenter[m]:             center of height range
        OUTPUT:
            xrSpec: xarray-dataset containing the spectra data
    '''
    import xarray as xr
    from os import path
   
    #path of this script 
    this_dir = path.dirname(path.realpath(__file__))
    #define sample path (needed for both saving and reading a sample)
    sample_path = this_dir + "/sample_data/" + date + "_" + time[0:2] + time[3:5] + time[6:8] + "_tR" + str(tRange) + "min_h" + str(hcenter) + "m_hR" + str(hRange) + "m.nc"

    if loadSample:
        if createSample:
            print("ERROR: cant create sample when loading sample: set either loadSample or createSample to False")
            sys.exit(0)
        if path.exists(sample_path):
            #load a window with several times and heights
            xrSpec = xr.open_dataset(sample_path)
        else:
            fileListForDay = glob.glob(sample_path.split("_tR")[0] + "*")
            if len(fileListForDay)==0:
                print("ERROR: no data for date: " + date + " and time" + time + "in PSDretrievalPSDretrieval//sample_data yet")
                sys.exit(0)
            else:
                print("ERROR: check time Range (tRange", tRange, "), height (hcenter", hcenter, ") and height-Range (hRange", hRange, ") in call of loadSpectra()\n"
                        "files available are: ", fileListForDay)
                sys.exit(0)
    else:
        if dataPath is None:
            print("error: provide path to radar data if loadSample=False")
            sys.exit(0)
        else: 
            '''
            load a spectra file
            '''
            if "tripex-pol" in dataPath:
                print("load file:"  + dataPath)
                xrSpec = loadTripexPol(dataPath="/data/obs/campaigns/tripex-pol/processed/",date=date,time=time,tRange=tRange,hcenter=hcenter,hRange=hRange)
                if createSample:
                    xrSpec.to_netcdf(sample_path)
            else:
                print("error: the path is not found or the file is not in the right format")
                sys.exit(0)
            

    return xrSpec

def selectSingleTimeHeight(xrSpec,centered=True,pdTime=None,height=None):
    '''
    select data from a single time and height
    Arguments:
        INPUT:    
            xrSpec: xarray-dataset containing the spectra data
        optional
            centered [boolean]:
                True: get data from center in time and range space
                False: provide time and height (TODO: not implemented!)
            pDTime:     pandas time stamp
            hcenter:    height range
        OUTPUT:
            xrSpecSingle: Spectrum from a single time and height
    '''

    if centered:
        centeredTime    = xrSpec.time[0] + (xrSpec.time[-1] - xrSpec.time[0])/2.
        centeredHeight  = xrSpec.range[0] + (xrSpec.range[-1] - xrSpec.range[0])/2.
        xrSpecSingle    = xrSpec.sel(time=centeredTime,range=centeredHeight,method="nearest")
    else:
        if (time is None) or (height is None):
            print("ERROR: time or height cant be None when centered=False in selectSingleTimeHeight")
            sys.exit(0)
        else:
            print("ERROR: time-height selection not implemented in selectSingleTimeHeight")
            sys.exit(0)

    return xrSpecSingle

def shiftSpectraByVerticalWind(SpecWindow):
    '''
    Shift the spectra by the vertical wind diagnosed with the position of the cloud droplet peak
    Arguments:
        INPUT:    
        IN- & OUTPUT:
            SpecWindow: spectra from a time-height window
            wArray:     array (time-height window) with estimated vertical wind 
    '''

    return SpecWindow,wArray
