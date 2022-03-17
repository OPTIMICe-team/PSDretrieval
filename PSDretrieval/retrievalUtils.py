#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Functions used for the retrieval
'''
import numpy as np
import matplotlib.pyplot as plt
from PSDretrieval import scattering as sc
from IPython.terminal.debugger import set_trace
import snowScatt
from snowScatt.instrumentSimulator.radarMoments import Ze


def findBestFittingPartType(selectedParticleTypes, xrSpec,verbose=False,whichDWRsToUse="both"):
    '''
    plot Doppler velocity vs. spectral DWR of X-Ka and Ka-W band combinations (kind of a v-D plot) from the snowScatt simulations
    Arguments:
        INPUT: 
            (model - snowScatt)
            selectedParticleTypes: list of ParticleTypes from which to choose
            (observations)
            xrSpec [can be a single spectra or a time-height window]: xarray containing all spectra information (including DWRs with labels "DWR_X_Ka" and "DWR_Ka_W")
                if single spectra: just plot DV vs. DWRs
                if Window:         average over (vertical wind corrected) DV and plot vs DWR
            (optional)
            verbose: display some progress with print commands
            whichDWRsToUse: select which DWR-vel relation is used
                both: use DWR_X_Ka and DWR_Ka_W
                DWR_X_Ka: use only DWR_X_Ka
                DWR_Ka_W: use only DWR_Ka_W
        OUTPUT: 
            bestPartType: particle type which fits the observation the best
            orderedListPartType: list of particle types ordered by how well they fit the observations
    '''

    print("Start: find best matching particle type in DV-DWR spaces")

    #these two DWR variables are needed here
    DWRkeys = ["DWR_X_Ka","DWR_Ka_W"]
    xrSpecCopy = xrSpec.copy() #this seems necessary for the way the averaging is applied, but there might be a cleaner way

    if xrSpecCopy["KaSpecH"].values.ndim>1:
        print("plot average DV vs DWR for a time-height window")
        for key in DWRkeys:
            xrSpecCopy[key] = xrSpecCopy[key].mean(dim=["time","range"])

    #initialize RMSE variables to find best parameter
    RMSE = dict()
    RMSEall = dict()
    RMSEfinalMin = np.inf 
    for pType in selectedParticleTypes:
        #get particle properties
        Zx, Zk, Zw, Dmax, K2, velModel = sc.model3fOne(pType)

        #calculate spectral DWRs
        DWRmodel = dict()
        DWRmodel["DWR_X_Ka"]    = Zx-Zk
        DWRmodel["DWR_Ka_W"]    = Zk-Zw
        for key in DWRkeys:
            DWRobs = xrSpecCopy[key].copy()
    
            #get an unambiguous DWR array
            DWRmodelUnAmb,ax = sc.getUnambigousDWRdmax(Dmax,DWRmodel[key])
            
            #interpolate the Model DV values to the DWR-grid of the observation
            velModelAtObsDWRgrid = np.interp(DWRobs,DWRmodelUnAmb,velModel)

            #remove DWR values in obs and model, where vel is NaN in the model
            DVobsClean                = -DWRobs.doppler[~np.isnan(velModelAtObsDWRgrid)].values
            velModelAtObsDWRgridClean = velModelAtObsDWRgrid[~np.isnan(velModelAtObsDWRgrid)]

            #calculate the RMSE
            RMSE[key] = np.average((DVobsClean - np.ma.masked_invalid(velModelAtObsDWRgridClean))**2) #this is actually the MSE (mean squared error) and not RMSE, but this is needed in the next lines

        #sum up both RMSE's
        if whichDWRsToUse=="both":
            RMSEfinal         = np.sqrt((RMSE["DWR_X_Ka"] + RMSE["DWR_Ka_W"])/2)
        elif whichDWRsToUse=="DWR_X_Ka":
            RMSEfinal         = np.sqrt(RMSE["DWR_X_Ka"])
        elif whichDWRsToUse=="DWR_Ka_W":
            RMSEfinal         = np.sqrt(RMSE["DWR_Ka_W"])

        RMSEall[pType]  = RMSEfinal

        if RMSEfinal<RMSEfinalMin:
            RMSEfinalMin = RMSEfinal
            bestPartType = pType
        if verbose:
            if whichDWRsToUse=="both":
                print(pType,"MSExk",RMSE["DWR_X_Ka"],"MSEkw",RMSE["DWR_Ka_W"],"RMSExk_kw",RMSEfinal)
            elif whichDWRsToUse=="DWR_X_Ka":
                print(pType,"RMSExk",RMSEfinal)
            elif whichDWRsToUse=="DWR_Ka_W":
                print(pType,"RMSEkw",RMSEfinal)

    orderedListPartType = {k: v for k, v in sorted(RMSEall.items(), key=lambda item: item[1])}.keys()
    if verbose:
        print("best Ptype:",bestPartType,"ordered list",orderedListPartType)
    else:
        print("best Ptype:",bestPartType)

    return bestPartType,orderedListPartType

def dB(x): #conversion: linear [mm**6/m**3] to logarithmic [dB]
    return 10.0*np.log10(x)

def Bd(x): #conversion: logarithmic [dB] to linear [mm**6/m**3]
    return 10.0**(0.1*x)

def removeDropletsFromSpectra(ZkObs,velObs,Nsmooth=10,vSpecLims=[0.0,1.0]):
    '''
    cut the cloud droplet peak from the spectra: 1. smooth spectra 2. find minimum 3. select reasonable minima 4. mask values with DV lower than minimum
    Arguments:
        INPUT & OUTPUT:
            ZkObs   [np.array, dB]:     observed spectral power of the Ka-Band
            velObs  [np.array, m/s]:    velocity from the observation corresponding to ZkObs (positive values are downward, towards the radar)
        INPUT:
            (optional)
            Nsmooth [integer]:          number of bins to smooth before finding minimum
            vSpecLims[float,float]:     Limits where minimum is allowed to be
    '''
    from scipy.signal import argrelextrema
    def moving_average(x, window):
        return np.convolve(x, np.ones(window), 'valid') / window

    # 1. smooth spectra
    ZkObsSmoothed = moving_average(ZkObs,10)
    velObsSmoothed = moving_average(velObs,10)

    # 2. find minimum
    i_minima = argrelextrema(ZkObsSmoothed, np.less)    

    # 3. select reasonable minima
    allMin = velObs[i_minima]
    reasonableMins = list(filter(lambda a: vSpecLims[0]<a<vSpecLims[1], allMin))
    if len(reasonableMins)>1:
        print("Error: more than one minimum (separating the cloud droplet peak and the ice peak) found in interval []" 
          + str(vLowSpecMin) +"," + str(vHighSpecMin) + "]; exit here; you may change the vSpecLims")
    elif len(reasonableMins)==0: #no superdroplet peak found

        velObsWithoutDroplets = velObs.copy()
        ZkObsWithoutDroplets  = ZkObs.copy()

    else: #same as len()=1 -> here we can get rid of the superdroplet peak

        # 4. mask values with DV lower than minimum
        velObsWithoutDroplets = velObs[velObs>reasonableMins].copy()
        ZkObsWithoutDroplets  = ZkObs[velObs>reasonableMins].copy()

    return ZkObsWithoutDroplets,velObsWithoutDroplets

def calculateNumberForEachDVbin(ZkModel,ZkObs,velModel,velObs,DmaxModel=None,removeDroplets=True):
    '''
    calculate number concentration for each doppler velocity (DV) bin
    Arguments:
        INPUT:
            ZkModel [np.array, dB]:     single particle reflectivity of the Ka-Band
            ZkObs   [np.array, dB]:     observed spectral power of the Ka-Band
            velModel[np.array, m/s]:    velocity from the model corresponding to ZkModel (positive values are downward, towards the radar)
            velObs  [np.array, m/s]:    velocity from the observation corresponding to ZkObs (positive values are downward, towards the radar)
            (optional)
            DmaxModel[np.array, m]:     if given this array is also converted to the obs-grid
            removeDroplets[boolean]:    remove the part of the spectra with the supercooled droplets
        OUTPUT: 
            velObs [np.array, m/s]:     part of DV array where number concentration can be retrieved
            NumConNorm [np.array, #/m^3/(m/s)]: Normalized number concentration
            DmaxModelAtObsDVgrid [np.array, m]: maximum dimension corresponding to velObs
    '''

    if removeDroplets:
        ZkObs,velObs = removeDropletsFromSpectra(ZkObs,velObs)

    #convert reflectivities from dB to linear units
    ZkModelLin  = Bd(ZkModel)
    ZkObsLin    = Bd(ZkObs)

    #interpolate the Model Ze values to the DV-grid of the observations
    ZkModelLinAtObsDVgrid = np.interp(velObs,velModel,ZkModelLin)
    if not DmaxModel is None:
        DmaxModelAtObsDVgrid = np.interp(velObs,velModel,DmaxModel)
    else:
        DmaxModelAtObsDVgrid = None

    #calculate the number of particles in each DV bin
    NumCon      = ZkObsLin/ZkModelLinAtObsDVgrid

    #normalize with velocity
    NumConNormV = NumCon[:-1]/-np.diff(velObs)
    NumConNormV[np.isinf(NumConNormV)] = np.nan

    #normalize with Dmax
    NumConNormD = NumCon[:-1]/-np.diff(DmaxModelAtObsDVgrid)
    #NumConNormD[NumCon<1] = np.nan
    NumConNormD[np.isinf(NumConNormD)] = np.nan

    return velObs[:-1],NumConNormV,NumConNormD,DmaxModelAtObsDVgrid[:-1]


def integrateSpectrum(spectrum):
    '''
    Integrate the Doppler Spectrum to get bulk reflectivity
        INPUT:
            Spectrum [np.array, dB]:    spectral power array
        OUTPUT: 
            Ze [float,dB]:               integrated reflectivity
    '''
  
    #convert to linear
    spectrumLin = Bd(spectrum.values)
    #sum up
    ZeLin = np.nansum(spectrumLin)
    #convert to dBz
    Ze  = dB(ZeLin)

    return Ze

def crossCheckIntegratedProp(Dmax,NumConNormD,spectrumX,PartType):
    '''
    Display some cross checks of integrated properties
        INPUT:
            NumConNorm [np.array, #/m^3/(m/s)]:     Normalized number concentration
            Spectrum [np.array, dB]:                spectral power array
            PartType:                               selected particle type 
    '''
   
    #remove nans from Dmax and N(Dmax) 
    DmaxClean           = Dmax[~np.isnan(NumConNormD)]
    NumConNormDclean    = NumConNormD[~np.isnan(NumConNormD)]
    #somehow they are also in the wrong order -> revert
    DmaxClean           = DmaxClean[::-1]
    NumConNormDclean    = NumConNormDclean[::-1]

    ###calculate integrated variables in physical space
    ##from retrieval
    Nretrieval          = np.nansum(NumConNormDclean*np.gradient(DmaxClean)) #total number concentration [#/m^ 3]
    massArray, __,__    = snowScatt.snowMassVelocityArea(DmaxClean, PartType) #particle masses at different sizes [m]
    massDistribution    = NumConNormDclean*massArray #mass distribution [kg/m^3/m]
    IWCretrieval        = np.nansum(massDistribution*np.gradient(DmaxClean)) #[kg/m^ 3] integrate to get IWC 
    massDistCDF         = np.cumsum(massDistribution*np.gradient(DmaxClean))/IWCretrieval #cumulative distribution function of mass [kg/m^3]
    MassMedianDiam      = float(DmaxClean[np.argwhere(massDistCDF>0.5)[0]])

    print("Nretrieval",Nretrieval,"1/m^3","IWCretrieval",IWCretrieval*1e3,"g/m^3","MassMedianDiam",MassMedianDiam*1e3,"mm")

    ##calculate the observed ZeX by integrating the spectrum
    ZeXobs = integrateSpectrum(spectrumX)
    ##estimate ZeX from the retrieved spectra with snowscatt
    wl = snowScatt._compute._c/13.6e9 #X-Band
    ZxFromRetrievedPSD = Ze(DmaxClean, NumConNormDclean, wl, PartType, temperature=273.15) #actual temperature is not considered here

    print("ZeXobs: ",ZeXobs,"ZxFromRetrievedPSD",ZxFromRetrievedPSD)

    return None
