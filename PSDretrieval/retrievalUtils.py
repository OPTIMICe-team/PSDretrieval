#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Functions used for the retrieval
'''
import numpy as np
import matplotlib.pyplot as plt
from PSDretrieval import scattering as sc
from IPython.terminal.debugger import set_trace

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
        DWRmodel["DWR_X_Ka"]  = Zx-Zk
        DWRmodel["DWR_Ka_W"] = Zk-Zw
        for key in DWRkeys:
            DWRobs = xrSpecCopy[key].copy()
            
            #interpolate the Model DV values to the DWR-grid of the observation
            velModelAtObsDWRgrid = np.interp(DWRobs,DWRmodel[key],velModel)

            #remove DWR values in obs and model, where vel is NaN in the model
            DWRObsClean               = DWRobs.doppler[~np.isnan(velModelAtObsDWRgrid)]
            velModelAtObsDWRgridClean = velModelAtObsDWRgrid[~np.isnan(velModelAtObsDWRgrid)]

            #calculate the RMSE
            RMSE[key] = np.average((DWRObsClean.doppler - np.ma.masked_invalid(velModelAtObsDWRgridClean))**2) #this is actually the MSE (mean squared error) and not RMSE, but this is needed in the next lines

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
            
                print(pType,"MSExk",RMSEfinal)
            elif whichDWRsToUse=="DWR_Ka_W":
                print(pType,"MSEkw",RMSEfinal)

    orderedListPartType = {k: v for k, v in sorted(RMSEall.items(), key=lambda item: item[1])}.keys()
    if verbose:
        print("best Ptype:",bestPartType,"ordered list",orderedListPartType)
    else:
        print("best Ptype:",bestPartType)

    return bestPartType,orderedListPartType
