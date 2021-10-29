#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Functions used for the retrieval
'''
import numpy as np
import matplotlib.pyplot as plt

def findNearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def getDmaxFromSDWR(sDWRobs,sDWRmod,Dmax,showIllus=True,ax=None):
    '''
    return a Dmax for DWR(DV)
        sDWRobs: spectral DWR from the spectrum
        sDWRmod: spectral DWR from the scattering database (snowScatt)
        Dmax:    Dmax corresponding to sDWRmod
    
        showIllus: create a plot to illustrate this function
        ax: axes (has to be given if showIllus=True)
    '''

    DmaxfromDWR = np.ones_like(sDWRobs)*np.nan #initialize Dmax(DWR)
    for i,(sDWRnow,DV) in enumerate(zip(sDWRobs.values,sDWRobs.doppler)):
        if sDWRnow<np.nanmax(sDWRmod) and sDWRnow>np.nanmin(sDWRmod): #are we in the unambigious range?
            idx = findNearest(sDWRmod,sDWRnow)
            DmaxfromDWR[i] = Dmax[idx]
        else:
            DmaxfromDWR[i] = np.nan    
        
    if showIllus:
        fig,ax = plt.subplots(nrows=1,ncols=1)
        ax.plot(-sDWRobs.doppler[::-1],DmaxfromDWR[::-1]*1e3)
        plt.xlabel("DV [m/s]")
        plt.ylabel("Dmax [mm]")
    else:
        ax = None

    return DmaxfromDWR

