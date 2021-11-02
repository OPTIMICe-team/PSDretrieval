#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Functions used for the retrieval
'''
import numpy as np
import matplotlib.pyplot as plt
from PSDretrieval import scattering as sc

def findNearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def getDmaxFromSDWR(sDWRobs,sDWRmod,Dmax,showIllus=True,ax=None):
    '''
    return a Dmax for DWR(DV)

    ARGUMENTS:
        sDWRobs: spectral DWR from the spectrum
        sDWRmod: spectral DWR from the scattering database (snowScatt)
        Dmax:    Dmax corresponding to sDWRmod
    OPTIONAL: 
        showIllus: create a plot to illustrate this function
        ax: axes (has to be given if showIllus=True)
    REURN:
        Dmax corresponding to sDWRobs
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

    return DmaxfromDWR

def calcNumberConcFromSpectrumAndZOne(Zobs,Zone,showIllus=False,ax=None):
    '''
    calculate the number concentration from the spectrum and the single particle scattering properties at the given Doppler velocity bin
    ARGUMENTS:
        Zobs: xarray vector with observed reflectivity and doppler velocity in Zobs.doppler
        Zone: [mm^6/m^3] single particle reflectivity (must be from same frequency as Zobs!!)
    OPTIONAL:
        showIllus: create a plot to illustrate this function
        ax: axes (has to be given if showIllus=True)
    RETURN:
        N: [m^-3/(m/s)] number concentration normalized with Doppler velocity 
        ax: axis handle
    '''

    ## divide the observed z (power at a specific DV) by Zone to get the number concentration N
    N = sc.Bd(Zobs)/ Zone
    Nnorm = N/(Zobs.doppler[1] - Zobs.doppler[0]) 

    if showIllus:
        ax.semilogy(-Zobs.doppler,Nnorm)
        ax.set_xlabel("DV [m/s]")
        ax.set_ylabel("N [m$^{-3}$/(m/s)]")

    return Nnorm,ax

def histDWRandDmaxVsDv(xrDWR,Spec,SpecNoise,aboveNoiseThreshold=30,showIllus=False,ax=None,fig=None):
    '''
        plot DWR and Dmax vs Doppler velocity
        and fit a power-law to vterm(Dmax)
    ARGUMENTS:
        xrDWR: x-array variable which contains the DWR values and the doppler velocity (xrDWR.doppler)
        Spec: spectral power (carefully choose the frequency - probably the one with the lowest sensitivity is the best)
        SpecNoise: Noise level (same frequency as Spec!)
    OPTIONAL:
        aboveNoiseThreshold: dismiss all DV sections where the reflectivity is less than 30dBz above the noise threshold
        showIllus: create a plot to illustrate this function
        ax: axes (has to be given if showIllus=True)
    RETURN:
        ax: axes handle
    '''

    #get the DV array #np.tile() flattens it and the "-" in front changes the convention so that positive DV is upward
    x_hist = -np.tile(xrDWR.doppler.values,xrDWR.values.shape[0]*xrDWR.values.shape[1])
    #flatten the dwr-array to be able to make a histogram
    y_hist = xrDWR.values.flatten(order="C")
    #remove the data below the noise level
    KaSpecFlat = Spec.values.flatten(order="F")
    KaGTnoiseFlag = (Spec>(SpecNoise+aboveNoiseThreshold)).values.flatten()
    y_hist[~KaGTnoiseFlag] = np.nan

    if showIllus:
        import matplotlib as mpl
        import copy
        my_cmap = copy.copy(mpl.cm.get_cmap('jet')) # copy the default cmap
        my_cmap.set_under("w",1)
        hist,xedge,yedge,im = ax.hist2d(x_hist,y_hist,bins=30,range=[[-0.5,2],[-2,17]],cmap=my_cmap,vmin=1) #somehow both: "set_under" and "vmin" is necessary
        cbar = fig.colorbar(im,ax=ax)
        cbar.set_label("counts")
        ax.set_xlabel("DV [m/s]")
        ax.set_ylabel("sDWR$_{Ka,W}$ [kg]")

    return ax
