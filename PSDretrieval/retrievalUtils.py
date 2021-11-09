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

def getDmaxFromSDWR(sDWRobs,sDWRmod,Dmax,showIllus=False,ax=None):
    '''
    return a Dmax for DWR(DV)

    ARGUMENTS:
        sDWRobs: spectral DWR from the spectrum (must be xarray for plotting; can also be np.darray for calculation only)
        sDWRmod: spectral DWR from the scattering database (snowScatt)
        Dmax:    Dmax corresponding to sDWRmod
    OPTIONAL: 
        showIllus: create a plot to illustrate this function
        ax: axes (has to be given if showIllus=True)
    REURN:
        Dmax corresponding to sDWRobs
    '''

    DmaxfromDWR = np.ones_like(sDWRobs)*np.nan #initialize Dmax(DWR)
    if not isinstance(sDWRobs, np.ndarray):
        try:
            sDWRobsvals = sDWRobs.values
        except: 
            print("sDWRobs in retrievalUtils.getDmaxFromSDWR must be either a np.ndarray of xarray variable")
            sys.exit(0)
    else:
        sDWRobsvals = sDWRobs

    for i,sDWRnow in enumerate(sDWRobsvals):
        if sDWRnow<np.nanmax(sDWRmod) and sDWRnow>np.nanmin(sDWRmod): #are we in the unambigious range?
            idx = findNearest(sDWRmod,sDWRnow)
            DmaxfromDWR[i] = Dmax[idx]
        else:
            DmaxfromDWR[i] = np.nan    
        
    if showIllus:
        if isinstance(sDWRobs, np.ndarray):
            print("plotting currently only possible with xarray -> skip plotting")
        else:
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


def func(x, a, b):
    return a * x ** b

def fitting2D(xx, yy,p0=[8.57,0.393]): #,p0):
    from scipy.optimize import curve_fit
    xx[np.isnan(yy)] = np.nan
    yy[np.isnan(xx)] = np.nan
    xx = xx[np.isfinite(xx)]
    yy = yy[np.isfinite(yy)]
    popt, pcov = curve_fit(func, xx, yy) #,p0[0],p0[1])
    [coeff_a,coeff_b] = popt[0],popt[1]
    #[coeff_a,coeff_b] = p0[0],p0[1]
    return coeff_a,coeff_b

def histDWRandDmaxVsDv(xrDWR,Spec,SpecNoise,DWRUnamb,Dmax,aboveNoiseThreshold=30,showIllus=False,ax=None,fig=None):
    '''
        plot DWR and Dmax vs Doppler velocity
        and fit a power-law to vterm(Dmax)
    ARGUMENTS:
        xrDWR: x-array variable which contains the DWR values and the doppler velocity (xrDWR.doppler)
        Spec: spectral power (carefully choose the frequency - probably the one with the lowest sensitivity is the best)
        SpecNoise: Noise level (same frequency as Spec!)
        DWRUnamb: modeled unambigious DWR array
        Dmax: maximum dimension array corresponding to DWRUnamb
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

        hist,xedge,yedge,im = ax[0].hist2d(x_hist,y_hist,bins=30,range=[[-0.5,2],[-2,17]],cmap=my_cmap,vmin=1) #somehow both: "set_under" and "vmin" is necessary
        cbar = fig.colorbar(im,ax=ax[0])
        cbar.set_label("counts")
        ax[0].set_xlabel("DV [m/s]")
        ax[0].set_ylabel("sDWR [dB]")
    else:
        
        hist,xedge,yedge = np.histogram2d(x_hist,y_hist,bins=30,range=[[-0.5,2],[-2,17]])

    yedgeDmax = getDmaxFromSDWR(yedge,DWRUnamb,Dmax) #calculate Dmax corresponding to the yedge-DWR array
    yhistDmax = getDmaxFromSDWR(y_hist,DWRUnamb,Dmax) #calculate Dmax corresponding to the y_hist-DWR array

    if showIllus:
        yedgeDmaxMatrix             = np.reshape(np.tile(yedgeDmax,xedge.shape[0]),(xedge.shape[0],yedge.shape[0])) #create array with dimension [xedge,yedge]
        yedgeDmaxValid              = yedgeDmax[yedgeDmax>0] #remove nans from yedgeDmax
        #yedgeDmaxMatrixValidOnly    = np.reshape(np.tile(yedgeDmaxValid,xedge.shape[0]),(xedge.shape[0],yedgeDmaxValid.shape[0])) #same as yedgeDmaxMatrix but only with valid values

        histValidDmax = hist[:,np.where(~np.isnan(yedgeDmax))[0]] #hist[~np.isnan(yedgeDmaxMatrix[:-1,:-1])]
        #histMaskedDmaxValid = np.ma.masked_where(np.isnan(yedgeDmaxMatrix[:-1,:-1]),hist)
        #histMaskedDmaxValid = histMaskedDmaxValid[histMaskedDmaxValid>0]
        print(xedge.shape,yedgeDmaxValid.shape,histValidDmax.shape)
        ax[1].pcolor(xedge[:-1],yedgeDmaxValid*1e3,np.transpose(histValidDmax),cmap=my_cmap,vmin=1)
        cbar.set_label("counts")
        ax[1].set_xlabel("DV [m/s]")
        ax[1].set_ylabel("D$_{max}$ [mm]")
        #ax[1].set_ylim([1e3*min(yedgeDmaxValid),1e3*max(yedgeDmaxValid)]) #TODO: why doesnt this work??
        ax[1].set_ylim(bottom=1e3*min(yedgeDmaxValid),top=12.)
        #ax[1].set_yscale('log')
        plt.tight_layout()

        [a_fit,b_fit] = fitting2D(yhistDmax,x_hist)
        print(a_fit,b_fit)
        #ax[1].plot(xedge[:-1],1./a_fit*xedge[:-1]**(1./b_fit))
        DmaxFitEval = np.linspace(1e-4,1e-2,100) #simple Dmax-array to evaluate the fit 
        ax[2].plot(DmaxFitEval*1e3,a_fit*DmaxFitEval**b_fit)
        ax[2].set_xlabel("D$_{max}$ [mm]")
        ax[2].set_ylabel("v$_{term}$ [m/s]")
        #import pdb; pdb.set_trace()

    return ax,hist,xedge,yedge
