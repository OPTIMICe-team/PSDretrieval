#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Interface to snowScatt 
'''
import numpy as np
import pandas as pd
import snowScatt

def dB(x): #conversion: linear [mm**6/m**3] to logarithmic [dB]
    return 10.0*np.log10(x)

def Bd(x): #conversion: logarithmic [dB] to linear [mm*6/m**3]
    return 10.0**(0.1*x)

def singlePsd(Ds, i): #monodisperse PSD
    psd = np.zeros_like(Ds)
    psd[i]=1
    return psd

def model3fOne(particleName,Dmax=np.linspace(0.3e-3, 20.0e-3, 2000),lindB="dB",Kfreq=None):
    '''
    get single-particle reflectivity, velocity and dielectric factor

    ARGUMENTS:
        particleName: name of particle types in snowscatt (see snowScatt.snowLibrary.info() for list of all available particle types)
    OPTIONAL:
        Dmax: [m] array of maximum dimensions
        lindB: ["lin","dB"] return Ze either in linear units or in dB
        Kfreq: frequency used to calculate dielectric factor K [None,X,Ka,W]; if None, then the consistent one is used
    RETURNS:
        reflectivity in X-,Ka- and W-Band
        ZxOne, ZkOne, ZwOne, Dmax, K2, vel #ssvel is not wavelength-dependent
    '''

    frequencies =  np.array([9.4e9, 35.6e9, 94.0e9])
    temperature = 270.0

    particle = particleName

    bck = pd.DataFrame(index=Dmax, columns=frequencies)
    for fi, freq in enumerate(frequencies):
        wl = snowScatt._compute._c/freq
        if Kfreq is None:
            eps = snowScatt.refractiveIndex.water.eps(temperature, freq, 'Turner')
        elif Kfreq=="X":
            eps = snowScatt.refractiveIndex.water.eps(temperature, frequencies[0], 'Turner')
        elif Kfreq=="Ka":
            eps = snowScatt.refractiveIndex.water.eps(temperature, frequencies[1], 'Turner')
        elif Kfreq=="W":
            eps = snowScatt.refractiveIndex.water.eps(temperature, frequencies[2], 'Turner')
        K2 = snowScatt.refractiveIndex.utilities.K2(eps)
        #get backscatter and velocity from database
        ssCbck, ssvel = snowScatt.backscatVel(diameters=Dmax,
                                              wavelength=wl,
                                              properties=particle,
                                              temperature=temperature)
        bck[freq] = wl**4*ssCbck/(K2*np.pi**5)

    ZxOne =  np.array([(1.0e18*bck.iloc[:, 0]*singlePsd(Dmax, i)).sum() for i in range(len(Dmax))])
    ZkOne =  np.array([(1.0e18*bck.iloc[:, 1]*singlePsd(Dmax, i)).sum() for i in range(len(Dmax))])
    ZwOne =  np.array([(1.0e18*bck.iloc[:, 2]*singlePsd(Dmax, i)).sum() for i in range(len(Dmax))])
   
    if lindB=="dB": 
        return dB(ZxOne), dB(ZkOne), dB(ZwOne), Dmax, K2, ssvel #ssvel is not wavelength-dependent
    elif lindB=="lin":
        return ZxOne, ZkOne, ZwOne, Dmax, K2, ssvel #ssvel is not wavelength-dependent
    else:
        print("lindB must be either lin or dB")
        sys.exit()

def getDWRs(particleType,Dmax=np.linspace(0.3e-3, 20.0e-3, 2000),Kfreq=None):
    '''
    simple wrapper around model3fOne

    ARGUMENTS:
        particleType: particle type in snowScatt (see snowScatt.snowLibrary.info() to get a list)
    OPTIONAL: 
        Dmax: array of maximum dimensions (does not need to be increasing)
        Kfreq: see model3fOne
    RETURN:
        DWRs: Zx-Zk,Zk-Zw
        Dmax array: Dmax
    '''

    Zx,Zk,Zw,Dmax,__,__    = model3fOne(particleType,Dmax=Dmax,Kfreq=Kfreq)

    return Zx-Zk,Zk-Zw,Dmax

def getUnambigousDWRdmax(Dmax,DWR,DmaxRetr=5e-3,DWRlowDetect=1,showIllus=False,ax=None,verbose=False):
    '''
    get an Unambigous DWR-Dmax relation
        Dmax: [m] maximum dimensions
        DWR:  [dB] dual-wavelength ratio
        (optional)
        DmaxRetr:       #[m] maximum size considered in retrieval; this inexplicitly assumes that larger particles are not relevant
        DWRlowDetect:   #[dB] DWRs smaller than this are disregarded (detection limit)
    
        showIllus: create a plot to illustrate this function
        ax: axes (has to be given if showIllus=True)
        verbose: display some progress with print commands
    RETURN:
        DWRUnamb: range of DWRs which are unambiguous
        ax: axes handle (None if showIllus=False)
    '''

    #mask large sizes
    DWRMaskedLarge = np.ma.masked_where(Dmax>DmaxRetr,DWR)
    DmaxMaskedLarge  = np.ma.masked_where(Dmax>DmaxRetr,Dmax)
    ##find unambiguous DWR-range
    #1. find Dmax where DWR is maximal
    DmaxAtDWRmax = Dmax[np.argmax(DWRMaskedLarge)]
    #2. find min(DWR) in [Dmax(DWRmax)<Dmax<Dmax(DmaxRetr)]
    DWRmin = np.min(DWRMaskedLarge[Dmax>DmaxAtDWRmax])
    #3. find largest unambiguous Dmax
    DmaxLowUnamb = Dmax[DWRMaskedLarge>DWRmin][0]
    if verbose:
        print("DmaxRetr",DmaxRetr*1e3,"DmaxAtDWRmax",DmaxAtDWRmax*1e3,"DWRmin",DWRmin)

    #mask ambiguous sizes
    DWRUnamb = np.ma.masked_where(Dmax>DmaxLowUnamb,DWRMaskedLarge)
    #mask low DWR
    DWRUnamb = np.ma.masked_where(DWRUnamb<DWRlowDetect,DWRUnamb)

    if showIllus:
        ax.semilogx(Dmax*1e3,DWRMaskedLarge,ls="--",c="k")
        ax.semilogx(Dmax*1e3,DWR,ls=":",c="k")
        ax.semilogx(Dmax*1e3,DWRUnamb,c="k")
        ax.axhline(y=DWRmin,c="gray",ls="-")
        ax.axhline(y=DWRlowDetect,c="gray",ls="-")
        ax.set_xlim([5e-1,3e1])
        ax.set_ylim([-1,np.max(DWR)*1.1])
        ax.set_xlabel("Dmax [mm]")
        ax.set_ylabel("sDWR [dB]")
    
    return DWRUnamb,ax

def getSinglePartRefl(particleType,Dmax,freq="k"):
    ''' 
    get the Single particle reflectivity (wrapper around model3fOne)

    ARGUMENTS:
        particleType: particle type in snowScatt (see snowScatt.snowLibrary.info() to get a list)
        Dmax: array of maximum dimensions (does not need to be increasing)
    OPTIONAL: 
        freq: [x,k,w] short name of the frequency
    RETURN:
        ZoneMak0: single particle reflectivity
    '''

    #get the single particle reflectivity
    ZxOne,ZkOne,ZwOne,Dmax,K2,ssvel = model3fOne(particleType,Dmax=Dmax,lindB="lin")
    if freq=="x":
        Zone = ZxOne
    elif freq=="k":
        Zone = ZkOne
    elif freq=="w":
        Zone = ZwOne
    ZoneMask0 = np.ma.masked_where(Zone==0,Zone)

    return ZoneMask0
