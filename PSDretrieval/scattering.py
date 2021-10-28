#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Interface to snowScatt 
'''
import numpy as np
import pandas as pd
import snowScatt

def dB(x):
    return 10.0*np.log10(x)

def Bd(x):
    return 10.0**(0.1*x)

def singlePsd(Ds, i):
    psd = np.zeros_like(Ds)
    psd[i]=1
    return psd

def model3fOne(particleName):
    '''
    get single-particle reflectivity, velocity and dielectric factor

    particleName: name of particle types in snowscatt (see snowScatt.snowLibrary.info() for list of all available particle types)
    '''
    frequencies =  np.array([9.4e9, 35.6e9, 94.0e9])
    temperature = 270.0

    Dmax = np.linspace(0.3e-3, 20.0e-3, 2000) # list of sizes
    particle = particleName

    bck = pd.DataFrame(index=Dmax, columns=frequencies)
    for fi, freq in enumerate(frequencies):
        wl = snowScatt._compute._c/freq
        eps = snowScatt.refractiveIndex.water.eps(temperature, freq, 'Turner')
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
    
    return dB(ZxOne), dB(ZkOne), dB(ZwOne), Dmax, K2, ssvel #ssvel is not wavelength-dependent

def getDWRs(particleType):
    '''
    simple wrapper around model3fOne
    '''

    Zx,Zk,Zw,Dmax,__,__    = model3fOne(particleType)

    return Zx-Zk,Zk-Zw,Dmax

def getUnambigousDWRdmax(Dmax,DWR,DmaxRetr=5e-3,DWRlowDetect=1,showIllus=False,ax=None):
    '''
    get an Unambigous DWR-Dmax relation
        Dmax: [m] maximum dimensions
        DWR:  [dB] dual-wavelength ratio
        
        DmaxRetr:       #[m] maximum size considered in retrieval; this inexplicitly assumes that larger particles are not relevant
        DWRlowDetect:   #[dB] DWRs smaller than this are disregarded (detection limit)
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
        ax.set_xlim([5e-1,1e1])
        ax.set_ylim([0,np.max(DWR)*1.1])
        ax.set_xlabel("Dmax [mm]")
        ax.set_ylabel("DWR [dB]")
    else:
        ax = None
    
    return DWR,ax
