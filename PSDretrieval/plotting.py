#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
all plotting routines
'''
from IPython.terminal.debugger import set_trace

def plotObsSpectra(xrSpec,ax):
    '''
    plot observed spectra
    '''

    #plot simple spectrum
    colors      = ["b","g","r"]
    labels      = ["X","Ka","W"]
    
    for i_var,var in enumerate(["XSpecH","KaSpecH","WSpecH"]):
        ax.plot(-xrSpec[var]["doppler"],xrSpec[var],color=colors[i_var],label=labels[i_var],lw=1)
        ax.legend()
        ax.set_xlabel("v [m/s]")
        ax.set_ylabel("z [dBz]")
        ax.set_xlim([-0.5,3.0])

    return ax

def plotSpectralDWR(xrDWR,ax):
    '''
    plot the spectral DWR
    xrDWR: xarray array containing only the DWR (no matter which frequency)
    '''

    ax.plot(-xrDWR.doppler,xrDWR)
    ax.set_xlim([-0.5,2])
    ax.set_xlabel("DV [m/s]")
    ax.set_ylabel("DWR$_{Ka,W}$ [dB]")

    return ax 
