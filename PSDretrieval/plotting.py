#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
all plotting routines
'''

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
