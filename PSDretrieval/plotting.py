#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
all plotting routines
'''
from IPython.terminal.debugger import set_trace
from PSDretrieval import scattering as sc

def plotObsSpectra(xrSpec,ax):
    '''
    plot observed spectra
    Arguments:
        INPUT:
            xrDWR: xarray containing all spectra information (including DWRs with labels "DWR_X_Ka" and "DWR_Ka_W")
        IN- & OUTPUT: 
            ax: axis handle
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
    Arguments:
        INPUT:
            xrDWR: xarray array containing only the DWR (no matter which frequency)
        IN- & OUTPUT: 
            ax: axis handle
    '''

    ax.plot(-xrDWR.doppler,xrDWR)
    ax.set_xlim([-0.5,2])
    ax.set_xlabel("DV [m/s]")
    ax.set_ylabel("DWR$_{Ka,W}$ [dB]")

    return ax 

def plotSDWRvsDVobs(xrSpec,axes):
    '''
    plot Doppler velocity vs. spectral DWR of X-Ka and Ka-W band combinations (kind of a v-D plot)
    Arguments:
        INPUT:
            xrSpec [can be single spectra or a time-height window]: xarray containing all spectra information (including DWRs with labels "DWR_X_Ka" and "DWR_Ka_W")
                if single spectra: just plot DV vs. DWRs
                if Window:         average over (vertical wind corrected) DV and plot vs DWR
        IN- & OUTPUT: 
            axes: axes handles
    '''

    #these two DWR variables are needed here
    DWRkeys = ["DWR_X_Ka","DWR_Ka_W"]

    xrSpecPlot = xrSpec.copy() #this seems necessary for the way the averaging is applied, but there might be a cleaner way

    if xrSpecPlot["KaSpecH"].values.ndim>1:
        print("plot average DV vs DWR for a time-height window")
        for key in DWRkeys:
            xrSpecPlot[key] = xrSpecPlot[key].mean(dim=["time","range"])
    else:
        print("plot DV vs DWR for a single spectrum")
    
    for i_ax,(ax,key) in enumerate(zip(axes,DWRkeys)):
        axes[i_ax].plot(xrSpecPlot[key],-xrSpecPlot.doppler,c="k",lw=5,ls="--",label="obs.")
        ax.set_ylabel("DV [m/s]")
    axes[0].legend()
    axes[0].set_xlabel("DWR$_{X,Ka}$ [dB]")
    axes[1].set_xlabel("DWR$_{Ka,W}$ [dB]")

    return axes

def plotSDWRvsDVmodel(vel,DWRxk,DWRkw,axes,pType):
    '''
    plot Doppler velocity vs. spectral DWR of X-Ka and Ka-W band combinations (kind of a v-D plot) from the snowScatt simulations
    Arguments:
        INPUT: (all spectral variables)
            vel     [numpy-array]:  velocity [m/s]
            DWRxk   [numpy-array]:  dual-wavelength ratio (X-Ka band) [dB]
            DWRkw   [numpy-array]:  dual-wavelength ratio (Ka-W band) [dB]
            pType   [str]:          particle-type name
        IN- & OUTPUT: 
            axes: axes handles
    '''

    for i_ax,(ax,DWR) in enumerate(zip(axes,[DWRxk,DWRkw])):
        axes[i_ax].plot(DWR,vel,label=pType)
        ax.set_ylabel("DV [m/s]")

    axes[0].set_xlabel("DWR$_{X,Ka}$ [dB]")
    axes[1].set_xlabel("DWR$_{Ka,W}$ [dB]")

    return axes

def plotSinglePartZe(pType,ax,freq="Ka"):
    '''
    plot single particle backscattering coefficient (needed for retrieval) vs. velocity
    Arguments:
        INPUT:
            pType[str]:     particle type from snowScatt
            (optional)
            freq [str]:     frequency (options are ["X","Ka","W"])
        IN- & OUTPUT: 
            ax: axes handles
    '''

    #get properties from snowScatt
    Zx, Zk, Zw, Dmax, K2, vel = sc.model3fOne(pType)

    if freq=="X":
        Zone = Zx
    elif freq=="Ka":
        Zone = Zk
    elif freq=="W":
        Zone = Zw

    ax.plot(vel,Zone)
    ax.set_xlabel("DV [m/s]")
    ax.set_ylabel("Zone [dB]")

    return ax

    
def plotNumCon(NumConNormV,NumConNormD,axes,vVec,DmaxVec):
    '''
    plot number concentration against some other values
    Arguments:
        INPUT: (all spectral variables)
            NumConNormV [numpy-array, #/m^3]:  number concentration normed by velocity
            NumConNormD [numpy-array, #/m^3]:  number concentration normed by diameter
            vVec [numpy-array, m/s]:           velocity array corresponding to above NumConNormV
            Dmax [numpy-array, m]:             diameter array corresponding to above NumConNormD
        IN - & OUTPUT: 
            axes: axes handles
    '''

    #loop over all x-variables
    for xVar,xLabel,ax,NumConNorm,ylabel in zip([vVec,DmaxVec],["vel [m/s]","Dmax [mm]"],axes,[NumConNormV,NumConNormD],["N [#/m$^3 / (m/s)$]","N [#/m$^3 / (m)$]"]):
        #plot
        ax.semilogy(xVar,NumConNorm)
        #set labels
        ax.set_xlabel(xLabel)
        ax.set_ylabel(ylabel)

    return axes
def plotSpectraObsAllHeights(data,xlim=[-5,1],ylim=[0,8000]):
    import matplotlib.pyplot as plt
    '''
    plot spectrogramms (for all heights)
    Arguments:
        INPUT: 
            xrSpec: dataarray which contains height,doppler series of spectrum (already one time step has to be selected)
            xlim: Doppler velocity limits
            ylim= range limits
        OUTPUT: 
            ax: returned are the axes with the plot
            
    '''        
    #- plot the data
    fig,axes = plt.subplots(figsize=(15,5),ncols=4,sharey=True)
    radData = {'XSpecH':{'data':data['XSpecH'],'title':'X','axis':axes[0],'lim':(-40,10)},
               'KaSpecH':{'data':data['KaSpecH'],'title':'Ka','axis':axes[1],'lim':(-40,10)},
               'WSpecH':{'data':data['WSpecH'],'title':'W','axis':axes[2],'lim':(-40,10)},
               'KaW':{'data':data['KaSpecH']-data['WSpecH'],'title':'DWR-KaW','axis':axes[3],'lim':(-5,10)}}
    for rad in radData.keys():
        
        plot = radData[rad]['axis'].pcolormesh(data.doppler.values,
                              data.range.values,
                              radData[rad]['data'].values[0],
                              vmin=radData[rad]['lim'][0],vmax=radData[rad]['lim'][1],
                              cmap='nipy_spectral')
        if rad == 'KaW':
            fig.colorbar(plot,ax=radData[rad]['axis'],label='[dB]')
        if rad == 'XSpecH':
            radData[rad]['axis'].set_ylabel(r'range [m]')
        
        radData[rad]['axis'].grid()
        radData[rad]['axis'].set_title(radData[rad]['title'])
        radData[rad]['axis'].set_xlabel(r'Doppler velocity [ms$^{-1}]$')
        radData[rad]['axis'].set_xlim(xlim)
        radData[rad]['axis'].set_ylim(ylim)
    plt.tight_layout()
    #plt.show()
    #return fig
