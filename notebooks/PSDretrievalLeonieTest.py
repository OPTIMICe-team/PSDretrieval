

from PSDretrieval import processRadar as pR
from PSDretrieval import plotting as pl
from PSDretrieval import scattering as sc
from PSDretrieval import retrievalUtils as rU
import snowScatt
from snowScatt.instrumentSimulator.radarSpectrum import dopplerSpectrum
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import warnings
warnings.filterwarnings('ignore')
    
heights = [1800, 1700]                                      # in [m]
date = "20181124"
time = "17:11:00"
westimate = [ -0.2, -0.22]                             # in [m/s] downwards
zeThreshold=-35
##################################################
# load the Spectrograms
##################################################

SpecOneTime = pR.loadSpectra(date=date,time=time,tRange=0,
                             dataPath="/data/obs/campaigns/tripex-pol/processed/",
                             loadAllHeights=True,loadSample=False)
###############################################
#- calculate mean in time-height window
###############################################
SpecWindows = xr.Dataset()
mean_SpecWindows = xr.Dataset()
#fig,ax = plt.subplots(nrows=len(heights),ncols=1, figsize=(5,5), sharex=True)
for i,h in enumerate(heights):
    SpecWindow  = pR.loadSpectra(date=date,
                                 time=time,
                                 tRange=0.25,
                                 hRange=80,
                                 hcenter=h,
                                 loadSample=False,
                                 createSample=False,
                                 dataPath="/data/obs/campaigns/tripex-pol/processed/")
    
    SpecWindows = xr.merge([SpecWindows,SpecWindow])#.expand_dims(dim='height').assign_coords(range=[SpecSingle.range]])
    
    # mean spectrum in window
    mean_SpecWindow = SpecWindow.mean(dim=["time","range"])
    mean_SpecWindow = mean_SpecWindow.expand_dims(dim='height').assign_coords(height=[h])
    mean_SpecWindows = xr.merge([mean_SpecWindows,mean_SpecWindow])

################################################
#- adjust by manual estimate for vertical wind
################################################

#fig,ax = plt.subplots(nrows=len(heights),ncols=1, figsize=(5,5), sharex=True)
SpecWindowsMeanShifted = xr.Dataset()
for i,h in enumerate(heights):
    SpecWindowWshifted  = rU.shiftSpecWindow(mean_SpecWindows.sel(height=h), westimate[i])
    #ax[i] = pl.plotObsSpectra(SpecWindowWshifted,ax[i],ls='--')
    #if i != len(heights)-1:
   # 	ax[i].set_xlabel('')
    #ax[i].grid(True)
    #ax[i].set_title('height: {0}'.format(h))
    SpecWindowWshifted = SpecWindowWshifted.expand_dims(dim='height').assign_coords(height=[h])
    SpecWindowsMeanShifted = xr.merge([SpecWindowsMeanShifted,SpecWindowWshifted])
#plt.tight_layout()
#plt.savefig('Spectra_shifted.png')    
#plt.close()
#quit()    

#################################################
#- cut noisy part and plot spectra again
#################################################

SpecWindowsMeanShiftedCut = xr.Dataset()
for h in heights:
    data = pR.cutLowZe(SpecWindowsMeanShifted.sel(height=h),zeThreshold=zeThreshold)
    data = data.expand_dims(dim='height').assign_coords(height=[h])
    
    SpecWindowsMeanShiftedCut = xr.merge([SpecWindowsMeanShiftedCut,data])

'''
#- Compare the snowscatt models of different particle types to the observed DV-DWR
allParticleTypes        = [*snowScatt.snowLibrary._fileList.keys()]
#get a list of all particle-type names from this aggregate type with different riming degrees
allUnrimVT = [k for k in allParticleTypes if ("vonTerzi_" in k and not "rimed" in k)]

fig,ax = plt.subplots(nrows=len(heights),ncols=2, figsize=(7,7), sharex=True, sharey=True)
for i,h in enumerate(heights):
    ax[i,:] = pl.plotSDWRvsDVobs(SpecWindowsMeanShiftedCut.sel(height=h),ax[i,:])
    for pType in allUnrimVT:
        #get spectral-resolved particle properties
        Zx, Zk, Zw, Dmax, K2, vel = sc.model3fOne(pType)
        #calculate spectral DWRs
        DWRxk = Zx-Zk; DWRkw = Zk-Zw
        ax[i,:] = pl.plotSDWRvsDVmodel(vel,DWRxk,DWRkw,ax[i,:],pType)
    if i != len(heights)-1:
        ax[i,0].set_xlabel('')
        ax[i,1].set_xlabel('')
    ax[i,1].set_ylabel('')
    ax[i,0].grid(True)
    ax[i,1].grid(True)
    ax[0,0].legend()
plt.tight_layout()
plt.savefig('all_vonTerzi_unrimed_comp.png')
plt.close()
'''
########################################
# choose best fitting particle type
########################################
#fig,ax = plt.subplots(nrows=len(heights),ncols=2, figsize=(15,15), sharex=False, sharey=False)

whichDWRsToUse="DWR_Ka_W" # choose from ["both","DWR_Ka_W",""DWR_X_Ka]
#find best fitting particle type
ParticleTypeList = [*snowScatt.snowLibrary._fileList.keys()] # read https://www.python.org/dev/peps/pep-0448/ 
ChosenParticles = []
# for the [*...] formalism

for h in heights:
    bestPartType,orderedListPartType = rU.findBestFittingPartType(ParticleTypeList,
                                                                SpecWindowsMeanShiftedCut.sel(height=h),
                                                                whichDWRsToUse=whichDWRsToUse)
    ChosenParticles.append(bestPartType)
    #plot sDWR vs DV for best fitting particle type
    #Zx, Zk, Zw, Dmax, K2, vel = sc.model3fOne(bestPartType)
    #DWRxk = Zx - Zk; DWRkw = Zk - Zw
#    ax[i,:] = pl.plotSDWRvsDVobs(SpecWindows[i],ax[i,:])
#    ax[i,:] = pl.plotSDWRvsDVmodel(vel,DWRxk,DWRkw,ax[i,:],bestPartType)
#plt.show()

##########################
# derive PSD
#########################
#ChosenParticlesVec = [['vonTerzi_plate','vonTerzi_plate'],['vonTerzi_dendrite','vonTerzi_dendrite']]
#ChosenParticles = ['vonTerzi_plate','vonTerzi_plate']
#fig,axes = plt.subplots(nrows=len(heights),ncols=2,figsize=(7,7), sharey=True)


print(ChosenParticles)
PSDsAndStuff = []
for i,h in enumerate(heights):
	Zx, Zk, Zw, Dmax, K2, vel = sc.model3fOne(ChosenParticles[i])
	velObs,NumConNormV,NumConNormD,DmaxAtObsDVgrid = rU.calculateNumberForEachDVbin(
	    Zk,
	    SpecWindowsMeanShiftedCut.sel(height=h).KaSpecH.values,
	    vel,
	    -SpecWindowsMeanShiftedCut.sel(height=h).KaSpecH.doppler.values,
	    DmaxModel=Dmax)
	
	
	PSDsAndStuff.append([velObs,NumConNormV,NumConNormD,DmaxAtObsDVgrid])
	#axes[i,:] = pl.plotNumCon(NumConNormV,NumConNormD,axes[i,:],velObs,DmaxAtObsDVgrid*1e3,label=ChosenParticles[i])
	#axes[i,0].grid()
	#axes[i,1].grid()
	#if i != len(heights)-1:
	#    axes[i,0].set_xlabel('')
	#    axes[i,1].set_xlabel('')
#axes[0,0].legend()	 
#plt.tight_layout()
#plt.savefig('derived_PSD_PVD.png')
#plt.show()
#quit()

##################################################
#- total particle number
##################################################
tpnv_array = []
tpnd_array = []
for i in range(len(heights)):
    # total particle number by doppler velocity
    NumConNormV = PSDsAndStuff[i][1]
    velObs = PSDsAndStuff[i][0]
    NumInVBin = NumConNormV[:-1]*abs(np.diff(velObs))
    tpnv = np.nansum(NumInVBin)
    tpnv_array.append(tpnv)
    # total particle number by diameter
    NumConNormD = PSDsAndStuff[i][2]
    DmaxAtObsDVgrid = PSDsAndStuff[i][3]
    NumInDBin = NumConNormD[:-1]*abs(np.diff(DmaxAtObsDVgrid))
    tpnd = np.nansum(NumInDBin)
    tpnd_array.append(tpnd)

    print(int(tpnv), int(tpnd))
'''
fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(15,8))

fig.suptitle('Total Number of Particles', fontsize=16)

ax[0].plot(tpnv_array,heights)
ax[1].plot(tpnd_array,heights)

ax[0].set_title('total number of particles from PVD')
ax[1].set_title('total number of particles from PSD')
ax[0].set_xlabel('number of particles in [m⁻³]')
ax[0].set_ylabel('height in [m]')
ax[0].grid()
ax[1].set_xlabel('number of particles in [m⁻³]')
ax[1].set_ylabel('height in [m]')
ax[1].grid()

plt.tight_layout()
plt.savefig('number_concentration.png')
plt.show()
'''
#####################################################
#IWC
####################################################
# TO DO: Use particle-type-specific parameters here instead.
a_m = 0.038
b_m = 2.
IWC_array = []
for i in range(len(heights)):
    NumConNormD = PSDsAndStuff[i][2]
    D = PSDsAndStuff[i][3]
    x = np.array([a_m*d**b_m for d in D])
    IWC = np.nansum(NumConNormD[:-1]*x[:-1]*abs(np.diff(D)))    # This means: integral N(D)*m(D) dD.
    IWC_array.append(IWC)
    print(IWC)
'''
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(8,5))
ax.plot(np.multiply(IWC_array,1e3),heights)
ax.set_title('Ice Water Content')
ax.set_xlabel('IWC in [g/m**3]')
ax.set_ylabel('height in [m]')
ax.grid()
plt.tight_layout()
plt.savefig('IWC.png')
plt.show()
'''
######################################################
# crosscheck moments
######################################################
for i,h in enumerate(heights):
    DmaxAtObsDVgrid = PSDsAndStuff[i][3]
    NumConNormD = PSDsAndStuff[i][2]
    bestPartType = ChosenParticles[i]
    SpecSingles_i = SpecWindowsMeanShiftedCut.sel(height=h)
    velObs = PSDsAndStuff[i][0]
    rU.crossCheckIntegratedProp(DmaxAtObsDVgrid,NumConNormD,SpecSingles_i.XSpecH,bestPartType,velObs=velObs)

#####################################################
# crosscheck spectra
#####################################################

fig,ax = plt.subplots(nrows=len(heights),ncols=1, figsize=(5,5), sharex=True)
for i,h in enumerate(heights):
    D = PSDsAndStuff[i][3]
    Nd = PSDsAndStuff[i][2]*np.abs(np.gradient(D)) # need to "renormalise" as it was normalised during calculation of particle concentration. Since we are going to a velocity spectrum, it needs to be multiplied by the gradient to get the total number concentration.
    vel = PSDsAndStuff[i][0]    
    
    ZxOne, ZkOne, ZwOne, Dmax, K2, vel = sc.model3fOne(ChosenParticles[i],Dmax=D,lindB="lin")
    ax[i] = pl.plotObsSpectra(SpecWindowsMeanShiftedCut.sel(height=h),ax[i]) # plot obs
    ax[i].plot(vel,10*np.log10(Nd*ZxOne),color='b',ls='--')
    ax[i].plot(vel,10*np.log10(Nd*ZkOne),color='g',ls='--')
    ax[i].plot(vel,10*np.log10(Nd*ZwOne),color='r',ls='--')
    #ax[i].plot(dopvel[3*i+0],rU.dB(ref_per_height[3*i+0]), color='blue', label='X_retr', linestyle='--')
    ax[i].grid()
    ax[i].legend()
    #print(Nd)
    
plt.tight_layout()
plt.savefig('Spectra_obs_retr_comp.png')
plt.show()


    
    
    
