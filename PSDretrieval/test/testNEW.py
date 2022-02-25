from PSDretrieval import processRadar as pR
from PSDretrieval import plotting as pl
from PSDretrieval import scattering as sc
from PSDretrieval import retrievalUtils as rU
import matplotlib.pyplot as plt
import numpy as np
from snowScatt import snowMassVelocityArea
import snowScatt
from IPython.terminal.debugger import set_trace


loadDefault = True
###define time
date    = "20190113"
if loadDefault:
    #unrimed 
    time    = "06:18:04"
    hRange  = 1600
else:
    #rimed
    time    = "06:57:00"
    hRange  = 2000


###load Data
#load spectra
if loadDefault:
    SpecWindow  = pR.loadSpectra()
    PeaksWindow  = pR.loadPeaks()
else:
    SpecWindow  = pR.loadSpectra(loadSample=False,dataPath="/data/obs/campaigns/tripex-pol/processed/",createSample=True,date=date,time=time,tRange=1,hRange=180,hcenter=hRange)

    PeaksWindow  = pR.loadPeaks(loadSample=False,dataPath="/data/obs/campaigns/tripex-pol/spectralPeaks/",createSample=True,date=date,time=time,tRange=1,hRange=180,hcenter=hRange)

#get vertical wind information from the Spectral data
SpecWindow = pR.addVerticalWindToSpecWindow(SpecWindow,PeaksWindow)
SpecSingle  = pR.selectSingleTimeHeight(SpecWindow)
SpecSingleWshifted  = pR.shiftSpectra(SpecSingle)

#plot Spectra and sDWR
fig,ax = plt.subplots(nrows=1,ncols=1)
#fig2,axes2 = plt.subplots(nrows=1,ncols=2)
__ = pl.plotObsSpectra(SpecSingleWshifted,ax)
#__ = pl.plotSDWRvsDVobs(SpecSingle,axes2)

#apply stricter noise threshold (needed for particleType selection of snowScatt models)
SpecWindow = pR.cutLowZe(SpecWindow,zeThreshold=-20)
#__ = pl.plotSDWRvsDVobs(SpecWindow,axes2)

#get names of all particle types
AllParticleTypes   = [*snowScatt.snowLibrary._fileList.keys()] #read https://www.python.org/dev/peps/pep-0448/ for the [*...] formalism
ParticleTypesList   = AllParticleTypes
#ParticleTypesList   = [k for k in allParticleTypes if 'vonTerzi_mixcoldend' in k]
#ParticleTypesList   = ["vonTerzi_mixcoldend","vonTerzi_column"]
#ParticleTypesList   = ["vonTerzi_mixcoldend_rimed02","vonTerzi_mixcoldend_rimed04"]

#find best fitting particle type
#[bestPartType,orderedListPartType] = rU.findBestFittingPartType(ParticleTypesList,SpecWindow,verbose=True,whichDWRsToUse="DWR_Ka_W")
bestPartType = "vonTerzi_column"

#plot sDWR vs DV for best fitting particle type
ZxModel, ZkModel, ZwModel, Dmax, K2, velModel = sc.model3fOne(bestPartType)
DWRxk = ZxModel - ZkModel; DWRkw = ZkModel - ZwModel
DWRxkModelUnAmb,ax = sc.getUnambigousDWRdmax(Dmax,DWRxk)
DWRkwModelUnAmb,ax = sc.getUnambigousDWRdmax(Dmax,DWRkw)
#axes2 = pl.plotSDWRvsDVmodel(velModel,DWRxkModelUnAmb,DWRkwModelUnAmb,axes2,bestPartType)
#axes2[0].legend()

#display the single particle reflectivity
fig3,ax = plt.subplots(nrows=1,ncols=1)
ax = pl.plotSinglePartZe(bestPartType,ax,freq="Ka")

##no, finally calculate the size distribution
velObs,NumConNorm,DmaxAtObsDVgrid = rU.calculateNumberForEachDVbin(ZkModel,SpecSingleWshifted.KaSpecH.values,velModel,-SpecSingleWshifted.KaSpecH.doppler.values,DmaxModel=Dmax)
#plot the number concentration vs. velocity and Dmax
fig4,axes3 = plt.subplots(nrows=1,ncols=2)
axes3 = pl.plotNumCon(NumConNorm,axes3,[velObs,DmaxAtObsDVgrid*1e3],["vel [m/s]","Dmax [mm]"])

plt.show()
