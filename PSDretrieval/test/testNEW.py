from PSDretrieval import processRadar as pR
from PSDretrieval import plotting as pl
from PSDretrieval import scattering as sc
from PSDretrieval import retrievalUtils as rU
import matplotlib.pyplot as plt
import numpy as np
from snowScatt import snowMassVelocityArea
import snowScatt
from IPython.terminal.debugger import set_trace


#define time
date    = "20190113"
time    = "06:18:04"
hRange  = 1600

###load Data
#load spectra
SpecWindow  = pR.loadSpectra()
#SpecWindow  = pR.loadSpectra(loadSample=False,dataPath="/data/obs/campaigns/tripex-pol/processed/",createSample=True,date=date,time=time,tRange=1,hRange=180,hcenter=hRange)

PeaksWindow  = pR.loadPeaks()
#PeaksWindow  = pR.loadPeaks(loadSample=False,dataPath="/data/obs/campaigns/tripex-pol/spectralPeaks/",createSample=True,date=date,time=time,tRange=1,hRange=180,hcenter=hRange)

#get vertical wind information from the Spectral data
SpecWindow = pR.addVerticalWindToSpecWindow(SpecWindow,PeaksWindow)
SpecSingle  = pR.selectSingleTimeHeight(SpecWindow)
SpecSingleWshifted  = pR.shiftSpectra(SpecSingle)

#plot Spectra and sDWR
#fig,ax = plt.subplots(nrows=1,ncols=1)
fig2,axes2 = plt.subplots(nrows=1,ncols=2)
#__ = pl.plotObsSpectra(SpecSingleWshifted,ax)
#__ = pl.plotSDWRvsDVobs(SpecSingle,axes2)

#apply stricter noise threshold (needed for particleType selection of snowScatt models)
SpecWindow = pR.cutLowZe(SpecWindow,zeThreshold=-20)
__ = pl.plotSDWRvsDVobs(SpecWindow,axes2)

#get names of all particle types
allParticleTypes        = [*snowScatt.snowLibrary._fileList.keys()] #read https://www.python.org/dev/peps/pep-0448/ for the [*...] formalism
allRimDegr = [k for k in allParticleTypes if 'vonTerzi_mixcoldend' in k]
selectedParticleTypes   = ["vonTerzi_mixcoldend","vonTerzi_mixcoldend_rimed05"]

#find best fitting particle type
[bestPartType,orderedListPartType] = rU.findBestFittingPartType(allParticleTypes,SpecWindow)

#plot sDWR vs DV for best fitting particle type
Zx, Zk, Zw, Dmax, K2, vel = sc.model3fOne(bestPartType)
DWRxk = Zx - Zk; DWRkw = Zk - Zw
axes2 = pl.plotSDWRvsDVmodel(vel,DWRxk,DWRkw,axes2,bestPartType)

##for pType in allParticleTypes:
##for pType in allRimDegr:
#for pType in selectedParticleTypes:
#    #get particle properties
#    Zx, Zk, Zw, Dmax, K2, vel = sc.model3fOne(pType)
#    #calculate spectral DWRs
#    DWRxk = Zx-Zk; DWRkw = Zk-Zw
#    #plot sDWR vs DV
#    axes2 = pl.plotSDWRvsDVmodel(vel,DWRxk,DWRkw,axes2,pType)
    

axes2[0].legend()

plt.show()
