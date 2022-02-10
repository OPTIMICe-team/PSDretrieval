from PSDretrieval import processRadar as pR
from PSDretrieval import plotting as pl
from PSDretrieval import scattering as sc
from PSDretrieval import retrievalUtils as rU
import matplotlib.pyplot as plt
import numpy as np
from snowScatt import snowMassVelocityArea
from IPython.terminal.debugger import set_trace

#define time
#date = "20190122"
#time = "14:55"

date    = "20190113"
time    = "06:18:04"
hRange  = 1600


#load spectra
#SpecWindow  = pR.loadSpectra()
SpecWindow  = pR.loadSpectra(loadSample=False,dataPath="/data/obs/campaigns/tripex-pol/processed/",createSample=True,date=date,time=time,tRange=1,hRange=180,hcenter=hRange)
SpecSingle  = pR.selectSingleTimeHeight(SpecWindow)
DWRkey="DWR_X_Ka"

#plot Spectra and sDWR
fig,ax = plt.subplots(nrows=1,ncols=1)
fig2,axes = plt.subplots(nrows=1,ncols=3)
__ = pl.plotObsSpectra(SpecSingle,ax)
__ = pl.plotSpectralDWR(SpecSingle[DWRkey],ax)
plt.show()

#get scattering properties
particleType = "vonTerzi_mixcoldend"
DWRxk,DWRkw,Dmax = sc.getDWRs(particleType)
if DWRkey=="DWR_Ka_W":
    DWR = DWRkw
    DmaxRetr = 5e-3 #[m] maximum size considered in retrieval; this inexplicitly assumes that larger particles are not relevant
elif DWRkey=="DWR_X_Ka":
    DWR = DWRxk
    DmaxRetr = 1e-2 #[m] maximum size considered in retrieval; this inexplicitly assumes that larger particles are not relevant


DWRUnamb,__ = sc.getUnambigousDWRdmax(Dmax,DWR,DmaxRetr=DmaxRetr,DWRlowDetect=1.,showIllus=False)

#get Dmax from sDWR (spectral resolved)
DmaxfromDWR = rU.getDmaxFromSDWR(SpecSingle[DWRkey],DWRUnamb,Dmax,showIllus=False)

#get the single particle reflectivity
ZkOne = sc.getSinglePartRefl(particleType,DmaxfromDWR,freq="k")

#get the number concentration from the spectrum and the single particle scattering properties at the given Doppler velocity bin
Nnorm,__ = rU.calcNumberConcFromSpectrumAndZOne(SpecSingle.KaSpecH,ZkOne,showIllus=False)

__ = rU.histDWRandDmaxVsDv(SpecWindow[DWRkey],SpecWindow.KaSpecH,SpecWindow.KaSpecNoiseH,DWRUnamb,Dmax,aboveNoiseThreshold=15,showIllus=True,ax=axes,fig=fig)
