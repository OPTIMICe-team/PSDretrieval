from PSDretrieval import processRadar as pR
from PSDretrieval import plotting as pl
from PSDretrieval import scattering as sc
from PSDretrieval import retrievalUtils as rU
import matplotlib.pyplot as plt
import numpy as np
from snowScatt import snowMassVelocityArea

#load spectra
SpecWindow  = pR.loadSpectra()
SpecSingle  = pR.selectSingleTimeHeight(SpecWindow)

#plot Spectra and sDWR
fig,ax = plt.subplots(nrows=1,ncols=1)
ax = pl.plotObsSpectra(SpecSingle,ax)
ax = pl.plotSpectralDWR(SpecSingle.DWR_Ka_W,ax)

#get scattering properties
particleType = "vonTerzi_mixcoldend"
DWRxk,DWRkw,Dmax = sc.getDWRs(particleType)
DWRkwUnamb,__ = sc.getUnambigousDWRdmax(Dmax,DWRkw,DmaxRetr=5e-3,DWRlowDetect=1.,showIllus=False,ax=ax)

#get Dmax from sDWR (spectral resolved)
DmaxfromDWR = rU.getDmaxFromSDWR(SpecSingle.DWR_Ka_W,DWRkwUnamb,Dmax,showIllus=True,ax=None)

#get the single particle reflectivity
ZkOne = sc.getSinglePartRefl(particleType,DmaxfromDWR,freq="k")

#get the number concentration from the spectrum and the single particle scattering properties at the given Doppler velocity bin
Nnorm,ax = rU.calcNumberConcFromSpectrumAndZOne(SpecSingle.KaSpecH,ZkOne,showIllus=False,ax=None)

ax = rU.histDWRandDmaxVsDv(SpecWindow.DWR_Ka_W,SpecWindow.KaSpecH,SpecWindow.KaSpecHspecNoise,aboveNoiseThreshold=30)
