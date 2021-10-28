from PSDretrieval import processRadar as pR
from PSDretrieval import plotting as pl
from PSDretrieval import scattering as sc
import matplotlib.pyplot as plt

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
DWRkwUnamb = sc.getUnambigousDWRdmax(Dmax,DWRkw,DmaxRetr=5e-3,DWRlowDetect=1.,showIllus=False,ax=ax)


