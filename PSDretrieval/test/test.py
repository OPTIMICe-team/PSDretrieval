from PSDretrieval import processRadar as pR
import matplotlib.pyplot as plt

#load spectra
SpecWindow  = pR.loadSpectra()
SpecSingle  = pR.selectSingleTimeHeight(SpecWindow)

#plot Spectra
fig,ax = plt.subplots(nrows=1,ncols=1)
ax = pR.plotObsSpectra(SpecSingle,ax)

