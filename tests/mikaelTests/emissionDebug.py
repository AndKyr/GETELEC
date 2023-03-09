import os
import sys
import numpy as np

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/src/'
sys.path.append(pythonpath)

from getelec import GETELECModel as GTM

gtm = GTM()

field = np.array([1, 2, 3, 4, 5])
radius = np.full(5, 50)
gamma = np.full(5, 10)
workFunction = np.full(5, 4)
temperature = np.full(5, 300)

conductionBandBottom = 4.05
fermiLevel = 4.61
bandGap = 1.12
effectiveMassConduction = 0.7
effectiveMassValence = 0.5

numberOfSpectrumPoints = 8

gtm = GTM(emitterType = 'metal', field = field, radius = radius, gamma = gamma, workFunction = workFunction, temperature = temperature, numberOfSpectrumPoints = numberOfSpectrumPoints)


gtm.run(calculateCurrent=True, calculateNottinghamHeat=True, calculateSpectrum=True)

print(gtm.currentDensity)