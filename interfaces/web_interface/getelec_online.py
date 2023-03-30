import sys
import numpy as np
import os
from pathlib import Path
import math

from multiprocessing import Pool

getelecRootPath = str(Path(__file__).parents[2].absolute())
sys.path.insert(0,getelecRootPath + "/src/")

import getelec as gt

def current_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):

    model = gt.GETELECModel(emitterType='metal', field=field, radius=radius, gamma=gamma, workFunction=workFunction, temperature=temperature)

    return model.calculateCurrentDensity().currentDensity

def current_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):

    model = gt.GETELECModel(emitterType='semiconductor', field=field, radius=radius, gamma=gamma, ec=ec, ef=ef, eg=eg, temperature=temperature, me=me, mp=mp)

    return model.calculateCurrentDensity().currentDensity

def heat_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):

    model = gt.GETELECModel(emitterType='metal', field=field, radius=radius, gamma=gamma, workFunction=workFunction, temperature=temperature)

    return model.calculateNottinghamHeat().nottinghamHeat

def heat_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):

    model = gt.GETELECModel(emitterType='semiconductor', field=field, radius=radius, gamma=gamma, ec=ec, ef=ef, eg=eg, temperature=temperature, me=me, mp=mp)
    
    return model.calculateNottinghamHeat().nottinghamHeat

def spectrum_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):

    model = gt.GETELECModel(emitterType='metal', field=field, radius=radius, gamma=gamma, workFunction=workFunction, temperature=temperature)

    return model.calculateElectronSpectrum().electronSpectrum

def spectrum_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):

    model = gt.GETELECModel(emitterType='semiconductor', field=field, radius=radius, gamma=gamma, ec=ec, ef=ef, eg=eg, temperature=temperature, me=me, mp=mp)

    return model.calculateElectronSpectrum().electronSpectrum

###CODE BELOW NEEDS UPDATING USES OLD STUFF

def fit_data(xML: np.array, yML, workFunction, mode = "simple"):

    if mode == "simple":
        gt._setTabulationPath(getelecRootPath + "/tabulated/1D_1024")
    else:
        gt._setTabulationPath(getelecRootPath + "tabulated/2D_512x256")
    
    voltage = 1./xML
    current = np.exp(yML)

    fitter = gt.IVDataFitter()

    fitter.setIVcurve(voltageData=voltage, currentData=current)
    fitter.setParameterRange()
    fitter.fitIVCurve()


    xth = np.linspace(min(xML),max(xML),32)
    yth = fitter.getOptCurrentCurve(1./xth)


    xplot = xML / fitter.parameters["fieldConversionFactor"]
    xplot_th = xth / fitter.parameters["fieldConversionFactor"]
        
    
    return xplot, xplot_th, yth, fitter.parameters["fieldConversionFactor"], fitter.parameters["radius"], fitter.prefactor
