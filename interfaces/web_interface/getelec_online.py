import sys
import numpy as np
import os
from pathlib import Path
import math

from multiprocessing import Pool

getelecRootPath = str(Path(__file__).parents[2].absolute())
sys.path.insert(0,getelecRootPath + "/src/")

import getelec as gt
from getelecModel import GETELECModel

def current_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):

    model = GETELECModel(emitterType='metal', field=field, radius=radius, gamma=gamma, workFunction=workFunction, temperature=temperature)

    return model.calculateCurrentDensity().currentDensity

def current_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):

    model = GETELECModel(emitterType='semiconductor', field=field, radius=radius, gamma=gamma, ec=ec, ef=ef, eg=eg, temperature=temperature, me=me, mp=mp)

    return model.calculateCurrentDensity().currentDensity

def heat_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):

    model = GETELECModel(emitterType='metal', field=field, radius=radius, gamma=gamma, workFunction=workFunction, temperature=temperature)

    return model.calculateNottinghamHeat().nottinghamHeat

def heat_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):

    model = GETELECModel(emitterType='semiconductor', field=field, radius=radius, gamma=gamma, ec=ec, ef=ef, eg=eg, temperature=temperature, me=me, mp=mp)
    
    return model.calculateNottinghamHeat().nottinghamHeat

def spectrum_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):

    model = GETELECModel(emitterType='metal', field=field, radius=radius, gamma=gamma, workFunction=workFunction, temperature=temperature)

    return model.calculateElectronSpectrum().electronSpectrum

def spectrum_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):

    model = GETELECModel(emitterType='semiconductor', field=field, radius=radius, gamma=gamma, ec=ec, ef=ef, eg=eg, temperature=temperature, me=me, mp=mp)

    return model.calculateElectronSpectrum().electronSpectrum

###CODE BELOW NEEDS UPDATING USES OLD STUFF

def getArgument(arg, idx):

    if isinstance(arg, (np.ndarray, list)):
        if idx >= len(arg):
            return arg[-1]
        else:
            return arg[idx]
    else:
        return arg


def OLDspectrum_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Workfunction [eV] - Material's workfunction
    Temperature [K] - Emitter's temperature
    energy [eV] - energy space (x-axis)
    electron count [#] - number of electrons per energy unit (y-axis)

    For more info refer to GETELEC TABULATOR's documentation
    """

    #emitter = gt.ConductionBandEmitter()
    emitter = gt.MetalEmitter()

    kT = gt.Globals.BoltzmannConstant * np.array(temperature)
    
    energy = []
    electron_count = []

    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(workFunction, i), getArgument(kT, i))
    
        _energy, _electron_count = emitter.totalEnergySpectrumArrays(numberOfPoints=257)

        energy.append(_energy)
        electron_count.append(_electron_count)

    return energy, electron_count

def OLDspectrum_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Ec [eV] - Bottom of the conduction band
    Ef [eV] - Fermi level
    Eg [eV] - band gap
    Temperature [K] - Emitter's temperature
    me [kg] - electron relative mass
    mp [kg] - hole relative mass
    energy_c [eV] - conduction band energy space (x-axis)
    count_c [#] - number of electrons per energy unit from conduction band (y-axis)
    energy_v [eV] - valence band energy space (x-axis)
    count_v [#] - number of electrons per energy unit from valence band (y-axis)

    For more info refer to GETELEC TABULATOR's documentation
    """

    kT = gt.Globals.BoltzmannConstant * np.array(temperature)

    emitter = gt.SemiconductorEmitter()


    #list of arrays

    energy = []
    electron_count = []



    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i), ec, ef, eg, getArgument(kT, i), me, mp)

        _energy, _electron_count = emitter.totalEnergyDistribution()

        energy.append(_energy)
        electron_count.append(_electron_count)

    return energy, electron_count

def fit_data(xML: np.array, yML, workFunction, mode = "simple"):

    if mode == "simple":
        gt.setTabulationPath(getelecRootPath + "/tabulated/1D_1024")
    else:
        gt.setTabulationPath(getelecRootPath + "tabulated/2D_512x256")
    
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

###EXPERIMENTAL

def metal_emitter_worker(input_list):
    """ Worker function for calculating current density of a single element from the field array
    """

    field, radius, gamma, workFunction, temperature = input_list

    emitter = gt.MetalEmitter()
    kT = gt.Globals.BoltzmannConstant * np.array(temperature)

    emitter.barrier.setParameters(field, radius, gamma)
    emitter.setParameters(workFunction, kT)
    currentDensity = emitter.currentDensity()

    return currentDensity

def current_metal_emitter_threaded(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray, threads = 8):
    """ Calculates the current density for an numpy arrays of inputs
        Uses multithreading via Pool
    """
    with Pool(processes=threads) as process_pool:
        input_list = [(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i), getArgument(workFunction, i), getArgument(temperature, i)) for i in range(len(field))]
        results = process_pool.map(metal_emitter_worker, input_list)
    return results