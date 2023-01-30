import sys
import numpy as np
import os
from pathlib import Path
import math

import concurrent.futures

getelecRootPath = str(Path(__file__).parents[2].absolute())
sys.path.insert(0,getelecRootPath + "/src/")
import getelec as gt

def getArgument(arg, idx):

    if isinstance(arg, (np.ndarray, list)):
        if idx >= len(arg):
            return arg[-1]
        else:
            return arg[idx]
    else:
        return arg

def current_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):
    """ Calculates the current density for an numpy arrays of inputs
    """

    #emitter = gt.ConductionBandEmitter()
    emitter = gt.MetalEmitter()

    currentDensity = np.copy(field)
    kT = gt.Globals.BoltzmannConstant * np.array(temperature)
    
    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(workFunction, i), getArgument(kT, i))
        currentDensity[i] = emitter.currentDensity()
        
    return currentDensity

def heat_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Workfunction [eV] - Material's workfunction
    Temperature [K] - Emitter's temperature
    nh_metal [W/nm^2] - Deposited Nottigham heat

    For more info refer to GETELEC TABULATOR's documentation
    """
    emitter = gt.ConductionBandEmitter()
    
    kT = gt.Globals.BoltzmannConstant * np.array(temperature)
    
    nh_metal = np.copy(field)

    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(workFunction, i), getArgument(kT, i))
    
        nh_metal[i] = emitter.nottinghamHeat()

    return nh_metal

def spectrum_metal_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):
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

def current_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):
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
    j_metal [A/nm^2] - Emitted current density

    For more info refer to GETELEC TABULATOR's documentation
    """

    kT = gt.Globals.BoltzmannConstant * np.array(temperature)

    emitter = gt.SemiconductorEmitter()

    currentDensity = np.copy(field)

    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i), ec, ef, eg, getArgument(kT, i), me, mp)

        currentDensity[i] = emitter.currentDensity()

    return currentDensity

def heat_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):
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
    nh_total [W/nm^2] - Deposited Nottigham heat

    For more info refer to GETELEC TABULATOR's documentation
    """

    kT = gt.Globals.BoltzmannConstant * np.array(temperature)

    emitter = gt.SemiconductorEmitter()

    nh_total = np.copy(field)

    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i), ec, ef, eg, getArgument(kT, i), me, mp)
        
        nh_total[i] = emitter.nottinghamHeat()

    return nh_total

def spectrum_semiconductor_emitter(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, ec: float, ef:float, eg: float, temperature: np.ndarray, me: float, mp: float):
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


#EXPERIMENTAL

def current_density_metal_beta(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.ndarray):
    """ Calculates the current density for an numpy arrays of inputs
        Inputs must be same length
    """

    #emitter = gt.ConductionBandEmitter()
    emitter = gt.MetalEmitter()

    currentDensity = np.copy(field)
    kT = gt.Globals.BoltzmannConstant * np.array(temperature)
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(lambda f, r, g, wf, kt: emitter.barrier.setParameters(f, r, g) or emitter.setParameters(wf, kt) or emitter.currentDensity(), getArgument(field, i),
         getArgument(radius, i), getArgument(gamma, i), getArgument(workFunction, i), getArgument(kT, i)) for i in range(len(field))]
        currentDensity = [f.result() for f in concurrent.futures.as_completed(futures)]
    return currentDensity

def current_density_metal_beta_chunks(field: np.ndarray, radius: np.ndarray, gamma: np.ndarray, workFunction: np.ndarray, temperature: np.array):
    """ Calculates the current density for an numpy arrays of inputs
        Inputs must be same length.
        Uses chunks and multithreading
    """

    def calculate_current_density_chunk(start_index, end_index):

        #emitter = gt.ConductionBandEmitter()
        emitter = gt.MetalEmitter()
        current_density_chunk = np.copy(field[start_index:end_index])
        kT = gt.Globals.BoltzmannConstant * np.array(temperature[start_index:end_index])
        
        for i in range(start_index, end_index):

            emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
            emitter.setParameters(getArgument(workFunction, i), getArgument(kT, i))
            current_density_chunk[i - start_index] = emitter.currentDensity()

        return current_density_chunk

    currentDensity = np.copy(field)
    max_chunk_size = 8
    chunks = range(0, len(field), max_chunk_size)

    with concurrent.futures.ThreadPoolExecutor() as executor:

        results = []

        for i in range(len(chunks)-1):

            start_index = chunks[i]
            end_index = chunks[i+1]
            results.append(executor.submit(calculate_current_density_chunk, start_index, end_index))

        currentDensity = [r.result() for r in concurrent.futures.as_completed(results)]

    return np.concatenate(currentDensity)


#field = np.linspace(3., 6., 32)

#print(current_density_metal(field, 100., 10., 4.5, 300.))