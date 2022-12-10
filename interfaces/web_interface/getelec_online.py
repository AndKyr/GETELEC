import sys
import numpy as np
import os
from pathlib import Path

getelecRootPath = str(Path(__file__).parents[2].absolute())
sys.path.insert(0,getelecRootPath + "/src/")
import getelec as gt

def getArgument(array, index):
    """
    Array - input data array
    index - position of value in array

    returns value at array[index] if it exists, else returns whole array
    
    """
    try:
        return(array[index])
    except(TypeError, IndexError):
        return array

def currentDensityMetal(field, radius, gamma, workFunction, temperature):
    """ Calculates the current density
        Inputs must be same length
    """

    emitter = gt.MetalEmitter()

    currentDensity = np.copy(field)

    kT = gt.Globals.BoltzmannConstant * temperature
    
    for i in range(len(field)):

        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))

        emitter.setParameters(getArgument(workFunction, i), getArgument(kT, i))

        currentDensity[i] = emitter.currentDensityFast()
        
    return currentDensity

def nottinghamHeatMetal(field, radius, gamma, workfunction, temperature):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Workfunction [eV] - Material's workfunction
    Temperature [K] - Emitter's temperature
    nh_metal [W/nm^2] - Deposited Nottigham heat

    For more info refer to GETELEC TABULATOR's documentation
    """

    metal_emitter = gt.MetalEmitter()
   
    kT = gt.Globals.BoltzmannConstant * temperature
    
    nh_metal = np.copy(field)

    for i in range(len(field)):

        metal_emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
    
        metal_emitter.setParameters(getArgument(workfunction, i), getArgument(kT, i))
    
        nh_metal[i] = metal_emitter.nottinghamHeatFast()

    return nh_metal

def electronSpectrumMetal(field, radius, gamma, workfunction, temperature):
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
    metal_emitter = gt.MetalEmitter()

    kT = gt.Globals.BoltzmannConstant * temperature
    
    energy = np.copy(field)
    electron_count = np.copy(field)

    for i in range(len(field)):

        metal_emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
    
        metal_emitter.setParameters(getArgument(workfunction, i), getArgument(kT, i))
    
        energy[i], electron_count[i] = metal_emitter.normalEnergyDistribution()

    return energy, electron_count

def current_semiconductor_emitter(field, radius, gamma, ec, ef, eg, temperature, me, mp):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Ec [eV] - Bottom of the conduction band
    Ef [eV] - Fermi level
    Eg [eV] - band gap
    Temperature [K] - Emitter's temperature
    me [kg] - electron effective mass
    mp [kg] - hole effective mass
    j_metal [A/nm^2] - Emitted current density

    For more info refer to GETELEC TABULATOR's documentation
    """
    semiconductor_emitter = gt.Semiconductor_Emitter()

    kT = gt.Globals.BoltzmannConstant * temperature

    j_total = np.copy(field)
    j_c = np.copy(field)
    j_v = np.copy(field)
    m = np.ones(field) * gt.Globals.ElectronMass

    for i in range(len(field)):

        semiconductor_emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(getArgument(ec, i), getArgument(ef, i), getArgument(eg, i), getArgument(kT, i), getArgument(m, i), getArgument(me, i), getArgument(mp, i))
        
        j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()

    return j_total

def heat_semiconductor_emitter(field, radius, gamma, ec, ef, eg, temperature, me, mp):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Ec [eV] - Bottom of the conduction band
    Ef [eV] - Fermi level
    Eg [eV] - band gap
    Temperature [K] - Emitter's temperature
    me [kg] - electron effective mass
    mp [kg] - hole effective mass
    nh_total [W/nm^2] - Deposited Nottigham heat

    For more info refer to GETELEC TABULATOR's documentation
    """

    semiconductor_emitter = gt.Semiconductor_Emitter()

    kT = gt.Globals.BoltzmannConstant * temperature

    nh_total = np.copy(field)
    nh_c = np.copy(field)
    nh_v = np.copy(field)
    m = np.ones(field) * gt.Globals.ElectronMass

    for i in range(len(field)):

        semiconductor_emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(getArgument(ec, i), getArgument(ef, i), getArgument(eg, i), getArgument(kT, i), getArgument(m, i), getArgument(me, i), getArgument(mp, i))
        
        nh_c[i], nh_v[i], nh_total[i] = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

    return nh_total

def spectrum_semiconductor_emitter(field, radius, gamma, ec, ef, eg, temperature, me, mp):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Ec [eV] - Bottom of the conduction band
    Ef [eV] - Fermi level
    Eg [eV] - band gap
    Temperature [K] - Emitter's temperature
    me [kg] - electron effective mass
    mp [kg] - hole effective mass
    energy_c [eV] - conduction band energy space (x-axis)
    count_c [#] - number of electrons per energy unit from conduction band (y-axis)
    energy_v [eV] - valence band energy space (x-axis)
    count_v [#] - number of electrons per energy unit from valence band (y-axis)

    For more info refer to GETELEC TABULATOR's documentation
    """
    semiconductor_emitter = gt.Semiconductor_Emitter()

    kT = gt.Globals.BoltzmannConstant * temperature

    energy_c = np.copy(field)
    count_c = np.copy(field)
    energy_v = np.copy(field)
    count_v = np.copy(field)
    m = np.ones(field) * gt.Globals.ElectronMass

    for i in range(len(field)):

        semiconductor_emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(getArgument(ec, i), getArgument(ef, i), getArgument(eg, i), getArgument(kT, i), getArgument(m, i), getArgument(me, i), getArgument(mp, i))
        
        energy_c[i], count_c[i], energy_v[i], count_v[i] = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

    return energy_c, count_c, energy_v, count_v

def fit_data(xML, yML, workFunction, mode = "simple"):

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