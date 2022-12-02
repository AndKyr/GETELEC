import sys
import numpy as np
import os
from pathlib import Path

getelecRootPath = str(Path(__file__).parents[2].absolute())
sys.path.insert(0,getelecRootPath + "/src/")
import getelec as gt

def getArgument(arg, index):
    try:
        return(arg[index])
    except(TypeError, IndexError) as error:
        return arg

def currentDensityMetal(field, radius, gamma, workFunction, temperature):
    """ Calculates the current density for an array of inputs 
    """

    emitter = gt.MetalEmitter(gt.Barrier(), gt.Supply())

    currentDensity = np.copy(field)
    kT = gt.BoltzmannConstant * temperature
    
    for i in range(len(field)):
        emitter.barrier.setParameters(getArgument(field, i), getArgument(radius, i), getArgument(gamma, i))
        emitter.setParameters(getArgument(workFunction, i), getArgument(kT, i))
        currentDensity[i] = emitter.currentDensityFast()
        
    return currentDensity

def heat_metal_emitter(Field, Radius, Gamma, Workfunction, Temperature):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Workfunction [eV] - Material's workfunction
    Temperature [K] - Emitter's temperature
    nh_metal [W/nm^2] - Deposited Nottigham heat

    For more info refer to GETELEC TABULATOR's documentation
    """
    tab = gt.Interpolator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = gt.Metal_Emitter(tab)
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter._emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter._emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
    
        nh_metal[i] = metal_emitter.Nottingham_Heat_from_Metals()

    return nh_metal

def spectrum_metal_emitter(Field, Radius, Gamma, Workfunction, Temperature):
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
    tab = gt.Interpolator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = gt.Metal_Emitter(tab)
    
    energy = np.copy(Field)
    electron_count = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter._emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter._emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
    
        energy[i], electron_count[i] = metal_emitter.Energy_Distribution_for_Metals()

    return energy, electron_count

def current_semiconductor_emitter(Field, Radius, Gamma, Ec, Ef, Eg, Temperature, me, mp):
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
    tab = gt.Interpolator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    semiconductor_emitter = gt.Semiconductor_Emitter(tab)

    j_total = np.copy(Field)
    j_c = np.copy(Field)
    j_v = np.copy(Field)
    m = np.ones(Field)*9.1093837015e-31 

    for i in range(len(Field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
        j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()

    return j_total

def heat_semiconductor_emitter(Field, Radius, Gamma, Ec, Ef, Eg, Temperature, me, mp):
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
    tab = gt.Interpolator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    semiconductor_emitter = gt.Semiconductor_Emitter(tab)

    nh_total = np.copy(Field)
    nh_c = np.copy(Field)
    nh_v = np.copy(Field)
    m = np.ones(Field)*9.1093837015e-31 

    for i in range(len(Field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
        nh_c[i], nh_v[i], nh_total[i] = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

    return nh_total

def spectrum_semiconductor_emitter(Field, Radius, Gamma, Ec, Ef, Eg, Temperature, me, mp):
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
    tab = gt.Interpolator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    semiconductor_emitter = gt.Semiconductor_Emitter(tab)

    energy_c = np.copy(Field)
    count_c = np.copy(Field)
    energy_v = np.copy(Field)
    count_v = np.copy(Field)
    m = np.ones(Field)*9.1093837015e-31 

    for i in range(len(Field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
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



field = np.linspace(3., 6., 32)

print(currentDensityMetal(field, 100., 10., 4.5, 300.))