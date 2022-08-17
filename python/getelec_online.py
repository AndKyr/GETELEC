import getelec_tabulator as gt_tab
import getelec_mod as gt_mod
import numpy as np

#emitted current
def current_metal_emitter(Field, Radius, Gamma, Workfunction, Temperature):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Workfunction [eV] - Material's workfunction
    Temperature [K] - Emitter's temperature
    j_metal [A/nm^2] - Emitted current density

    For more info refer to GETELEC TABULATOR's documentation
    """
    tab = gt_tab.Tabulator()

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = gt_tab.Metal_Emitter(tab)

    j_metal = np.copy(Field)
    
    for i in range(len(Field)):

        metal_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Emitter_Parameters(Workfunction[i], kT[i])
    
        j_metal[i] = metal_emitter.Current_Density()
        
    return j_metal

#nottingham heat
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
    tab = gt_tab.Tabulator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = gt_tab.Metal_Emitter(tab)
    
    nh_metal = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Emitter_Parameters(Workfunction[i], kT[i])
    
        nh_metal[i] = metal_emitter.Nottingham_Heat()

    return nh_metal

#electric spectrum - WIP
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
    tab = gt_tab.Tabulator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature
    
    metal_emitter = gt_tab.Metal_Emitter(tab)
    
    energy = np.copy(Field)
    electron_count = np.copy(Field)

    for i in range(len(Field)):
        metal_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter.emitter.Interpolate_Gammow()
    
        metal_emitter.Define_Emitter_Parameters(Workfunction[i], kT[i])
    
        energy[i], electron_count[i] = metal_emitter.Energy_Distribution()

    return energy, electron_count

def current_semiconductor_emitter(Field, Radius, Gamma, Ec, Ef, Eg, Temperature, me, mp):
    """
    Field [nm] - Electric field
    Radius [nm] - Emitter's tip radius
    Gamma [int] - Math parameter
    Ec [eV] - Bottom of the conduction band NB! on the backend make it into arrays
    Ef [eV] - Fermi level
    Eg [eV] - band gap
    Temperature [K] - Emitter's temperature
    me [kg] - electron effective mass
    mp [kg] - hole effective mass
    j_metal [A/nm^2] - Emitted current density

    For more info refer to GETELEC TABULATOR's documentation
    """
    tab = gt_tab.Tabulator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    semiconductor_emitter = gt_tab.Semiconductor_Emitter(tab)

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
    tab = gt_tab.Tabulator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    semiconductor_emitter = gt_tab.Semiconductor_Emitter(tab)

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
    tab = gt_tab.Tabulator()
    
    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    semiconductor_emitter = gt_tab.Semiconductor_Emitter(tab)

    energy_c = np.copy(Field)
    count_c = np.copy(Field)
    energy_v = np.copy(Field)
    count_v = np.copy(Field)
    m = np.ones(Field) * 9.1093837015e-31 

    for i in range(len(Field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m[i], me[i], mp[i])
        
        energy_c[i], count_c[i], energy_v[i], count_v[i] = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

    return energy_c, count_c, energy_v, count_v

#1st analyze iv button
def fit_data(xML, yML, F0, W0, R0, Gamma0, Temp0):
    fit_data = gt_mod.fitML(xML, yML, F0, W0, R0, Gamma0, Temp0)
    return fit_data

#notice i made approx=1
def plot_data(xfn, beta, W0, R0, Gamma0, Temp0, approx=1):
    plot_data = gt_mod.MLplot(xfn, beta, W0, R0, Gamma0, Temp0, approx)
    return plot_data