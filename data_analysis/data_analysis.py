from array import array
import getelec_tabulator as getelec_new
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mb

def Load_Data(path_to_data:str):

    data = np.loadtxt(path_to_data)

    field = data[:,0]
    radius = data[:,1]
    gamma = data[:,2]
    ec = data[:,3]
    ef = data[:,4]
    eg = data[:,5]
    temp = data[:,6]

    return field, radius, gamma, ec, ef, eg, temp

def FWHW(x_axis:array, y_axis:array, frac:int):

    d = y_axis - (max(y_axis) / frac) 
    indexes = np.where(d > 0)[0] 

    return abs(x_axis[indexes[-1]] - x_axis[indexes[0]])

def Spectrum_Semiconductor_Emitter(field:array, radius:array, gamma:array, ec:array, ef:array, eg:array, temperature:array):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * temperature

    m = 9.1093837015e-31 
    me = 1.59*m
    mp = 0.33*m

    tab = getelec_new.Interpolator()
    semiconductor_emitter = getelec_new.Semiconductor_Emitter(tab)

    energies = np.zeros([len(field)*2,512])
    counts = np.zeros([len(field)*2,512])


    for i in range(len(field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(field[i], radius[i], gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(ec[i], ef[i], eg[i], kT[i], m, me, mp)
        
        energy_c, count_c, energy_v, count_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()    

        energies[i,] = energy_c
        energies[i+int((len(field)/2)),] = energy_v
        counts[i,] = count_c
        counts[i+int((len(field)/2)),] = count_v

    energy = np.linspace(np.min(energies), np.max(energies), 512)
    electrons = np.zeros(len(energy))

    for i in range(len(field)):
            f = interpolate.interp1d(energies[i], counts[i], kind="quadratic", bounds_error=False, fill_value=0.)
            electrons = electrons + f(energy)
    
    fwhm = FWHW(energy,electrons,2)
    
    font = 60
    x = 40
    y = 17
   
    t = "FWHM = "+str(np.round(fwhm,5))+" eV"

    mb.rcParams["font.size"] = font
    mb.rcParams["axes.labelsize"] = font
    mb.rcParams["xtick.labelsize"] = font
    mb.rcParams["ytick.labelsize"] = font+5
    mb.rcParams["legend.fontsize"] = font/1.5
    mb.rcParams["lines.linewidth"] = 7
    fig, ax = plt.subplots(figsize=(x,y))
    ax.ticklabel_format(scilimits=[-1,1])
    plt.plot(energy, electrons, '-',color = "orange", label = t)
    plt.legend()
    plt.grid("True")
    plt.xlabel("Energy (eV)")
    plt.ylabel("$J_{eN}$(E) (A$nm^{-2}$$eV^{-1}$)")
    plt.title("Electron spectrum_SIMULATION@600V")
    plt.savefig("Electron spectrum_SIMULATION@600V.png")
    plt.show()

    return energy, electrons

def Spectra_Metal_Emitter(field:array, radius:array, gamma:array, workfunction:array, temperature:array):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * temperature

    tab = getelec_new.Interpolator()
    metal_emitter = getelec_new.Metal_Emitter(tab)

    energies = np.zeros([len(field)*2,512])
    counts = np.zeros([len(field)*2,512])


    for i in range(len(field)):

        metal_emitter._emitter.Define_Barrier_Parameters(field[i], radius[i], gamma[i])
        metal_emitter._emitter.Interpolate_Gammow()

        metal_emitter.Define_Metal_Emitter_Parameters(workfunction[i], kT[i])
        
        energy, count = metal_emitter.Energy_Distribution_for_Metals()    

        energies[i,] = energy
        counts[i,] = count

    energy = np.linspace(np.min(energies), np.max(energies), 512)
    electrons = np.zeros(len(energy))

    for i in range(len(field)):
            f = interpolate.interp1d(energies[i], counts[i], kind="quadratic", bounds_error=False, fill_value=0.)
            electrons = electrons + f(energy)
    
    fwhm = FWHW(energy,electrons,2)
    
    font = 60
    x = 40
    y = 17
   
    t = "FWHM = "+str(np.round(fwhm,5))+" eV"

    mb.rcParams["font.size"] = font
    mb.rcParams["axes.labelsize"] = font
    mb.rcParams["xtick.labelsize"] = font
    mb.rcParams["ytick.labelsize"] = font+5
    mb.rcParams["legend.fontsize"] = font/1.5
    mb.rcParams["lines.linewidth"] = 7
    fig, ax = plt.subplots(figsize=(x,y))
    ax.ticklabel_format(scilimits=[-1,1])
    plt.plot(energy, electrons, '-',color = "orange", label = t)
    plt.legend()
    plt.grid("True")
    plt.xlabel("Energy (eV)")
    plt.ylabel("$J_{eN}$(E) (A$nm^{-2}$$eV^{-1}$)")
    plt.title("Electron spectrum_SIMULATION@600V")
    plt.savefig("Electron spectrum_SIMULATION@600V.png")
    plt.show()

    return energy, electrons

def Emitting_Area_Forbes(emitted_current:float, max_current_density:float):
    arc_length = emitted_current/max_current_density
    return 2*np.pi*arc_length**2

def Emitting_Area_Salva(current_density:array, arc_length:array):
    indexes = np.where(current_density>=1e-15)[0]

    return np.sum(2*np.pi*arc_length[indexes[:-1]]**2)

def Emitting_Site_Dynamics(potential:array, current_density:array):
    
    x_axis = 1/potential
    y_axis = np.log10(abs(current_density)/potential**2)

    font = 60
    x = 40
    y = 17

    mb.rcParams["font.size"] = font
    mb.rcParams["axes.labelsize"] = font
    mb.rcParams["xtick.labelsize"] = font
    mb.rcParams["ytick.labelsize"] = font+5
    mb.rcParams["legend.fontsize"] = font/1.5
    mb.rcParams["lines.linewidth"] = 7
    fig, ax = plt.subplots(figsize=(x,y))
    ax.ticklabel_format(scilimits=[-1,1])
    plt.plot(x_axis, y_axis, '-')
    plt.legend()
    plt.grid("True")
    plt.xlabel("1/V ($V^{-1}$)")
    plt.ylabel("$J_{eN}$(E) (A$nm^{-2}$)")
    plt.title("Emitting site dynamics")
    plt.savefig("Emitting site dynamics.png")
    plt.show()

    return True

path = "/home/salva/Documents/getelec_priv/python/spectra"

field, radius, gamma, ec, ef, eg, temp = Load_Data(path)

current_density = getelec_new.current_semiconductor_emitter(field, radius, gamma, ec, ef, eg, temp)


"""field = np.array([5,5.05,5.1,5.15,5.2,5.3])
radius = np.ones(len(field))*20
gamma =  np.ones(len(field))*10
ec = np.array([4.05,4.06,4.07,4.08,4.09,4.1])
ef = np.ones(len(field))*4.5
eg = np.ones(len(field))*1.12
temp = np.ones(len(field))*300
#me = 1.59
#mp = 0.33

spectra_semiconductor_emitter(field, radius, gamma, ec, ef, eg, temp)"""