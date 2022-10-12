import getelec_tabulator as gt_tab
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mb

def spectra_semiconductor_emitter(Field, Radius, Gamma, Ec, Ef, Eg, Temperature):

    kBoltz = 8.6173324e-5 
    kT = kBoltz * Temperature

    m = 9.1093837015e-31 
    me = 1.59*m
    mp = 0.33*m

    tab = gt_tab.Tabulator()
    semiconductor_emitter = gt_tab.Semiconductor_Emitter(tab)

    energies = np.zeros([len(Field)*2,512])
    counts = np.zeros([len(Field)*2,512])


    for i in range(len(Field)):

        semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        semiconductor_emitter.emitter.Interpolate_Gammow()

        semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[i], Ef[i], Eg[i], kT[i], m, me, mp)
        
        energy_c, count_c, energy_v, count_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()    

        energies[i,] = energy_c
        energies[i+int((len(Field)/2)),] = energy_v
        counts[i,] = count_c
        counts[i+int((len(Field)/2)),] = count_v

    energy = np.linspace(np.min(energies), np.max(energies), 512)
    electrons = np.zeros(len(energy))

    for i in range(len(Field)):
            f = interpolate.interp1d(energies[i], counts[i], kind="quadratic", bounds_error=False, fill_value=0.)
            electrons = electrons + f(energy)
    

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
    plt.plot(energy, electrons, '-',color = "orange", label = "total emitted electrons")
    plt.legend()
    plt.grid("True")
    plt.xlabel("Energy (eV)")
    plt.ylabel("$J_{eN}$(E) (A$nm^{-2}$$eV^{-1}$)")
    plt.title("Electron spectrum_SIMULATION@600V")
    plt.savefig("Electron spectrum_SIMULATION@600V.png")
    plt.show()

    return energy, electrons


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

data = np.loadtxt("/home/salva/Documents/getelec_priv/python/spectra")

field = np.zeros(len(data))
radius = np.zeros(len(data))
gamma = np.zeros(len(data))
ec = np.zeros(len(data))
ef = np.zeros(len(data))
eg = np.zeros(len(data))
temp = np.zeros(len(data))

for i in range(len(data)):
    field[i] = data[i,0]
    radius[i] = data[i,1]
    gamma[i] = data[i,2]
    ec[i] = data[i,3]
    ef[i] = data[i,4]
    eg[i] = data[i,5]
    temp[i] =data[i,6]

spectra_semiconductor_emitter(field, radius, gamma, ec, ef, eg, temp)


"""TO DO

1) Implement routine to calculate and show half-width
2) Implement routine to calculate emitter area
3) Implemetn routine to visualise dynamics of emitting side
4) Implement routine to calculate band bending slope and comparison with D(E)
5) Calculate e- distribution of metal, semiconductor and non-bending semiconductor emitting same I

"""