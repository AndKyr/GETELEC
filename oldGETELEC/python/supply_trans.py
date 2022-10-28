#! /usr/bin/python
import numpy as np
import getelec_mod as getelec_old

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mb

font = 30
# mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.5
fsize = (15,18)


def Nsupply(Ez, kT):
    
    if (Ez > 10 * kT):
        return np.exp(-Ez/kT)
    elif(Ez < -10 * kT):
        return   - Ez / kT
    else:
        return  np.log(1. + np.exp(-Ez/kT))
    

kBoltz = 8.6173324e-5

T = [300., 1000., 2000., 3000., 3000., 3000., 3000., 3000.]
F = [5., 5., 5., 5., 4., 3., 2., 1.]




Ez = np.linspace(-3., 5., 256)
Ns = np.copy(Ez)
sfile = 'output/spectra.csv'
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

Zs = 1.618e-4

for k in range(len(F)):
    fig = plt.figure(figsize=fsize)


    ax, ax2 = fig.subplots(2,1, sharex = True)

    # ax = fig.gca()
    # ax2 = ax.twinx()



    for j in range(len(Ez)):
        Ns[j] = Nsupply(Ez[j], kBoltz * T[k])  
    ax.semilogy(Ez, Ns, c = colors[0], ls = '--', label = r'$L(E_z)$, T = %d K'%T[k])

    # ax2.tick_params(axis='y', labelcolor='red')
    # ax2.ticklabel_format(axis='y', style='sci')
    # ax2.set_ylabel(r'$10^4 \times N(E_z)$ [A nm$^{-2}$ eV$^{-1}$]', color = 'red')


    this = getelec_old.emission_create(R = 5000.)
    this.approx = 2
    ax2.set_xlabel(r"$E_z-E_F$ [eV]")
    ax.set_ylabel(r"$D(E_z), L(E_z) $")
    ax2.set_ylabel(r"$j(E_z) / j_{max}$")

    kk = 0



       
    this.F = F[k] 

    this.Temp = T[k]
    this.cur_dens()
    E,JE,J, Nj, G = np.loadtxt(sfile, delimiter = ',' ,unpack=True)
    ax.semilogy(E, 1./(1. + np.exp(G)), c = colors[0], label = r'$D(E_z)$, F = %g V/nm'%(F[k]))
    ax2.plot(E, Nj / max(Nj),c = colors[0])
    imax = np.argmax(Nj)
    ax.set_title(r'$n \equiv 1/\beta k_B T$ = %.1f'%(1/(-kBoltz * this.Temp * np.gradient(G)[imax] / np.gradient(E)[imax])))
    ax2.set_title(r'$j_{max} =$ %.1e A nm$^{-2}$ eV$^{-1}$'%(max(Nj) * kBoltz * this.Temp* Zs))
    this.Temp = 2500.
    this.cur_dens()
    E,JE,J, Nj, G = np.loadtxt(sfile, delimiter = ',' ,unpack=True)
    ax.semilogy(E, 1./(1. + np.exp(G)), c = colors[0])
        # ax2.plot(E, Nj,c = colors[k+3], label = 'F = %d V/nm, T = %d K'%(F[k], this.Temp))
        

    ax.set_ylim([min(1./(1. + np.exp(G))), max(Ns)])  
    # ax2.legend()
    ax2.grid()

    ax.legend()
    ax.grid()

    fig.tight_layout()

    plt.savefig('/home/ad/home/k/kyritsak/Documents/sup_trans_f=%dT=%d.png'%(F[k],T[k]))
    plt.close('all')
    # print '/home/ad/home/k/kyritsak/Documents/sup_trans_f=%dT=%d.png'%(F[k],T[k])
