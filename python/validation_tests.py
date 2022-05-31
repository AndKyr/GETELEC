from pkgutil import get_data
import getelec_tabulator as gtab
import getelec_mod as gt
import numpy as np
import matplotlib as mb
import matplotlib.pyplot as plt
import datetime

# region One data point calculation routine
    # This routine calculates the current density from semiconductors (two methods) and metals, as well as the plotting of the energy distributions
"""  
Npoly = 5
NGi = 512
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gtab.Tabulator()

Workfunction = 4.5
Ec = 4.4195
Ef = 4.75
Eg = 0.661
m = 9.1093837015e-31 
me = 1.59*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.33*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
Temp = 300.
kT = kBoltz * Temp

metal_emitter = gtab.Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

metal_emitter.Define_Emitter_Parameters(Workfunction, kT)

j_metal = metal_emitter.Current_Density()
pn_metal = metal_emitter.Nottingham_Heat()
energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()


semiconductor_emitter = gtab.Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT, m, me, mp)

j_c, j_v, j_total = semiconductor_emitter.Current_Density_from_Semiconductors()

pn_c, pn_v, pn_total = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

j_c2, j_v2, j_total2 = semiconductor_emitter.Distributed_Current_Density_from_Semiconductors()

c_abs_error = abs(j_c2-j_c)
c_rel_error = c_abs_error/j_c

v_abs_error = abs(j_v2-j_v)
v_rel_error = v_abs_error/j_v

total_abs_error = abs(j_total2-j_total)
total_rel_error = total_abs_error/j_total

print("current_c", j_c, j_c2, "relative error:", c_rel_error)
print("current_v", j_v, j_v2, "relative error:", v_rel_error)
print("current_semi", j_total, j_total2, "relative error:", total_rel_error)
print("current_metal", j_metal)


font = 60
x = 40
y = 17


mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font+5
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 7

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_space_metal-Workfunction, distribution_metal/max(distribution_metal), color = "steelblue", label = "electrons from metal") 
plt.plot(energy_c-Ef+Eg/2, distribution_c/max(distribution_c), color = "orange", label = "electrons from conduction band") 
plt.plot(energy_v-Ef+Eg/2, distribution_v/max(distribution_v), color = "green", label = "electrons from valence band") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_{eN}$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Normalised electron energy distribution")
plt.savefig("Normalised electron energy distribution.svg")
plt.savefig("Normalised electron energy distribution.png")

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_space_metal-Workfunction, distribution_metal, color = "steelblue", label = "electrons from metal") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Electron energy distribution from a metal")
plt.savefig("Electron energy distribution from a metal.svg")
plt.savefig("Electron energy distribution from a metal.png")

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_c-Ef+Eg/2, distribution_c, color = "orange", label = "electrons from conduction band") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Electron energy distribution from conduction band")
plt.savefig("Electron energy distribution from conduction band.svg")
plt.savefig("Electron energy distribution from conduction band.png")

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_v-Ef+Eg/2, distribution_v, color = "green", label = "electrons from valence band") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Electron energy distribution from valence band")
plt.savefig("Electron energy distribution from valence band.svg")
plt.savefig("Electron energy distribution from valence band.png")
"""
# endregion

# region Multiple data point calculation routine - metal

#Npoly = 5
#NGi = 128
#zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gtab.Tabulator()

Fmax = 1/tab.Finv[0]
Fmin = 1/tab.Finv[-1]
Rmax = 1/tab.Rinv[0]
Rmin = 1/tab.Rinv[-1]
gammax = 1/tab.gaminv[0]
gammin = 1/tab.gaminv[-1]

Np = 80960

Fi = np.random.rand(Np) * (Fmax - Fmin) + Fmin
Ri = np.random.rand(Np) * (Rmax - Rmin) + Rmin
gami = np.random.rand(Np) * (gammax - gammin) + gammin
Wi = np.random.rand(Np) * (7.5 - 2.5) + 2.5
Ti = np.random.rand(Np) * (3000 - 100) + 100
kT = Ti * kBoltz
Ji = np.copy(Fi)
Pi = np.copy(Fi)
Jget = np.copy(Ji)
Pget = np.copy(Ji)

metal_emitter = gtab.Metal_Emitter(tab)

print("calculating from tabulator")
tab_start = datetime.datetime.now()
for i in range(len(Fi)):
    metal_emitter.emitter.Define_Barrier_Parameters(Fi[i], Ri[i], gami[i])
    metal_emitter.emitter.Interpolate_Gammow()
    metal_emitter.Define_Emitter_Parameters(Wi[i], kT[i])
    Ji[i] = metal_emitter.Current_Density()
    Pi[i] = metal_emitter.Nottingham_Heat()
tab_end = datetime.datetime.now()

print("calculating from getelec")
get_start = datetime.datetime.now()
em = gt.emission_create(approx=2)
for i in range(len(Fi)):   
    em.F = Fi[i]
    em.W = Wi[i]
    em.Temp = Ti[i]
    em.gamma = gami[i]
    em.R = Ri[i]
    em.cur_dens()
    Jget[i] = em.Jem
    Pget[i] = em.heat
get_end = datetime.datetime.now()


abserr = abs(Ji - Jget)
relerr = abserr / Jget
bad = np.where(np.logical_and(relerr > 0.5, abserr > 1.e-25))[0]

Pn_abserr = abs(Pi - Pget)
Pn_relerr = Pn_abserr / Pget
Pbad = np.where(np.logical_and(Pn_relerr > 0.5, Pn_abserr > 1.e-25))[0]

print("bad = ", bad)
print("rms error in J = ", np.sqrt(np.mean(relerr[abserr > 1.e-25]**2)))
print("rms error in Pn = ", np.sqrt(np.mean(Pn_relerr[Pn_abserr > 1.e-25]**2)))
print("getelec running time =", get_end-get_start)
print("tabulat running time =", tab_end-tab_start)


#for i in bad:
#    print("Jget, Ji : ", Jget[i], Ji[i])
#    emit.set(Fi[i], Ri[i], gami[i])
#    emit.interpolate()
#    emit.get_lims(Wi[i], kT[i])
#    emit.integrate_quad(Wi[i], kT[i])
#    emit.integrate_quad_Nottingham(W[i], kT[i])

fig = plt.figure(figsize=(16,6))
plt.loglog(Ji, Jget, '.')
plt.loglog(abs(Pi), abs(Pget), '.')
plt.loglog([1.e-50, 1.], [1.e-50, 1.])
plt.grid()
plt.title("J comparison")
plt.savefig("J comparison.png")
#plt.show()

fig2 = plt.figure(figsize=(16,6))
plt.loglog(Pi, Pget, '.')
plt.loglog(abs(Pi), abs(Pget), '.')
plt.loglog([1.e-50, 1.], [1.e-50, 1.])
plt.grid()
plt.title("Pn comparison")
plt.savefig("Pn comparison.png")
#plt.show()

#endregion

# region Multiple data point calculation routine - semiconductor
"""
Npoly = 5
NGi = 128
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gt.Tabulator()

Ef = 4.75
Eg = 0.661
Temp = 300.
m = 9.1093837015e-31 
me = 1.59*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.33*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
kT = kBoltz * Temp

semiconductor_emitter = gt.Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

metal_emitter = gt.Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

resolution = 100

j_c = np.zeros(resolution) 
j_v = np.zeros(resolution) 
j_total = np.zeros(resolution) 
j_c2 = np.zeros(resolution) 
j_v2 = np.zeros(resolution) 
j_total2 = np.zeros(resolution) 
pn_c = np.zeros(resolution) 
pn_v = np.zeros(resolution) 
pn_total = np.zeros(resolution) 
Bottom_Ec = np.zeros(resolution) 
j_metal = np.zeros(resolution) 
pn_metal = np.zeros(resolution)
c_rel_error = np.zeros(resolution)
v_rel_error = np.zeros(resolution)
total_rel_error = np.zeros(resolution)

for i in range(resolution):
    Ec = i*(0.5/resolution)+4.5

    semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT, m, me, mp)

    metal_emitter.Define_Emitter_Parameters(Ec, kT)


    j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()
    
    j_c2[i], j_v2[i], j_total2[i] = semiconductor_emitter.Distributed_Current_Density_from_Semiconductors()
    
    pn_c[i], pn_v[i], pn_total[i] = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()
    
    energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()
    
    
    j_metal[i] = metal_emitter.Current_Density()
    
    pn_metal[i] = metal_emitter.Nottingham_Heat()
    
    energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()
    
    
    Bottom_Ec[i] = -Ec
    
    
    c_abs_error = abs(j_c2[i]-j_c[i])
    c_rel_error[i] = c_abs_error/j_c[i]

    v_abs_error = abs(j_v2[i]-j_v[i])
    v_rel_error[i] = v_abs_error/j_v[i]

    total_abs_error = abs(j_total2[i]-j_total[i])
    total_rel_error[i] = total_abs_error/j_total[i]
    
    
    title = "Normalised energy distrution, Ec = " + str(Ec)
    
    save_png = "Normalised energy distrution, Ec = " + str(Ec) + ".png"
    
    #fig = plt.figure(figsize=(16,6))
    #plt.plot(energy_space_metal, distribution_metal/max(distribution_metal))
    #plt.plot(energy_c, distribution_c/max(distribution_c))
    #plt.plot(energy_v, distribution_v/max(distribution_v))
    #plt.grid("True")
    #plt.title(title)
    #plt.savefig(save_png)

c_mean = np.mean(c_rel_error)
c_std = np.std(c_rel_error)
c_max = max(c_rel_error)

v_mean = np.mean(v_rel_error)
v_std = np.std(v_rel_error)
v_max = max(v_rel_error)


total_mean = np.mean(total_rel_error)
total_std = np.std(total_rel_error)
total_max = max(total_rel_error)

font = 60
x = 40
y = 17

mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font+5
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 5

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.semilogy(Bottom_Ec, j_c, color = "steelblue", label = "electrons from conduction band")
plt.semilogy(Bottom_Ec, j_v, color = "green", label = "electrons from valence band")
plt.semilogy(Bottom_Ec, j_total, color = "orange", label = "total emitted electrons")
#plt.plot(Bottom_Ec, j_metal)
#plt.yscale("symlog", linscale=-1)
#plt.legend()
plt.grid("True")
plt.xlabel("Bottom of conduction band (eV)")
plt.ylabel("$J_{e}$ (A$nm^{-2}$)")
plt.title("Detail")
plt.savefig("Detail.png")
plt.savefig("Detail.svg")

file = open("Field_Emission_Data.txt", "w+")

file.write("Simulation temperature = %f\r\n\r\n" % Temp)

for i in range(len(Bottom_Ec)):
    file.write("\r\nSemiconductor Ec = %d" % (Bottom_Ec[i]))
    file.write(" Metal Ef = %d\r\n" % (Bottom_Ec[i]))
    file.write("    Current from semiconductor %e" % (j_total[i]))
    file.write("    Current from Metal = %e\r\n" % (j_metal[i]))
    file.write("        Semiconductor current relative error = %e\r\n" % (total_rel_error[i]))
    file.write("    Heat from semiconductor = %e" % (pn_total[i]))
    file.write("    Heat from metal = %e\r\n" % (pn_metal[i]))
    
file.write("\r\nMax relative error = %e\r\n" % (total_max))
  
file.close()
"""
# endregion