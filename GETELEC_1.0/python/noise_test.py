#! /usr/bin/python3
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
fsize = (18,10)

# cMap = plt.cm.copper
Nrand = 2048
em = getelec_old.emission_create(W = 4.5, R = 4., gamma = 3.)

def J_noise_average(F, F_noise):
	fi = F + np.random.randn(Nrand) * F_noise
	Ji = np.copy(fi)
	for i in range(len(fi)):
		em.F = fi[i]
		em.cur_dens()
		Ji[i] = em.Jem
	
	return np.mean(Ji), np.std(Ji)

Npoints = 32

fnoise = 0.03
Fshift = 0.05

Fi = np.linspace(5., 7.5, Npoints)
Ji = np.copy(Fi)
J2i = np.copy(Fi)
J2noise = np.copy(Fi)
Jnoise = np.copy(Fi)
J_std = np.copy(Fi)
J2_std = np.copy(Fi)

for i in range(len(Fi)):
	em.F = Fi[i]
	em.cur_dens()
	Ji[i] = em.Jem

	em.F = Fi[i] + Fshift
	em.cur_dens()
	J2i[i] = em.Jem

	Jnoise[i], J_std[i] = J_noise_average(Fi[i], fnoise * Fi[i])
	J2noise[i], J2_std[i] = J_noise_average(Fi[i] + \
		Fshift, fnoise * (Fi[i] + Fshift))

# plt.plot(Fi, Ji, label = 'Jnormal')
# plt.errorbar(Fi, Jnoise,yerr = J_std, label = 'Javg')

# plt.plot(Fi, J2i, '--', label = 'JShifted')
# plt.errorbar(Fi, J2noise,yerr = J2_std, fmt = '--', label = 'J2avg_shifted')

plt.plot(Fi, J2i/Ji)
plt.plot(Fi,J2noise / Jnoise)

ax = plt.gca()
#ax.set_yscale('log')

plt.xlabel(r"$F[GV/m]$")
plt.ylabel(r"$J [A/nm^2]$")
plt.grid()
plt.legend()



plt.show()
