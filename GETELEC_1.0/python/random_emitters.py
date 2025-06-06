#! /usr/bin/python
import numpy as np
import getelec_mod as getelec_old
import matplotlib.pyplot as plt
import matplotlib as mb

font = 30
mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.5
fsize = (18,10)


Nemitters = 1000 # number of independent emitters

mean_height = 100 #nm
mean_radius = 3  #nm

sigma_height = 50
sigma_radius = .6

"""
mean_beta = 35
sigma_beta = 10

mean_area = 500
sigma_area = 150
"""


F0 = [2., 5., 15.]
R0 = [5., 20., 100.]
gamma0 = [1, 1.1, 10.]
Temp0 = [299., 300., 300.]
W0 = [4.49999, 4.5, 4.5]

heights = (np.random.randn(Nemitters)) * sigma_height + mean_height
radii = (np.random.randn(Nemitters)) * sigma_radius + mean_radius

areas = np.pi * radii**2
betas = abs(heights / radii)

Vexp, Iexp = np.loadtxt("IV_Yaroslava.txt", unpack=True)
Iexp *= 1.e-6
Efarexp = Vexp / 60000

Xfnfar = np.linspace(7,15, 32)

Efar = 1./Xfnfar
Itot = np.copy(Efarexp)

print "minF = %f, maxF = %f"%(min(Efarexp)* min(betas), max(Efarexp)* max(betas))


this = getelec_old.emission_create(W = 4.1, R = 100., gamma = 1.1, Temp = 300., approx = -1)

for i in range(len(Efarexp)):
    Itot[i] = 0.
    for j in range(len(betas)):
        this.F = betas[j] * Efarexp[i]
        this.R = radii[j]
        this.cur_dens()
        Itot[i] += this.Jem * areas[j]
        
plt.hist(betas)
plt.show()
        
xdata = 1./Efarexp
ydata = np.log(Itot)

# fit= gt.fitML(xdata,ydata, F0, W0, R0, gamma0, Temp0)
# popt = fit.x
# yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
# yshift = max(yopt) - max(ydata)
    
# print 'beta = %10.3e, W = %10.3f, R = %10.3f, gamma = %10.3f, Temp = %10.3f, sigmaAeff = %10.3e' \
        # % (popt[0], popt[1],  popt[2], popt[3], popt[4], 1e-9*np.exp(-yshift))        

# yth = np.exp(gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)  

sigma = max(Itot) / max(Iexp)

print "sigma = ", sigma
        
plt.semilogy(1/Efarexp,Itot / sigma)

plt.semilogy(1/Efarexp, Iexp, '.')


plt.show()
