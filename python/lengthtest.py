#!/usr/bin/python3

"""This plots the spectra as outputed in spectra.csv from getelec. Spectra has to be 
T in the GetelecPar.in""" 
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mb
import os

font = 30
mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5




fig1 = plt.figure()
ax1 = fig1.gca()
ax1.set_xlabel(r"$R [nm]$")
ax1.set_ylabel(r"$L [nm]$")
ax1.set_title(r'Barrier length (L) - Radius of curvature (R)')


R = np.linspace(5.,50., 100)
L = np.copy(R)
Lex = np.copy(R)



temp = 300.

kT = 8.6173324e-5 * temp

os.system('cp in/par2.in in/GetelecPar.in')

for i in range(len(L)):
    appstr = os.popen('./bin/current.exe 5. 4.5 300 %f | grep "L = " '%(R[i])).read().split()
    #os.system('./bin/current.exe %f 4.5 3000.'%field)
    L[i] = float(appstr[2])


    
os.system('cp in/par1.in in/GetelecPar.in')
    
for i in range(len(L)):
    appstr = os.popen('./bin/current.exe 5. 4.5 300 %f | grep "L = " '%(R[i])).read().split()
    #os.system('./bin/current.exe %f 4.5 3000.'%field)
    Lex[i] = float(appstr[2])
    

appstr = os.popen('./bin/current.exe 5. 4.5 300 3000. | grep "L = " '%(R[i])).read().split()
Lsn = float(appstr[2])

    
ax1.plot(R,L,'b-', label = 'Approximate - KX')
ax1.plot(R,Lex,'r-', label = 'Numerical')
ax1.plot(np.array([R[0], R[-1]]), np.array([Lsn, Lsn]), 'k-', label = 'Approximate - SN')

ax1.legend(loc = 'best')    

plt.show()
