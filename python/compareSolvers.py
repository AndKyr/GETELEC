import gamowCalculator as gc
import time
import matplotlib.pyplot as plt
import numpy as np

import getelec_wrap as gt
font = 25
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
figureSize = [10,10]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

calculator = gc.GamowCalculator(solverType="IVP", XCdataFile="", minimumPotential=15.)

start_time = time.time()
calculator.calculateGamow(4.5)
# calculator.calculateGamowCurve(64, minBarrierDepth=-1., maxGamow=50)
# ax1, ax2 = calculator.solver.plotWaveFunction()

calculator.solver.getSolutionVector()
fig, (ax1,ax2, ax3) = plt.subplots(3,1, sharex=True, figsize = figureSize, tight_layout=True)
ax1.plot(calculator.solver.xVector, calculator.solver.potentialFunction(calculator.solver.xVector), label=r"$U(z)$", color = colors[6])
ax1.plot([calculator.solver.xVector[0], calculator.solver.xVector[-1]], [calculator.solver.energy, calculator.solver.energy], label=r"$E$", color = colors[5])
ax1.grid()
ax1.legend()
ax1.set_ylabel(r"$\textrm{Energy [eV]}$")

ax2.plot(calculator.solver.xVector, np.real(calculator.solver.solutionVector), ".-", label=r"$\Re(\Psi) \textrm{ (%d steps)}$"%len(calculator.solver.solutionVector), color = colors[1])
ax2.plot(calculator.solver.xVector, np.imag(calculator.solver.solutionVector), ".-", label=r"$\Im(\Psi)$", color = colors[0])

ax2.legend()
ax2.grid()
ax2.set_ylabel(r"$\Psi$")
ax3.set_xlabel(r"$z \textrm{[nm]}$")
ax2.set_ylim(-2.5, 10)

emitter = gt.GetelecInterface()
emitter.calculateTransmissionCoefficients(np.array([-4.5]), np.array([12.]), writeFiles=True)

file = "odeSolution_i_0_E_-4.500000.dat"

y1 = np.loadtxt(file)

ax3.plot(y1[:,0], y1[:,1], '.-', label=r"$\Re(s')$", color = colors[2])
ax3.plot(y1[:,0], y1[:,2], '.-', label=r"$\Im(s')\textrm{ (%d steps)}$"%y1.shape[0], color = colors[3])
ax3.plot(y1[:,0], y1[:,3], '.-', label=r"$\Re(s)$", color = colors[4])


ax3.legend()
ax3.grid()
ax3.set_ylabel(r"$s$")

plt.savefig("solverComparison.pdf")
plt.show()

