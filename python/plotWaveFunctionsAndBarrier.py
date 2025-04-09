import numpy as np
import matplotlib.pyplot as plt
import matplotlib


font = 20
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["font.family"] = "a"
matplotlib.rcParams["font.size"] = font
matplotlib.rcParams["axes.labelsize"] = font
matplotlib.rcParams["xtick.labelsize"] = font
matplotlib.rcParams["ytick.labelsize"] = font
matplotlib.rcParams["legend.fontsize"] = font
# matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['figure.figsize'] = 12, 8
matplotlib.rcParams['lines.markersize'] = 10
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']



# Load the wave function data
x, barrier = np.loadtxt('nastyCaseBarrier.dat', unpack=True)

fig, (ax1, ax2) = plt.subplots(2,1)


ax1.plot(x, barrier, label="V(x) - E", color=colors[0])
ax1.plot([x[0], x[-1]], [0, 0], 'k--')

solutions = np.loadtxt('nastyCaseSolution.dat')

for i in range(1, 7):
    ax2.plot(solutions[:,0], solutions[:, i], label=f"Wave function {i}", color=colors[i])

solutions = np.loadtxt('nastyCaseSolution_neighborEnergy.dat')

for i in range(1, 7):
    ax2.plot(solutions[:,0], solutions[:, i],'--', label=f"Wave function {i}", color=colors[i])

ax2.legend()
ax2.set_xlabel(r"$x$ (nm)")
plt.show()