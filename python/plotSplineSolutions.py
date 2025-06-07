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

y1 = np.loadtxt('problematicSplineSolution.dat')


plt.plot(y1[:,0], y1[:,1], ".", label="Re(s')", color = colors[0])
plt.plot(y1[:,0], y1[:,4], "--", label="spline Re(s')", color = colors[0])

plt.plot(y1[:,0], y1[:,2], ".", label="Im(s')", color = colors[1])
plt.plot(y1[:,0], y1[:,5], "--", label="spline Im(s')", color = colors[1])

plt.plot(y1[:,0], y1[:,3], "," , label="Re(s)", color = colors[2])
plt.plot(y1[:,0], y1[:,6], "--", label="spline Re(s)", color = colors[2])


nodes = np.loadtxt('problematicSplineNodes.dat')
dx = 0.1
plt.plot(nodes[:,0], nodes[:,1], 'o', label="Re(s') nodes")
plt.plot(nodes[:,0], nodes[:,3], 'o', label="Im(s') nodes")
plt.plot(nodes[:,0], nodes[:,5], 'o', label="Re(s) nodes")

# plt.quiver(nodes[:,0], nodes[:,1], dx * nodes[:,2] / nodes[:,2], nodes[:,2] * dx, angles = "xy", scale_units='xy', scale=1, color=colors[0])
# plt.quiver(nodes[:,0], nodes[:,3], dx* nodes[:,2] / nodes[:,2], nodes[:,4] * dx, angles = "xy", scale_units='xy', scale=1, color=colors[1])
# plt.quiver(nodes[:,0], nodes[:,5], dx* nodes[:,2] / nodes[:,2], nodes[:,6] * dx, angles = "xy", scale_units='xy', scale=1, color=colors[2])



plt.legend()
plt.show()