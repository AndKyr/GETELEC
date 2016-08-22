import matplotlib
import numpy as np
import matplotlib.pyplot as plt

font = 32
fontstyle = 'Serif'
lw = 2

matplotlib.rcParams["font.family"] = fontstyle
matplotlib.rcParams["font.size"] = font
matplotlib.rcParams["axes.labelsize"] = font
matplotlib.rcParams["xtick.labelsize"] = font
matplotlib.rcParams["ytick.labelsize"] = font
matplotlib.rcParams["legend.fontsize"] = font

fig = plt.figure(figsize=(20,15))
ax = fig.gca()
#ax.grid()
ax.set_xlabel(r"$x [nm]$")
ax.set_ylabel(r"$V \ [ \ V \ ]$")
#ax.set_title(r"FN-plot test")

Nplot = 100
xmax = 10.
ls = ['b-', 'r--', 'm:', 'k-.']
gammas = [10.,100.]
Rs = [2., 5.]

def Vmod(x,F=5.,R=5.,g=10.):
    return (F*R*x*(g-1.)+F*x**2)/(g*x + R*(g-1))

x = np.linspace(0.,xmax, Nplot)
i=0
for R in Rs:
    for g in gammas:
        y = Vmod(x, 5., R, g)
        lab = r'$R=' + str(int(R)) + r'nm$, $\gamma=' + str(int(g)) + r'$'
        ax.plot(x,y,ls[i],linewidth=lw,label=lab)
        i = i+1

ax.legend(loc="best")
plt.show()



