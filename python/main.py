import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import data as da

matplotlib.rcParams["font.family"] = "Serif"
matplotlib.rcParams["font.size"] = 35
matplotlib.rcParams["axes.labelsize"] = 35
matplotlib.rcParams["xtick.labelsize"] = 35
matplotlib.rcParams["ytick.labelsize"] = 35
matplotlib.rcParams["legend.fontsize"] = 35

fig = plt.figure(figsize=(20,15))
ax = fig.gca()
fig2 =  plt.figure(figsize=(20,15))
#ax2 = ax.twinx()
ax2 = fig2.gca()

ax.grid()
ax.set_xlabel(r"$1/F [nm/V]$")
ax.set_ylabel(r"$J [A/nm^2]$")
ax2.set_ylabel(r"$P_N [W/nm^2]$")
ax2.set_xlabel(r"$1/F [nm/V]$")
ax2.grid()
#for tl in ax2.get_yticklabels():
    #tl.set_color('r')
    
def separate(J,reg,sharp):
    Jfb = J[np.where(np.logical_and( reg == ord('F'), sharp == ord('B')))]
    Jfs = J[np.where(np.logical_and( reg == ord('F'), sharp == ord('S')))]
    Jtb = J[np.where(np.logical_and( reg == ord('T'), sharp == ord('B')))]
    Jts = J[np.where(np.logical_and( reg == ord('T'), sharp == ord('S')))]
    Jib = J[np.where(np.logical_and( reg == ord('I'), sharp == ord('B')))]
    Jis = J[np.where(np.logical_and( reg == ord('I'), sharp == ord('S')))]
    return (Jfb,Jfs,Jtb,Jts,Jib,Jis)
    
    
reload(da)
(J,F,Japp,heat,heatapp,reg,sharp,regap,sharpap) = da.dummy()

#ax2.set_ylim([1.e-5 * np.min(abs(heat)),np.max(abs(heat))])
#ax.set_ylim([np.min(J),1.e5*np.max(J)])

(Jfb,Jfs,Jtb,Jts,Jib,Jis) = separate(J,reg,sharp)
(Ffb,Ffs,Ftb,Fts,Fib,Fis) = separate(F,reg,sharp)
(hfb,hfs,htb,hts,hib,his) = separate(abs(heat),reg,sharp)


ax.semilogy(1/Ffb,Jfb,'c--',label="Field-Blunt", linewidth = 2)  # connect points with a blue line
ax.semilogy(1/Fib,Jib,'k--',label="Intermediate-Blunt", linewidth = 2) 
ax.semilogy(1/Ftb,Jtb,'m--',label="Thermal-Blunt", linewidth = 2)

ax.semilogy(1/Ffs,Jfs,'c-',label="Field-Sharp", linewidth = 2)  # connect points with a blue line
ax.semilogy(1/Fis,Jis,'k-',label="Intermediate-Sharp", linewidth = 2) 
ax.semilogy(1/Fts,Jts,'m-',label="Thermal-Sharp", linewidth = 2)

ax2.semilogy(1/Ffb,hfb,'b--',label=r"P_FN-Field-Blunt", linewidth = 2)  # connect points with a blue line
ax2.semilogy(1/Fib,hib,'k--',label="Intermediate-Blunt", linewidth = 2) 
ax2.semilogy(1/Ftb,htb,'r--',label="Thermal-Blunt", linewidth = 2)

ax2.semilogy(1/Ffs,hfs,'b-',label="Field-Sharp", linewidth = 2)  # connect points with a blue line
ax2.semilogy(1/Fis,his,'k-',label="Intermediate-Sharp", linewidth = 2) 
ax2.semilogy(1/Fts,hts,'r-',label="Thermal-Sharp", linewidth = 2)

ax.semilogy(1/F,Japp,'k:',label="Field-Blunt", linewidth = 2)
ax2.semilogy(1/F,abs(heatapp),'k:',label="Field-Blunt", linewidth = 2)

#ax.legend(loc="best")
#ax2.legend(loc="best")
 
plt.show()
