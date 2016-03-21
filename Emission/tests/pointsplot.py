import numpy as np
import pylab as plt

#read data from file
x,y,reg, h= np.loadtxt("J-Fplot.csv",delimiter=',',unpack='True')
font=30

xf=x[np.where(reg==-1)]
xi=x[np.where(reg==0)]
xt=x[np.where(reg==1)]

yf=y[np.where(reg==-1)]
yi=y[np.where(reg==0)]
yt=y[np.where(reg==1)]

hf=h[np.where(reg==-1)]
hi=h[np.where(reg==0)]
ht=h[np.where(reg==1)]

#plot results
plt.figure(1)
plt.semilogy(xt,np.exp(yt),'r',label="ft")  # connect points with a blue line
plt.plot(xf,np.exp(yf),'b',label="ff")
plt.plot(xi,np.exp(yi),'k',label="fi")
plt.xlabel(r'$1/F (m/GV)$',fontsize=font)
plt.ylabel(r'$J (A/nm^2)$',fontsize=font)
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)

plt.figure(2)
plt.semilogy(xt,abs(ht),'r',label="ft")  # connect points with a blue line
plt.plot(xf,abs(hf),'b',label="ff")
plt.plot(xi,abs(hi),'k',label="fi")
plt.xlabel(r'$1/F (m/GV)$',fontsize=font)
plt.ylabel(r'$|heat| (W/nm^2)$',fontsize=font)
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)
plt.show()


