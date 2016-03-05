import numpy as np
import matplotlib.pyplot as plt

#read data from file
x,y,reg, h= np.loadtxt("J-Fplot.csv",delimiter=',',unpack='True')	

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
plt.plot(xt,yt,'r',label="ft")  # connect points with a blue line
plt.plot(xf,yf,'b',label="ff")
plt.plot(xi,yi,'k',label="fi")

plt.figure(2)
plt.semilogy(xt,abs(ht),'r',label="ft")  # connect points with a blue line
plt.plot(xf,abs(hf),'b',label="ff")
plt.plot(xi,abs(hi),'k',label="fi")
plt.show()


