import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func(x, F, R, beta):
	nom=(F*R*(beta-1.)*x + F*x*x)
	denom=(beta*x+R*(beta-1.))
	denom[denom==0.]=1e-100
	
	return nom/denom

xdata = np.linspace(0, 4, 50)
y = func(xdata, 2.5, 1.3, 10.)
ydata = y + 0.2 * np.random.normal(size=len(xdata))

popt, pcov = curve_fit(func, xdata, ydata)

print popt
print pcov
plt.plot(xdata,ydata,'*')
plt.plot(xdata,y)
plt.plot(xdata,func(xdata,popt[0],popt[1],popt[2]))
plt.show()ged
