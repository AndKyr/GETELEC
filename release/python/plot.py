import numpy as np
import matplotlib.pyplot as plt


x = np.loadtxt("xdata.csv",delimiter=',')
y = np.loadtxt("ydata.csv",delimiter=',')


for i in range(x.shape[1]):
	if ((i % 2) == 0):
		print 'plotting first'
		plt.plot(x[:,i],y[:,i])
	else:
		print 'plotting second'
		plt.scatter(x[:,i],y[:,i],color='r')

plt.show()
