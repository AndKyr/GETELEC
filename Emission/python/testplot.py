import numpy as np
import matplotlib.pyplot as plt

#read data from file
x,y= np.loadtxt("out.csv",delimiter=',',unpack='True')	



#plot results
plt.figure(1)
plt.plot(x,y,'b',label="ft")  # connect points with a blue line
plt.figure(2)
plt.semilogy(x,abs(y),'b',label="ft") 
plt.show()
