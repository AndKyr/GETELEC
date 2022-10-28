import numpy as np
from scipy.interpolate import interpn
def value_func_3d(x, y):
    return np.array([2 * x + 3 * y, x + y])

#x = np.logspace(0, 1, 5)
x = np.linspace(1., 10., 15)

y = np.linspace(0, 5, 16)
values = np.ones([15,16,2])


for i in range(len(x)):
    for j in range(len(y)):
        values[i,j,:] = value_func_3d(x[i], y[j])

#Evaluate the interpolating function at a point

point = np.array([2.21, 3.12])
print(interpn((x,y), values, point)[0])
print( value_func_3d(*point))