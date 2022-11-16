import numpy as np

def dummy(x):

    #out = np.float(400+5*x)
    #out = np.random.rand(x)
    #out = np.double(out1)
    #print(type(x))
    #print(type(x))

    out = x/1E3

    return out

def dummy_for(x):

    #out = np.float(400+5*x)
    #out = np.random.rand(x)
    #out = np.double(out1)
    #print(type(x))
    #print(type(x))
    
    out = np.zeros(len(x))
    for i in range(len(x)):
        out[i] = x[i]/1E3

    return out
