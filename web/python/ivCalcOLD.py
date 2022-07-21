#! /home/bitnami/.local/bin/hug -f 
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import json
from numba import jit
import time

#get data from parameters (json format!)

dataIn = sys.argv[1]

#deal with directory stuff

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

#getelec library import

import getelec_mod as gt

#main function that is called at the end!

#input: indata_arr is an array!

def fit_fun(indata_arr): 
   
   #default arrays of float values
    F0 = [1., 5., 20.]
    R0 = [1., 5., 50.]
    gamma0 = [1., 10., 100.]
    temp0 = [299.99999, 300., 300.01]
    
    #load the data
    xdata = 1./np.array(indata_arr[0])
    ydata = np.log(np.array(indata_arr[1]))
    W0 = np.array([1.-1e-4, 1., 1.+1e-4]) * float(indata_arr[2][0])

    #magic stuff with getelec
    fit = gt.fitML(xdata,ydata, F0, W0, R0, gamma0, temp0)
    popt = fit.x
    yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    yshift = max(yopt) - max(ydata)

    xth = np.linspace(min(xdata),max(xdata),32)
    yth = np.exp(gt.MLplot(xth, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)      
    
    xplot = xdata / popt[0]

    xplot_th = xth / popt[0]

    #converting getelec function outputs to a json manually
    
    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": indata_arr[1], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": popt[0], "Radius": popt[2], "sigma_Aeff": 1e-9*np.exp(-yshift), \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    #converts dictionary to json and returns it (as string!)
    return json.dumps(outdata)


#converts string input from sys.argv and returns array of float arrays
def convertInput(input):
    _dataIn = input.split("]")
    _data = []
    for i in range(len(_dataIn)):
        _ddataIn = _dataIn[i].split("[")
        if(len(_ddataIn) > 1):
            paramListStr = _ddataIn[1].split(",")
            paramListFloat = []
            for j in range(len(paramListStr)):
                paramListFloat.append(float(paramListStr[j]))
            _data.append(paramListFloat)
        else:
            param = _ddataIn[0].split(":")[1].split("}")[0]
            _data.append([float(param)])
    return _data

#this is where data is logged as std:out, and when server sees a message starting
#with {"type":"ivCalc"... it will know it got the values from this code

#cache=true for caching...
#parallel=true for parrallelizing...
# fit_fun_jit = jit()(fit_fun)

start = time.process_time()

fit_fun(convertInput(dataIn))

elapsedTime = time.process_time() - start

print("Time for regular code: " + str(elapsedTime))

# start = time.process_time()

# fit_fun_jit(convertInput(dataIn))

# elapsedTime = time.process_time() - start

# print("Time for optimized code, 1st run: " + str(elapsedTime))

# start = time.process_time()

# fit_fun_jit(convertInput(dataIn))

# elapsedTime = time.process_time() - start

# print("Time for optimized code, 2nd run: " + str(elapsedTime))

# start = time.process_time()

# fit_fun_jit(convertInput(dataIn))

# elapsedTime = time.process_time() - start

# print("Time for optimized code, 3rd run: " + str(elapsedTime))

print(fit_fun(convertInput(dataIn)))

