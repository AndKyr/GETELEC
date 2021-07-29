#! /home/bitnami/.local/bin/hug -f 
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import json
import hug

mainpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,mainfolder = os.path.split(mainpath)
pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

import getelec_mod as gt

@hug.post()
def fit_fun(body):
    print("the input is", body, "of type", type(body))
    
    
    F0 = [1., 5., 14.]
    R0 = [1., 5., 50.]
    gamma0 = [1., 10., 100.]
    Temp0 = [299.99999, 300., 300.01]
    
    #load the data
    indata_dict = body#json.loads(body)
    xdata = 1./np.array(indata_dict['Voltage'])
    ydata = np.log(np.array(indata_dict['Current']))
    W0 = np.array([1.-1e-4, 1., 1.+1e-4]) * float(indata_dict['Work_function'])

    fit= gt.fitML(xdata,ydata, F0, W0, R0, gamma0, Temp0)
    popt = fit.x
    yopt = gt.MLplot(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    yshift = max(yopt) - max(ydata)

    xth = np.linspace(min(xdata),max(xdata),32)
    yth = np.exp(gt.MLplot(xth, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)      
    
    xplot = xdata / popt[0]

    xplot_th = xth / popt[0]
    
    outdata = {'xplot_mrk': xplot.tolist(), 'yplot_mrk': indata_dict['Current'], \
                'xplot_line': xplot_th.tolist(), 'yplot_line': yth.tolist(), \
                'beta': popt[0], 'Radius': popt[2], 'sigma_Aeff': 1e-9*np.exp(-yshift), \
                'xAxisUnit': "1 / (Local Field [V/nm])", 'yAxisUnit': "Current [Amps]"}
    
    return json.dumps(outdata)


data  = {"Voltage": [2.413e+02, 2.511e+02, 2.622e+02, 2.706e+02, 2.803e+02, 2.915e+02, 2.999e+02, 3.096e+02, 3.208e+02, 3.305e+02, 3.403e+02, 3.515e+02, 3.612e+02, 3.710e+02, 3.808e+02, 3.891e+02, 4.003e+02, 4.100e+02, 4.184e+02, 4.338e+02, 4.491e+02, 4.644e+02, 4.826e+02, 4.993e+02], \
          "Current": [8.719e-01, 1.582e+00, 2.967e+00, 5.038e+00, 8.555e+00, 1.406e+01, 2.309e+01, 3.670e+01, 5.643e+01, 8.678e+01, 1.249e+02, 1.797e+02, 2.502e+02, 3.371e+02, 4.540e+02, 6.116e+02, 7.459e+02, 9.720e+02, 1.267e+03, 1.764e+03, 2.376e+03, 3.096e+03, 4.310e+03, 5.617e+03], "Work_function": 4.5
        }


outdata = fit_fun(data)



print(outdata)


