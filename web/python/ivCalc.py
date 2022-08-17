import cProfile
import pstats
import json
import os
import sys
import numpy as np

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

dataInArr = convertInput(sys.argv[1])

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

from getelec_online import fit_data, plot_data

outdata = {}

def main():

    print(json.dumps(dataInArr))

    F0 = [1., 5., 20.]
    R0 = [1., 5., 50.]
    gamma0 = [1., 10., 100.]
    temp0 = [299.99999, 300., 300.01]

    xdata = 1./np.array(dataInArr[0])
    ydata = np.log(np.array(dataInArr[1]))
    W0 = np.array([1.-1e-4, 1., 1.+1e-4]) * float(dataInArr[2][0])

    fit = fit_data(xdata, ydata, F0, W0, R0, gamma0, temp0)

    popt = fit.x
    yopt = plot_data(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    yshift = max(yopt) - max(ydata)

    xth = np.linspace(min(xdata),max(xdata),32)
    yth = np.exp(plot_data(xth, popt[0], popt[1], popt[2], popt[3], popt[4]) - yshift)

    xplot = xdata / popt[0]
    xplot_th = xth / popt[0]

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": dataInArr[1], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": popt[0], "Radius": popt[2], "sigma_Aeff": 1e-9*np.exp(-yshift), \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    print(json.dumps(outdata))

with cProfile.Profile() as pr:
    main()

stats = pstats.Stats(pr)
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats()
#stats.dump_stats(filename='ivCalcProfiling.prof')
