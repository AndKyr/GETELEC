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

from getelec_online import fit_data

outdata = {}

def main():

    print(json.dumps(dataInArr))

    xdata = 1./np.array(dataInArr[0])
    ydata = np.log(np.array(dataInArr[1]))
    workFunction = float(dataInArr[2][0])

    xplot, xplot_th, yth, beta, radius, sigmaAeff = fit_data(xdata, ydata, workFunction, mode = "simple")

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": dataInArr[1], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": beta, "sigma_Aeff": sigmaAeff, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    print(json.dumps(outdata))

def mainWithRadius():
    print(json.dumps(dataInArr))

    xdata = 1./np.array(dataInArr[0])
    ydata = np.log(np.array(dataInArr[1]))
    workFunction = float(dataInArr[2][0])

    xplot, xplot_th, yth, beta, radius, sigmaAeff = fit_data(xdata, ydata, workFunction, mode = "withRadius")

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": dataInArr[1], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": beta, "Radius": radius, "sigma_Aeff": sigmaAeff, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    print(json.dumps(outdata))

with cProfile.Profile() as pr:
    main()

#stats = pstats.Stats(pr)
#stats.sort_stats(pstats.SortKey.TIME)
#stats.print_stats()
#stats.dump_stats(filename='ivCalcProfiling.prof')
