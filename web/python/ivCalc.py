import cProfile
import pstats
import json
import os
import sys
import numpy as np

def convertInput():

    dataInArr = (sys.argv[1])

    data = json.dumps(dataInArr)
    data = data.split("}")
    data = data[0].split("{")[1]
    data = data.split("]")

    lines = []
    result = []

    for line in data:

        lined = line.split(":")[1:]

        for linedd in lined:

            linedd = linedd[1:]
            lines.append(linedd)

    for line in lines:

        if "," in line:

            _line = line.split(",")
            line = []

            for el in _line:

                if el[0] == "[":
                    el = el[1:]

                line.append(float(el))

            result.append(line)
            
        else:

            if len(line) > 0:

                if line[0] == "[":

                    line = line[1:]

                result.append([float(line)])

    return data

dataInArr = convertInput()

print(f'Converted input. Now doing calcs')

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/interfaces/web_interface'
sys.path.append(pythonpath)

print(pythonpath)

from getelec_online import fit_data

outdata = {}

def main():

    #print(json.dumps(dataInArr))

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

    #print(json.dumps(dataInArr))

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

    if float(dataInArr[3][0]) == 0:
        main()
    if float(dataInArr[3][0]) == 1:
        mainWithRadius()

#stats = pstats.Stats(pr)
#stats.sort_stats(pstats.SortKey.TIME)
#stats.print_stats()
#stats.dump_stats(filename='ivCalcProfiling.prof')
