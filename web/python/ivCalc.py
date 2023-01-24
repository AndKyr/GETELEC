import cProfile
import pstats
import json
import os
import sys
import numpy as np
import concurrent.futures

data = json.loads(sys.argv[1])

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

    xdata = 1./np.array(data['Voltage'])
    ydata = np.log(np.array(data['Current']))
    workFunction = np.full(ydata.shape, float(data['Work_function'][0]))

    xplot, xplot_th, yth, beta, radius, sigmaAeff = fit_data(xdata, ydata, workFunction, mode = "simple")

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": data['Current'], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": beta, "sigma_Aeff": sigmaAeff, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    print(json.dumps(outdata))

def main_beta():

    xdata = 1./np.array(data['Voltage'])
    ydata = np.log(np.array(data['Current']))
    workFunction = np.full(ydata.shape, float(data['Work_function'][0]))

    xplot = np.empty(ydata.shape)
    xplot_th = np.empty(ydata.shape)
    yth = np.empty(ydata.shape)
    beta = np.empty(ydata.shape)
    radius = np.empty(ydata.shape)
    sigmaAeff = np.empty(ydata.shape)

    with concurrent.futures.ThreadPoolExecutor() as executor:

        args_list = [(xdata[i], ydata[i], workFunction[i]) for i in range(len(xdata))]
        results = [executor.submit(fit_data, *args) for args in args_list]

        concurrent.futures.wait(results)

        for i, result in enumerate(results):
            xplot[i], xplot_th[i], yth[i], beta[i], radius[i], sigmaAeff[i] = result.result()

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": data['Current'], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": beta, "sigma_Aeff": sigmaAeff, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    print(json.dumps(outdata))


def mainWithRadius():

    #print(json.dumps(dataInArr))

    xdata = 1./np.array(data['Voltage'])
    ydata = np.log(np.array(data['Current']))
    workFunction = float(data['Work_function'][0])

    xplot, xplot_th, yth, beta, radius, sigmaAeff = fit_data(xdata, ydata, workFunction, mode = "withRadius")

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplot.tolist(), "yplot_mrk": data['Current'], \
                "xplot_line": xplot_th.tolist(), "yplot_line": yth.tolist(), \
                "beta": beta, "Radius": radius, "sigma_Aeff": sigmaAeff, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}

    print(json.dumps(outdata))

with cProfile.Profile() as pr:

    if data['CalculateR'] == 0:
        main()
    if data['CalculateR'] == 1:
        mainWithRadius()

#stats = pstats.Stats(pr)
#stats.sort_stats(pstats.SortKey.TIME)
#stats.print_stats()
#stats.dump_stats(filename='ivCalcProfiling.prof')
