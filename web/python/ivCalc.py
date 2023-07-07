import json
import os
import sys

import numpy as np

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/src/'

sys.path.append(pythonpath)

data = json.loads(sys.argv[1])

from IVDataFitter import IVDataFitter

def main():

    voltageData = 1./np.array(data['Voltage'])
    currentData = np.log(np.array(data['Current']))
    workFunction = np.full(currentData.shape, float(data['Work_function'][0]))

    ivFitter = IVDataFitter()

    ivFitter.setIVdata(voltageData=voltageData, currentData=currentData)
    ivFitter.setParameterRange()
    ivFitter.fitIVCurve()

    fittedCurrent = ivFitter.getOptCurrentCurve(voltageData)

    # Andreas please modify fields below?

    outdata = { "type": "ivCalc",
                "xplot_mrk": fittedCurrent.tolist(), "yplot_mrk": data['Current'], \
                "xplot_line": fittedCurrent.tolist(), "yplot_line": data['Current'], \
                "beta": -1, "sigma_Aeff": ivFitter.preFactor, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]"}
    
    print(json.dumps(outdata))

if __name__ == "__main__":
    main()

