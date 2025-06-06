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

    voltageData = np.array(data['Voltage'])
    currentData = np.array(data['Current'])
    workFunction = np.full(currentData.shape, float(data['Work_function'][0]))
    fitRadius = bool(data["CalculateR"])

    ivFitter = IVDataFitter()

    ivFitter.setIVdata(voltageData=voltageData, currentData=currentData)
    ivFitter.setParameterRange(workFunction=workFunction)
    if fitRadius:
        ivFitter.setParameterRange(radius=[1., 10. , 2000.])
        
    ivFitter.fitIVCurve()

    voltageCurveData = 1/np.linspace(1/max(voltageData), 1/min(voltageData), 128)

    fittedCurrent = ivFitter.getOptCurrentCurve(voltageCurveData)

    xplotMrk = 1 / (ivFitter.fittingParameters["fieldConversionFactor"] * voltageData)

    xplotLine = 1 / (ivFitter.fittingParameters["fieldConversionFactor"] * voltageCurveData)

    ivFitter.orthodoxyStatus()




    # Andreas please modify fields below?

    outdata = { "type": "ivCalc",
                "xplot_mrk": xplotMrk.tolist(), "yplot_mrk": data['Current'], \
                "xplot_line": xplotLine.tolist(), "yplot_line": fittedCurrent.tolist(), \
                "beta": ivFitter.fittingParameters["fieldConversionFactor"], "sigma_Aeff": ivFitter.preFactor, \
                "xAxisUnit": "1 / (Local Field [V/nm])", "yAxisUnit": "Current [Amps]", "orthodoxyMessage": ivFitter.orthodoxyMessage}
    
    if (fitRadius):
        outdata["Radius"] = ivFitter.fittingParameters["radius"]
    
    print(json.dumps(outdata))

if __name__ == "__main__":
    main()

