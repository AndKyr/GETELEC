import json
import os
import sys


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

from getelec_online import current_metal_emitter, heat_metal_emitter, spectrum_metal_emitter
from getelec_online import current_semiconductor_emitter, heat_semiconductor_emitter, spectrum_semiconductor_emitter

def main():

    F0 = [1., 5., 20.]
    R0 = [15., 25., 50.]
    gamma0 = [1., 10., 100.]
    temp0 = [299.99999, 300., 300.01]

    

    outdata1 = {"type": "current"}

    outdata2 = {"type": "nottinghamHeat"}

    outdata3 = {"type": "electronSpectrum"}


    print(json.dumps(outdata1))
    print(json.dumps(outdata2))
    print(json.dumps(outdata3))


main()