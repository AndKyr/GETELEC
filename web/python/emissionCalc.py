import json
import os
import sys

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

from getelec_online import current_metal_emitter, heat_metal_emitter, spectrum_metal_emitter
from getelec_online import current_semiconductor_emitter, heat_semiconductor_emitter, spectrum_semiconductor_emitter

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
                line.append(float(el))
            result.append(line)
            
        else:
            result.append(float(line))
    
    return result


def main():

    F0 = [1., 5., 20.]
    R0 = [15., 25., 50.]
    gamma0 = [1., 10., 100.]
    temp0 = [299.99999, 300., 300.01]

    data = convertInput()
    
    materialType = data[0]
    sweepParam = data[1]
    field = data[2]
    radius = data[3]
    wf = data[4]
    temp = data[5]
    ec = data[6]
    ef = data[7]
    eg = data[8]
    gammaMetal = data[9]
    gammaSemi = data[10]
    me = data[11]
    mp = data[11]
    
    # for i in range(len(lines)):
    #     if "," in lines[i]:
    #         _values = lines[i].split(",")
    #         for num in _values:
    #             result[i].append(float(num))

    #     else:
    #         result.append(float(lines[i]))


    # print(json.dumps(dataInArr))
    

    # outdata1 = {"type": "current"}

    # outdata2 = {"type": "nottinghamHeat"}

    # outdata3 = {"type": "electronSpectrum"}


    # print(json.dumps(outdata1))
    # print(json.dumps(outdata2))
    # print(json.dumps(outdata3))


main()