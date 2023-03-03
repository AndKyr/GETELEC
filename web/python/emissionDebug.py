from distutils.file_util import write_file
import json
from math import gamma
import os
import sys
import numpy as np

import timeit

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/interfaces/web_interface'
sys.path.append(pythonpath)

from getelec_online import current_metal_emitter, heat_metal_emitter, spectrum_metal_emitter
from getelec_online import current_semiconductor_emitter, heat_semiconductor_emitter, spectrum_semiconductor_emitter

# from getelec import Interpolator

# Interpolator.calculateAndSaveTable()

#example_Data = '{"materialType":[1],"sweepParam":[2],"field":[2,2.48,2.95,3.43,3.9,4.38,4.86,5.33,5.81,6.29,6.76,7.24,7.71,8.19,8.67,9.14,9.62,10.1,10.57,11.05,11.52,12],"radius":[50],"work_function":[4.5],"temperature":[300],"ec":[4.05],"ef":[4.61],"eg":[1.12],"gammaMetal":[10],"gammaSemi":[10],"me":[0.98],"mp":[0.5],"calculateEC":[1],"calculateNH":[1],"calculateES":[0]}'

example_Data = '{"materialType":[2],"sweepParam":[2],"field":[2,2.48,2.95,3.43,3.9,4.38,4.86,5.33,5.81,6.29,6.76,7.24,7.71,8.19,8.67,9.14,9.62,10.1,10.57,11.05,11.52,12],"radius":[50],"work_function":[4.5],"temperature":[300],"ec":[4.05],"ef":[4.61],"eg":[1.12],"gammaMetal":[10],"gammaSemi":[10],"me":[0.98],"mp":[0.49],"calculateEC":[1],"calculateNH":[1],"calculateES":[1]}'

def convertInput(toConvert):

    dataInArr = json.loads(toConvert)

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
            
            if line[0] == "[":
                line = line[1:]

            result.append([float(line)])
    
    return result

data = convertInput(example_Data)

field = data[2]
radius = data[3]

wf = data[4]
temp = data[5]

ec = data[6][0]
ef = data[7][0]

eg = data[8][0]
gammaMetal = data[9]

gammaSemi = data[10][0]
me = data[11][0]
mp = data[12][0]

calculateEC = (str(data[13][0]))[0]

data = spectrum_metal_emitter(field, radius, gammaMetal, wf, temp)

print(data.keys())

# print(f"Now measuring performance differences...")

# average_time = timeit.timeit(lambda: current_density_metal_beta(field, radius, gammaMetal, wf, temp), number = 1000)

# print("Average run time of current_density_metal_beta:", average_time / 10)

# average_time = timeit.timeit(lambda: current_density_metal(field, radius, gammaMetal, wf, temp), number = 1000)

# print("Average run time of current_density_metal:", average_time / 10)

#Average run time of current_density_metal_beta: 0.5402713211951777
#Average run time of current_density_metal: 0.18813487290171907