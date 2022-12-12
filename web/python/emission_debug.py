from distutils.file_util import write_file
import json
from math import gamma
import os
import sys
import numpy as np

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/interfaces/web_interface'
sys.path.append(pythonpath)

from getelec_online import currentDensityMetal

example_Data = '{"materialType":[1],"sweepParam":[2],"field":[2,2.48,2.95,3.43,3.9,4.38,4.86,5.33,5.81,6.29,6.76,7.24,7.71,8.19,8.67,9.14,9.62,10.1,10.57,11.05,11.52,12],"radius":[50],"work_function":[4.5],"temperature":[300],"ec":[4.05],"ef":[4.61],"eg":[1.12],"gammaMetal":[10],"gammaSemi":[10],"me":[0.98],"mp":[0.5],"calculateEC":[1],"calculateNH":[0],"calculateES":[0]}'


def forceSameLength(_data):

        maxlen = max(len(_data['field']), len(_data['radius']), len(_data['work_function']), len(_data['temperature']), len(_data['ec']), len(_data['ef']), len(_data['eg']), len(_data['gammaMetal']), len(_data['gammaSemi']), len(_data['me']), len(_data['mp']))

        shortFields = [var for var in _data.keys() if (var in ['field', 'radius', 'work_function', 'temperature', 'ec', 'ef', 'eg', 'gammaMetal', 'gammaSemi', 'me', 'mp'] and len(_data[f"{var}"]) < maxlen)]

        for _field in shortFields:

            _data[f"{_field}"] = np.full(maxlen, _data[f"{_field}"])

        return _data

data = forceSameLength(json.loads(example_Data))

field = data['field']
radius = data['radius']

work_function = data['work_function']
temperature = data['temperature']

gammaMetal = data['gammaMetal']

calculateEC = data['calculateEC']

data1 = currentDensityMetal(field, radius, gammaMetal, work_function, temperature).tolist()

print(data1)