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

from getelec_online import current_density_metal, heat_metal_emitter, spectrum_metal_emitter
from getelec_online import current_semiconductor_emitter, heat_semiconductor_emitter, spectrum_semiconductor_emitter

from getelec_online import current_density_metal_beta


def main():

    data = json.loads(sys.argv[1])

    data1 = []
    data2 = []
    data3 = []
    data4 = []
    data5 = []
    data6 = []
    data3e = []
    data6e = []

    # convert all relevant fields to np arrays

    for field in ['field', 'radius', 'work_function', 'temperature', 'ec', 'ef', 'eg', 'gammaMetal', 'gammaSemi', 'me', 'mp']:

        data[field] = np.array(data[field])
        
    if data['sweepParam'] == 2:

        data['sweepParam'] = "field"

    elif data['sweepParam'] == 3:

        data['sweepParam'] = "radius"

    elif data['sweepParam'] == 4:

        data['sweepParam'] = "work_function"

    elif data['sweepParam'] == 5:

        data['sweepParam'] = "temperature"

    if data['materialType'] == 1:

        if data['calculateEC'] == 1:
            
            data1 = current_density_metal(data['field'], data['radius'], data['gammaMetal'], data['work_function'], data['temperature']).tolist()
        
        if data['calculateNH'] == 1:

            data2 = heat_metal_emitter(data['field'], data['radius'], data['gammaMetal'], data['work_function'], data['temperature']).tolist()

        if data['calculateES'] == 1:

            energies, electronCounts = spectrum_metal_emitter(data['field'], data['radius'], data['gammaMetal'], data['work_function'], data['temperature'])

            data3 = []
            data3e = []

            for energy in energies:
                data3.append(energy.tolist())

            for electronCount in electronCounts:
                data3e.append(electronCount.tolist())
            
    elif data['materialType'] == 2:

        if data['calculateEC'] == 1:

            data4 = current_semiconductor_emitter(data['field'], data['radius'], data['gammaSemi'], data['ec'], data['ef'], data['eg'], data['temperature'], data['me'], data['mp']).tolist()

        if data['calculateNH'] == 1:  
              
            data5 = heat_semiconductor_emitter(data['field'], data['radius'], data['gammaSemi'], data['ec'], data['ef'], data['eg'], data['temperature'], data['me'], data['mp']).tolist()
  
        if data['calculateES'] == "1":

            energy, electronCount = spectrum_semiconductor_emitter(data['field'], data['radius'], data['gammaSemi'], data['ec'], data['ef'], data['eg'], data['temperature'], data['me'], data['mp']).tolist()
            data6 = energy.tolist()
            data6e = electronCount.tolist()

    outdata = {

        "materialType": data['materialType'], 
        "sweepParam": data['sweepParam'], 
        "field": data['field'].tolist(), 
        "radius": data['radius'].tolist(), 
        "work_function": data['work_function'].tolist(), 
        "temperature": data['temperature'].tolist(), 
        "ec": data['ec'].tolist(), 
        "ef": data['ef'].tolist(), 
        "eg": data['eg'].tolist(), 
        "gammaMetal": data['gammaMetal'].tolist(), 
        "gammaSemi": data['gammaSemi'].tolist(), 
        "me": data['me'].tolist(), 
        "mp": data['mp'].tolist(), 
        "metalEC": data1, 
        "metalNH": data2, 
        "metalESenergy": data3, 
        "metalESelcount": data3e, 
        "semiEC": data4, 
        "semiNH": data5, 
        "semiESenergy": data6, 
        "semiESelcount": data6e

    }

    print(json.dumps(outdata))


main()