import os
import sys
import numpy as np
import json

mainpath,filename = os.path.split(os.path.realpath(__file__))

emissionpath,mainfolder = os.path.split(mainpath)
emissionpath,mainfolder = os.path.split(emissionpath)

pythonpath = emissionpath + '/src/'
sys.path.append(pythonpath)
sys.path.append(emissionpath + '/interfaces/web_interface')

from getelec import GETELECModel as GTM

from getelec_online import current_metal_emitter, heat_metal_emitter, spectrum_metal_emitter
from getelec_online import current_semiconductor_emitter, heat_semiconductor_emitter, spectrum_semiconductor_emitter

# gtm = GTM()

# field = np.array([1, 2, 3, 4, 5])
# radius = np.full(5, 50)
# gamma = np.full(5, 10)
# workFunction = np.full(5, 4)
# temperature = np.full(5, 300)

# conductionBandBottom = 4.05
# fermiLevel = 4.61
# bandGap = 1.12
# effectiveMassConduction = 0.7
# effectiveMassValence = 0.5

# numberOfSpectrumPoints = 8

# gtm = GTM(emitterType = 'metal', field = field, radius = radius, gamma = gamma, workFunction = workFunction, temperature = temperature, numberOfSpectrumPoints = numberOfSpectrumPoints)


# gtm.run(calculateCurrent=True, calculateNottinghamHeat=True, calculateSpectrum=True)

# print(f'Current Density: ${gtm.getCurrentDensity()},\nNottingham Heat: ${gtm.getNottinghamHeat()},\nElectronSpectrum: ${gtm.getElectronSpectrum()}')

# print(gtm.field)
# gtm.setParameters(field = np.array([2, 3, 4, 5, 6]))
# print(gtm.field)

# gtm.getCurrentDensity()

# gtm.calculateCurrentDensity()
# print(gtm.getCurrentDensity())

jsonData = "{\"materialType\":1,\"sweepParam\":2,\"field\":[2,2.1,2.2,2.3,2.4,2.51,2.61,2.71,2.81,2.91,3.01,3.11,3.21,3.31,3.41,3.52,3.62,3.72,3.82,3.92,4.02,4.12,4.22,4.32,4.42,4.53,4.63,4.73,4.83,4.93,5.03,5.13,5.23,5.33,5.43,5.54,5.64,5.74,5.84,5.94,6.04,6.14,6.24,6.34,6.44,6.55,6.65,6.75,6.85,6.95,7.05,7.15,7.25,7.35,7.45,7.56,7.66,7.76,7.86,7.96,8.06,8.16,8.26,8.36,8.46,8.57,8.67,8.77,8.87,8.97,9.07,9.17,9.27,9.37,9.47,9.58,9.68,9.78,9.88,9.98,10.08,10.18,10.28,10.38,10.48,10.59,10.69,10.79,10.89,10.99,11.09,11.19,11.29,11.39,11.49,11.6,11.7,11.8,11.9,12],\"radius\":[50],\"work_function\":[4.5],\"temperature\":[300],\"ec\":[4.05],\"ef\":[4.61],\"eg\":[1.12],\"gammaMetal\":[10],\"gammaSemi\":[10],\"me\":[0.98],\"mp\":[0.5],\"calculateEC\":1,\"calculateNH\":1,\"calculateES\":0}"


def forceSameLength(data):

    max_len = max([len(data[field]) for field in ['field', 'radius', 'work_function', 'temperature', 'gammaMetal', 'gammaSemi']])

    for field in ['field', 'radius', 'work_function', 'temperature', 'gammaMetal', 'gammaSemi']:

        while len(data[field]) < max_len:
            
            data[field].append(data[field][-1])

    
    return data

def main():

    data = forceSameLength(json.loads(jsonData))

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
            
            data1 = current_metal_emitter(data['field'], data['radius'], data['gammaMetal'], data['work_function'], data['temperature'])
        
        if data['calculateNH'] == 1:

            data2 = heat_metal_emitter(data['field'], data['radius'], data['gammaMetal'], data['work_function'], data['temperature'])

        if data['calculateES'] == 1:

            res = spectrum_metal_emitter(data['field'], data['radius'], data['gammaMetal'], data['work_function'], data['temperature'])

            for i, arr in enumerate(res['energy']):
                res['energy'][i] = arr.tolist()

            for i, arr in enumerate(res['electronCount']):
                res['electronCount'][i] = arr.tolist()

            data3 = res['energy']
            data3e = res['electronCount']


            
    elif data['materialType'] == 2:

        if data['calculateEC'] == 1:

            data4 = current_semiconductor_emitter(data['field'], data['radius'], data['gammaSemi'], data['ec'], data['ef'], data['eg'], data['temperature'], data['me'], data['mp']).tolist()

        if data['calculateNH'] == 1:  
              
            data5 = heat_semiconductor_emitter(data['field'], data['radius'], data['gammaSemi'], data['ec'], data['ef'], data['eg'], data['temperature'], data['me'], data['mp']).tolist()
  
        if data['calculateES'] == "1":

            res = spectrum_semiconductor_emitter(data['field'], data['radius'], data['gammaSemi'], data['ec'], data['ef'], data['eg'], data['temperature'], data['me'], data['mp'])

            for i, arr in enumerate(res['energy']):
                res['energy'][i] = arr.tolist()

            for i, arr in enumerate(res['electronCount']):
                res['electronCount'][i] = arr.tolist()

            data6 = res['energy']
            data6e = res['electronCount']

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
