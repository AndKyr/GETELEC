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

from getelec_online import current_metal_emitter, heat_metal_emitter, spectrum_metal_emitter
from getelec_online import current_semiconductor_emitter, heat_semiconductor_emitter, spectrum_semiconductor_emitter

def main():

    data = json.loads(sys.argv[1])

    print(data)
    
    materialType = data['materialType']
    sweepParam = data['sweepParam']
    field = data['field']
    radius = data['radius']
    wf = data['work_function']
    temp = data['temperature']

    ec = data['ec']
    ef = data['ef']
    eg = data['eg']

    gammaMetal = data['gammaMetal']
    gammaSemi = data['gammaSemi']

    me = data['me']
    mp = data['mp']

    calculateEC = data['calculateEC']
    calculateES = data['calculateES']
    calculateNH = data['calculateNH']

    data1 = []
    data2 = []
    data3 = []
    data4 = []
    data5 = []
    data6 = []
    data3e = []
    data6e = []

    if sweepParam == 2:

        sweepParam = "field"

    elif sweepParam == 3:

        sweepParam = "radius"

    elif sweepParam == 4:

        sweepParam = "wf"

    elif sweepParam == 5:

        sweepParam = "temp"

    if materialType == 1:

        if calculateEC == 1:

            print(temp)
            print(field)
            
            #data1 = current_metal_emitter(field, radius, gammaMetal, wf, temp).tolist()
        
        if calculateNH == 1:

            data2 = heat_metal_emitter(field, radius, gammaMetal, wf, temp).tolist()

        if calculateES == 1:

            energies, electronCounts = spectrum_metal_emitter(field, radius, gammaMetal, wf, temp)

            data3 = []
            data3e = []

            for energy in energies:
                data3.append(energy.tolist())

            for electronCount in electronCounts:
                data3e.append(electronCount.tolist())
            
    elif materialType == 2:

        if calculateEC == 1:

            data4 = current_semiconductor_emitter(field, radius, gammaSemi, ec, ef, eg, temp, me, mp).tolist()

        if calculateNH == 1:  
              
            data5 = heat_semiconductor_emitter(field, radius, gammaSemi, ec, ef, eg, temp, me, mp).tolist()
  
        if calculateES == "1":

            energy, electronCount = spectrum_semiconductor_emitter(field, radius, gammaSemi, ec, ef, eg, temp, me, mp).tolist()
            data6 = energy.tolist()
            data6e = electronCount.tolist()

    outdata = {"materialType": materialType, "sweepParam": sweepParam,
        "field": field, "radius": radius, "work_function": wf, "temperature": temp,
        "ec": ec, "ef": ef, "eg": eg, "gammaMetal": gammaMetal, "gammaSemi": gammaSemi,
        "me": me, "mp": mp, "metalEC": data1, "metalNH": data2, "metalESenergy": data3,
        "metalESelcount": data3e, "semiEC": data4, "semiNH": data5, "semiESenergy": data6,
        "semiESelcount": data6e
    }

    print(json.dumps(outdata))

main()