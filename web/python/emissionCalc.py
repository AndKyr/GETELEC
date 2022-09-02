from distutils.file_util import write_file
import json
from math import gamma
import os
import sys
import numpy as np

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
            result.append([float(line)])
    
    return result

def main():

    data = convertInput()
    
    materialType = (str(data[0][0]))[0]
    sweepParam = str(data[1][0])[0]
    field = data[2]
    radius = data[3]
    wf = data[4]
    temp = data[5]
    ec = data[6][0]
    ef = data[7][0]
    eg = data[8][0]
    gammaMetal = data[9]
    gammaSemi = data[10]
    me = data[11][0]
    mp = data[12][0]
    calculateEC = (str(data[13][0]))[0]
    calculateES = (str(data[15][0]))[0]
    calculateNH = (str(data[14][0]))[0]

    data1 = []
    data2 = []
    data3 = []
    data4 = []
    data5 = []
    data6 = []

    if sweepParam == "2":
        sweepParam = "field"
    elif sweepParam == "3":
        sweepParam = "radius"
    elif sweepParam == "4":
        sweepParam = "wf"
    elif sweepParam == "5":
        sweepParam = "temp"

    if materialType == "1":

        if calculateEC == "1":
            data1 = current_metal_emitter(field, radius, gammaMetal, wf, temp).tolist()
        
        if calculateNH == "1":
            data2 = heat_metal_emitter(field, radius, gammaMetal, wf, temp).tolist()

        if calculateES == "1":
            # data3 = spectrum_metal_emitter(field, radius, gammaMetal, wf, temp)
            data3

    
    elif materialType == "2":

        if calculateEC == "1":
            data4 = current_semiconductor_emitter(field, radius, gammaSemi, ec, ef, eg, temp, me, mp).toList()

        if calculateNH == "1":    
            data5 = heat_semiconductor_emitter(field, radius, gammaSemi, ec, ef, eg, temp, me, mp).tolist()
  
        if calculateES == "1":
            # data6 = spectrum_semiconductor_emitter(field, radius, gammaSemi, ec, ef, eg, temp, me, mp)
            data6

    # print(json.dumps(dataInArr))
    

    outdata = {"materialType": materialType, "sweepParam": sweepParam,
        "field": field, "radius": radius, "work_function": wf, "temperature": temp,
        "ec": ec, "ef": ef, "eg": eg, "gammaMetal": gammaMetal, "gammaSemi": gammaSemi,
        "me": me, "mp": mp, "metalEC": data1, "metalNH": data2, "metalES": data3,
        "semiEC": data4, "semiNH": data5, "semiES": data6
    }


    # asdas
    print(json.dumps(outdata))

main()