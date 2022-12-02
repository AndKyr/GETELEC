import json
import numpy as np

def forceSameLength(_data):

        maxlen = max(len(_data['field']), len(_data['radius']), len(_data['wf']), len(_data['temp']), len(_data['ec']), len(_data['ef']), len(_data['eg']), len(_data['gammaMetal']), len(_data['gammaSemi']), len(_data['me']), len(_data['mp']))

        shortFields = [var for var in _data.keys() if (var in ['field', 'radius', 'wf', 'temp', 'ec', 'ef', 'eg', 'gammaMetal', 'gammaSemi', 'me', 'mp'] and len(_data[f"{var}"]) < maxlen)]

        for _field in shortFields:
            
            _data[f"{_field}"] = np.concatenate(_data[f"{_field}"], np.ones((maxlen - len(_data[f"{_field}"])) * _data[f'{_field}'][-1]))

        return _data


def main():

    res = json.loads('{"materialType":1,"sweepParam":2,"field":[2,2.48,2.95,3.43,3.9,4.38,4.86,5.33,5.81,6.29,6.76,7.24,7.71,8.19,8.67,9.14,9.62,10.1,10.57,11.05,11.52,12],"radius":[50],"work_function":[4.5],"temperature":[300],"ec":[4.05],"ef":[4.61],"eg":[1.12],"gammaMetal":[10],"gammaSemi":[10],"me":[0.98],"mp":[0.5],"calculateEC":1,"calculateNH":0,"calculateES":0}')

    res['field'] = np.array(res['field'])
    res['radius'] = np.array(res['radius'])
    res['wf'] = np.array(res['work_function'])
    res['temp'] = np.array(res['temperature'])

    res['ec'] = np.array(res['ec'])
    res['ef'] = np.array(res['ef'])
    res['eg'] = np.array(res['eg'])

    res['gammaMetal'] = np.array(res['gammaMetal'])
    res['gammaSemi'] = np.array(res['gammaSemi'])

    res['me'] = np.array(res['me'])
    res['mp'] = np.array(res['mp'])

    res = forceSameLength(res)

    print(res)

if __name__ == "__main__":
    main()

