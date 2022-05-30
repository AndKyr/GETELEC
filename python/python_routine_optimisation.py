#from tkinter import TRUE
from unittest import result
import numpy as np
from multiprocessing import Pool
import multiprocessing as mp
import time
import getelec_tabulator as gtab
import pandas as pd
import threading as th

def current(f,r,g,w,t):

    kBoltz = 8.6173324e-5 

    tab = gtab.Tabulator()

    kT = kBoltz * t

    metal_emitter = gtab.Metal_Emitter(tab)
    #j_metal = np.zeros(len(f))

    #for i in range(len(f)):
    metal_emitter.emitter.Define_Barrier_Parameters(f, r, g)
    metal_emitter.emitter.Interpolate_Gammow()

    metal_emitter.Define_Emitter_Parameters(w, kT)

    j_metal = metal_emitter.Current_Density()

    return j_metal

def current2(f,r,g,w,t):

    kBoltz = 8.6173324e-5 

    tab = gtab.Tabulator()

    kT = kBoltz * t

    metal_emitter = gtab.Metal_Emitter(tab)
    j_metal = np.zeros(len(f))

    for i in range(len(f)):
        metal_emitter.emitter.Define_Barrier_Parameters(f[i], r[i], g[i])
        metal_emitter.emitter.Interpolate_Gammow()

        metal_emitter.Define_Emitter_Parameters(w[i], kT[i])

        j_metal[i] = metal_emitter.Current_Density()

    return j_metal

def dummy(x,y):
    out = x*y/1E3
    print(out)
    return out

def dummy_for(x,y):
    out = np.zeros(len(x))
    for i in range(len(x)):
        out[i] = x[i]*y[i]/1E3
    return out

dimension = 1E2

field = np.ones(int(dimension))*10
radius = np.ones(int(dimension))*10
gamma = np.ones(int(dimension))*10
work = np.ones(int(dimension))*4.5
temp = np.ones(int(dimension))*300

start1 = time.time()
r = current2(field, radius, gamma, work, temp)
#g = current(field[1], radius[1],gamma[1],work[1],temp[1])
end1 = time.time()
#print("nunp",end1-start1)

start = time.time()
pool = Pool()
args=[(field, radius, gamma, work, temp)]
t = mp.Process(current2,args)
t = t.get()
end = time.time()
pool.close()
#print("pool",end-start)

start2 = time.time()
p1 = mp.Process(target=current2, args=(field, radius, gamma, work, temp))
p1.start()
#r = p1.join()
#print(r)
end2 = time.time()

start3 = time.time()
p2 = th.Thread(target=current2, args=(field, radius, gamma, work, temp))
p2.start()
r = p2.join()
#print(r.get())
end3 = time.time()

print("nunp",end1-start1)
print("nunp",end2-start2)
print("nunp",end3-start3)
print("pool",end-start)

"""I will get to include graphs about relative errors for each of the optimisatio steps"""
#print(df_out)

