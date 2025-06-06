#! /usr/bin/python

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib



mainpath,filename = os.path.split(os.path.realpath(__file__))
emissionpath,mainfolder = os.path.split(mainpath)
pythonpath = emissionpath + '/python'
sys.path.append(pythonpath)

import getelec_mod as getelec_old

x = np.linspace(0.,3.,32)
V = 5*x - 0.1 * x**2

this = getelec_old.emission_create(xr = x, Vr = V, mode = -10)
this.print_C_data()
this.print_data(True)
