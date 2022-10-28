#! /usr/bin/python

import getelec_mod as gt
import sys
import numpy as np


tab = gt.Tabulator(Nf=16, Nr=8, Ngamma=1)
tab.tabulateGamowTable(Nfield=int(sys.argv[1]), Nradius = int(sys.argv[2]), Ngamma=int(sys.argv[3]))

print("Tabulated successfully. Shape of GamowTable = ", np.shape(np.load("tabulated/GamowTable.npy")))
