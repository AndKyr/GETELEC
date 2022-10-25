#! /usr/bin/python
import numpy as np
import getelec_mod as gt

import matplotlib.pyplot as plt
import matplotlib as mb

font = 30
# mb.rcParams["font.family"] = "Serif"
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2.5
fsize = (18,10)


tab = gt.Tabulator(Nf=16, Nr=8, Ngamma=1)
tab.tabulateGamowTable(Nfield=16, Nradius=1, Ngamma=1)

print(np.shape(np.load("tabulated/GamowTable.npy")))
