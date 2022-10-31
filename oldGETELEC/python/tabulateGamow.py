#! /usr/bin/python

import getelec_mod as gt
import numpy as np


import argparse


parser = argparse.ArgumentParser(description = "Description for my parser")
parser.add_argument("-nf", "--Nfield", help = "Number of Field values", required = True, default = "256")
parser.add_argument("-nr", "--Nradius", help = "Number of Radius values", required = False, default = "1")
parser.add_argument("-ng", "--Ngamma", help = "Number of gamma values", required = False, default = "1")
parser.add_argument("-np", "--Npolynomial", help = "Number of polynomial terms for fitting G(E)", required = False, default = "4")

parser.add_argument("-nG", "--NGamow", help = "Number of Gamow points to be calculated for fitting", required = False, default = "128")
parser.add_argument("-nP", "--Nprocess", help = "Number of processes to be used for calculation", required = False, default = "1")
parser.add_argument("-o", "--output", help = "Output folder", required = False, default = "./")

args = parser.parse_args()

tab = gt.Tabulator(Nf=int(args.Nfield), Nr = int(args.Nradius), Ngamma=int(args.Ngamma), Npoly=int(args.Npolynomial), NGamow=int(args.NGamow) )
tab.tabulateGamowTable()

print("Tabulated successfully. Shape of GamowTable = ", np.shape(np.load("tabulated/GamowTable.npy")))
