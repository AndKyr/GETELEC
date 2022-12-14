import sys
import numpy as np
import os
from pathlib import Path

getelecRootPath = str(Path(__file__).parents[2].absolute())
sys.path.insert(0,getelecRootPath + "/src/")
import getelec as gt

import matplotlib.pyplot as plt



class ConductionBandTests:



    def __init__(self) -> None:
        bar = gt.Barrier(5, 1000, 10., tabulationFolder="tabulated/1D_1024")
        sup = gt.Supply()
        self.emitter = gt.ConductionBandEmitter(bar, sup)
        self.emitter.setParameters(workfunction=4., kT=gt.Globals.BoltzmannConstant * 800., Ec=0.1, effectiveMass=1.01)
    
    
    def effectiveMassTest(self, minMass = 0.01, maxMass = 10., Npoints =64, plotSpectra = False):

        masses = np.geomspace(minMass, maxMass, Npoints)
        currentDensity = np.copy(masses)
        for i in range(len(masses)):
            self.emitter.setParameters(effectiveMass=masses[i])
            currentDensity[i] = self.emitter.currentDensityFast()

            print("m* = %g, J = %g"%(masses[i], currentDensity[i]))

            if (plotSpectra):
                Energy, spectrum = self.emitter.normalEnergyDistribution()
                plt.plot(Energy, spectrum, label="m*=%.2g m"%masses[i])

        if (plotSpectra):
            plt.xlabel("energy[eV]")
            plt.ylabel("current density per energy [A/nm^2 / eV]")
            plt.grid()
            plt.legend()
            plt.savefig("spectraForDifferentEffectiveMasses.png")

        plt.figure()
        plt.semilogx(masses, currentDensity)
        plt.xlabel("effective mass / electron mass")
        plt.ylabel("current density [nA / nm^2]")
        plt.grid()
        plt.savefig("currentDensity-effectiveMass.png")

    def conductionBandBottomTest(self, minEc = -5, maxEc = 1., Npoints = 128, plotSpectra = False):

        arrayEc = np.linspace(minEc, maxEc, Npoints)
        currentDensity = np.copy(arrayEc)
        for i in range(len(arrayEc)):
            self.emitter.setParameters(Ec=arrayEc[i])
            currentDensity[i] = self.emitter.currentDensityFast()

            # print("Ec = %g, J = %g, lowEnergyLimit = %g"%(arrayEc[i], currentDensity[i], self.emitter._lowEnergyLimit))

            if (plotSpectra):
                Energy, spectrum = self.emitter.normalEnergyDistribution()
                plt.plot(Energy, spectrum, label="Ec=%.2g eV"%arrayEc[i])

        if (plotSpectra):
            plt.xlabel("energy[eV]")
            plt.ylabel("current density per energy [A/nm^2 / eV]")
            plt.grid()
            plt.legend()
            plt.savefig("spectraForDifferentEc.png")

        plt.figure()
        plt.semilogy(arrayEc, currentDensity)
        plt.xlabel("Ec [eV]")
        plt.ylabel("current density [nA / nm^2]")
        plt.grid()
        plt.savefig("currentDensity-Ec.png")
        
if (__name__ == "__main__"):
    tests = ConductionBandTests()

    # tests.effectiveMassTest()
    tests.conductionBandBottomTest()