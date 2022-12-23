import sys
import numpy as np
import os
from pathlib import Path
import scipy.integrate as ig

getelecRootPath = str(Path(__file__).parents[1].absolute())
sys.path.insert(0,getelecRootPath + "/src/")
import getelec as gt

import matplotlib.pyplot as plt

showFigures = True
tolerance = 1.e-3


class ConductionBandTests:

    def __init__(self) -> None:
        bar = gt.Barrier(5, 1000, 10., tabulationFolder= getelecRootPath + "/tabulated/1D_1024")
        self.emitter = gt.ConductionBandEmitter(bar)
        self.emitter.setParameters(workfunction=4.5, kT=gt.Globals.BoltzmannConstant * 1500., Ec=0.1, effectiveMass=1.)
    
    

    def effectiveMassTest(self, minMass = 0.001, maxMass = 1., Npoints = 8, plotIntegrand = False):

        masses = np.geomspace(minMass, maxMass, Npoints)
        currentDensity = np.copy(masses)
        nottinghamHeat = np.copy(masses)

        if plotIntegrand:
            fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)

        for i in range(len(masses)):
            self.emitter.setParameters(effectiveMass=masses[i])
            currentDensity[i] = self.emitter.currentDensity(mode = "slow", saveIntegrand=plotIntegrand)
            Jfast = self.emitter.currentDensity(mode = "fast")
            JfromTED = self.emitter.currentDensityFromTED()
            if (abs(1 - currentDensity[i] / JfromTED) > tolerance or abs(1 - currentDensity[i] / Jfast) > tolerance):
                print("current density from two methods not matching", currentDensity[i], " != ", JfromTED, Jfast)

            nottinghamHeat[i] = self.emitter.nottinghamHeat(mode="slow", saveIntegrand=plotIntegrand)
            Pfast = self.emitter.nottinghamHeat(mode = "fast")
            PfromTED = self.emitter.nottinghamHeatFromTED()
            if (abs(1 - nottinghamHeat[i] / PfromTED) > tolerance or abs(1 - nottinghamHeat[i] / Pfast) > tolerance):
                print("nottingham heat from two methods not matching", nottinghamHeat[i], " != ", PfromTED, Pfast)


            if plotIntegrand:
                ax1.plot(self.emitter.currentDensityIntegrandPoints, self.emitter.currentDensityIntegrandArray, \
                    label="m*=%.2g m"%masses[i])
                ax2.plot(self.emitter.nottinghamHeatIntegrandPoints, self.emitter.nottinghamHeatIntegrandArray)

        if plotIntegrand:
            ax1.set_xlabel("energy[eV]")
            ax1.set_ylabel("current density integrand [A/nm^2 / eV]")
            ax1.grid()
            ax2.set_ylabel("nottingham heat integrand [A/nm^2]")
            ax2.grid()
            # plt.legend()
            # plt.savefig("currentDensityIntegrandForMasses.png")
            if (showFigures):
                plt.show()

        plt.figure()
        plt.semilogx(masses, currentDensity)
        plt.xlabel("effective mass / electron mass")
        plt.ylabel("current density [nA / nm^2]")
        plt.grid()
        plt.savefig("currentDensity-effectiveMass.png")
        if (showFigures):
            plt.show()

    def conductionBandBottomTest(self, minEc = -5, maxEc = 1., Npoints = 4, plotSpectra = False):

        arrayEc = np.linspace(minEc, maxEc, Npoints)
        currentDensity = np.copy(arrayEc)
        for i in range(len(arrayEc)):
            self.emitter.setParameters(Ec=arrayEc[i], workfunction=4.8, kT = gt.Globals.BoltzmannConstant * 800)
            currentDensity[i] = self.emitter.currentDensity()

            if (plotSpectra):
                Energy, spectrum = self.emitter.calculateTotalEnergySpectrum()
                JfromTED = gt.Globals.SommerfeldConstant * np.trapz(spectrum, Energy)
                print("Ec = %g, J = %g,  J / JTED = %g, Npoints  = %g"%(arrayEc[i], currentDensity[i], currentDensity[i]  / JfromTED, len(Energy)))
                plt.plot(Energy, spectrum, "-", label="Ec=%.2g eV"%arrayEc[i])

        if (plotSpectra):
            plt.xlabel("energy[eV]")
            plt.ylabel("current density per energy [A/nm^2 / eV]")
            plt.grid()
            plt.legend()
            plt.savefig("spectraForDifferentEc.png")
            if (showFigures):
                plt.show()
            plt.close()

        plt.figure()
        plt.semilogy(arrayEc, currentDensity)
        plt.xlabel("Ec [eV]")
        plt.ylabel("current density [nA / nm^2]")
        plt.grid()
        plt.savefig("currentDensity-Ec.png")
        if (showFigures):
            plt.show()

        
if (__name__ == "__main__"):
    tests = ConductionBandTests()
    # tests.conductionBandBottomTest(plotSpectra=True)
    tests.effectiveMassTest(plotIntegrand=True)