import sys
import numpy as np
from pathlib import Path
import scipy.integrate as ig

getelecRootPath = str(Path(__file__).parents[1].absolute())
sys.path.insert(0,getelecRootPath + "/src/")
import getelec as gt

import matplotlib.pyplot as plt
import matplotlib as mb

font = 15
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
figureSize = [16,10]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

showFigures = True
tolerance = 1.e-2


class BandEmitterTests:

    def __init__(self) -> None:
        bar = gt.Barrier(5, 1000, 10., tabulationFolder= getelecRootPath + "/tabulated/1D_1024")
        self.conductionEmitter = gt.ConductionBandEmitter(bar)
        self.valenceEmitter = gt.ValenceBandEmitter(bar)
        self.conductionEmitter.setParameters(workfunction=4.5, kT=gt.Globals.BoltzmannConstant * 1500., energyBandLimit=0.1, effectiveMass=1.)
        self.valenceEmitter.setParameters(workfunction=4.5, kT=gt.Globals.BoltzmannConstant * 1500., energyBandLimit = -.3, effectiveMass=.1)
    
    

    def effectiveMassTest(self, minMass = 0.001, maxMass = 1., Npoints = 4, plotIntegrand = False, plotSpectra = False):

        masses = np.geomspace(minMass, maxMass, Npoints)
        currentDensityConduction = np.copy(masses)
        nottinghamHeatConduction = np.copy(masses)
        currentDensityValence = np.copy(masses)
        nottinghamHeatValence = np.copy(masses)

        if plotIntegrand:
            fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=tuple(figureSize))
            fig2, (ax3, ax4) = plt.subplots(2,1, sharex=True, figsize=tuple(figureSize))
        else:
            ax1 = None
            ax2 = None
            ax3 = None
            ax4 = None

        for i in range(len(masses)):
            self.conductionEmitter.setParameters(effectiveMass=masses[i])
            self.valenceEmitter.changeParameters(effectiveMass=masses[i])
            currentDensityConduction[i], nottinghamHeatConduction[i] = self.conductionEmitter.runFullTest(ax1, ax2, label=r"$, m=%.2g$"%masses[i])
            currentDensityValence[i], nottinghamHeatValence[i] = self.valenceEmitter.runFullTest(ax3, ax4, label=r"$, m=%.2g$"%masses[i])


        if plotIntegrand:
            plt.figure(fig)
            ax4.set_xlabel(r"E - E_F [eV]")
            ax3.set_ylabel(r"integrand [A / nm$^2$ eV]")
            ax3.grid()
            ax4.set_ylabel(r"P_N integrand  [A/nm$^2$]")
            ax2.grid()
            ax1.legend()
            ax2.legend()
            plt.savefig("conductionCurrentDensityIntegrandForMasses.png")

            plt.figure(fig2)
            ax4.set_xlabel(r"E - E_F [eV]")
            ax3.set_ylabel(r"integrand [A / nm$^2$ eV]")
            ax3.grid()
            ax4.set_ylabel(r"P_N integrand  [A/nm$^2$]")
            ax4.grid()
            ax3.legend()
            ax4.legend()

            plt.savefig("valenceCurrentDensityIntegrandForMasses.png")
            if (showFigures):
                plt.show()

        plt.figure()
        plt.semilogx(masses, currentDensityConduction, '.-')
        plt.xlabel(r"$m^* / m_e$")
        plt.ylabel(r"current density [nA / nm$^2$]")
        plt.grid()
        plt.savefig("currentDensity-effectiveMass.png")
        if (showFigures):
            plt.show()

    def conductionBandBottomTest(self, minEc = -5, maxEc = 1., Npoints = 4, plotSpectra = False):

        arrayEc = np.linspace(minEc, maxEc, Npoints)
        currentDensity = np.copy(arrayEc)
        for i in range(len(arrayEc)):
            self.conductionEmitter.setParameters(energyBandLimit=arrayEc[i], workfunction=4.8, kT = gt.Globals.BoltzmannConstant * 800)
            currentDensity[i] = self.conductionEmitter.currentDensity()

            if (plotSpectra):
                self.conductionEmitter.calculateTotalEnergySpectrum()
                Energy, spectrum = self.conductionEmitter.totalEnergySpectrumArrays()
                JfromTED = gt.Globals.SommerfeldConstant * np.trapz(spectrum, Energy)
                print("Ec = %g, J = %g,  J / JTED = %g, Npoints  = %g"%(arrayEc[i], currentDensity[i], currentDensity[i]  / JfromTED, len(Energy)))
                plt.plot(Energy, spectrum, "-", label="Ec=%.2g eV"%arrayEc[i])

        if (plotSpectra):
            plt.xlabel("energy[eV]")
            plt.ylabel(r"current density per energy [A/nm$^2$ / eV]")
            plt.grid()
            plt.legend()
            plt.savefig("spectraForDifferentEc.png")
            if (showFigures):
                plt.show()
            plt.close()

        plt.figure()
        plt.semilogy(arrayEc, currentDensity)
        plt.xlabel("Ec [eV]")
        plt.ylabel(r"current density [nA / nm$^2$]")
        plt.grid()
        plt.savefig("currentDensity-Ec.png")
        if (showFigures):
            plt.show()


class tabulatorTests:

    def __init__(self) -> None:
        pass

    def compareOldAndNew(self):
        model = gt.GETELECModel(emitterType="metal")
        model.setTabulationPath("tabulated/1D_512_new")
        model.emitter.barrier.calculateAndSaveTable(NField=64, NRadius=1, dataFolder="tabulated/1D_512_new")
        # model.emitter.barrier._loadTablesFromFileIfPossible()
        model.calculateCurrentDensity()
        print("from new = ", model.getCurrentDensity())


        gt._setTabulationPath("tabulated/1D_1024")
        model = gt.GETELECModel(emitterType="metal")
        model.calculateCurrentDensity()
        print("from old = ", model.getCurrentDensity())


if (__name__ == "__main__"):
    # tests = BandEmitterTests()
    # tests.conductionBandBottomTest(plotSpectra=True)
    # tests.effectiveMassTest(plotIntegrand=True, plotSpectra=True)
    # model = gt.GETELECModel(emitterType="semiconductor", conductionBandBottom=0.1, bandGap=0.5, numberOfSpectrumPoints=256)
    # model.calculateElectronSpectrum()
    # spectra = model.getElectronSpectrum()
    # plt.plot(spectra["energy"][0], spectra["electronCount"][0])
    # plt.savefig("spectra.png")

    # tbTest = tabulatorTests()


    # tbTest.compareOldAndNew()

    tabulator = gt.GamowTabulator(Nf=16, Nr=1, XCdataFile="tabulated/XCdata_W110.npy")
    gt._setTabulationPath("tabulated/1D_512_new")

    tabulator.tabulateGamowTable()

