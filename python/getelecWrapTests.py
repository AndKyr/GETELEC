import getelec_wrap as gt

import numpy as np
import matplotlib.pyplot as plt

import unittest as tst

class TestGetelecInterface(tst.TestCase):
    def testDensityAndSpectra(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg", numberOfThreads=1)

        Ntests = 4

        # fields = np.random.rand(Ntests) * 10. + 1.
        # radii = 1./np.random.rand(Ntests)
        # gammas = np.random.random(Ntests) * 10. + 1.
        # kTs = np.random.random(Ntests) * 0.22
        # workFunctions = np.random.random(Ntests) * 3. + 2.
        # bandDepths = np.random.random(Ntests) * 3. + 10.
        # effectiveMasses = np.random.random(Ntests) * 1. + 1.
        getelec.setRandomParams(Ntests)

        
        
        # getelec.setParameters(field=fields, radius=radii, gamma=gammas, kT=kTs, workFunction=workFunctions, bandDepth=bandDepths, effectiveMass=effectiveMasses)
        getelec.run(["CurrentDensity", "TotalEnergyDistribution", "NormalEnergyDistribution", "ParallelEnergyDistribution", "NottinghamHeat"])

        allSpectra = getelec.extractAllSpectra()
        densities = getelec.getCurrentDensity()
        nottinghams = getelec.getNottinghamHeat()

        

        plotRows = 2
        plotColumns = 2
        plotsPerFigure = plotRows * plotColumns
        Nfigures = int(np.ceil(Ntests / plotsPerFigure))
        counter = 0
        for i in range(Nfigures):
            fig, axes = plt.subplots(plotRows, plotColumns, squeeze=True)

            for ax in axes:
                for spectraTypeKey in allSpectra:
                    if (spectraTypeKey == "TotalEnergyDistributionSplines"):
                        spline = allSpectra[spectraTypeKey][counter]
                        xPlot = np.linspace(min(spline.x), max(spline.x), 512)
                        yPlot = spline(xPlot)
                        ax.plot(xPlot, yPlot, label=spectraTypeKey)
                    else:
                        ax.plot(allSpectra[spectraTypeKey][counter][0], allSpectra[spectraTypeKey][counter][1], label=spectraTypeKey)
                        
                ax.set_xlabel("Energy [eV]")
                ax.set_ylabel("Spectra [A/nm^2/eV]")
            plt.show()
    
    def testIVcurve(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg")
        xFN = np.linspace(.1, .5, 512)
        getelec.setField(1./xFN)
        getelec.run()
        currentDensities = getelec.getCurrentDensity()
        poly = np.polyfit(xFN, np.log(currentDensities), 2)

        plt.semilogy(xFN, currentDensities, ".")
        plt.semilogy(xFN, np.exp(np.polyval(poly, xFN)))
        plt.xlabel("1/Field (nm/V)")
        plt.ylabel("Current Density (A/nm^2)")
        plt.show()

        # expected_poly = [12.40708813, -75.37623784, -4.9849858]
        expected_poly = [ 7.209221, -73.143162,  -5.085605]

        np.testing.assert_allclose(poly, expected_poly, rtol=1e-3)  


    

if (__name__ == "__main__"):
    # Example usage:
    
    tst.main()