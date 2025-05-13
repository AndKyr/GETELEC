import getelec_wrap as gt

import numpy as np
import matplotlib.pyplot as plt

import unittest as tst

class TestGetelecInterface(tst.TestCase):
    def testDensityAndSpectra(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg", numberOfThreads=1)

        Ntests = 8

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

        results = getelec.getResultDictionary()

        

        plotRows = 2
        plotColumns = 2
        plotsPerFigure = plotRows * plotColumns
        Nfigures = int(np.ceil(Ntests / plotsPerFigure))
        counter = 0
        for i in range(Nfigures):
            fig, axes = plt.subplots(plotRows, plotColumns, squeeze=True, figsize = (28,20))

            for ax in axes.flatten():
                title = ""
                for key in results:
                    if ("Distribution" in key):
                        ax.plot(results[key][counter][0], results[key][counter][1], label=key)
                    else:
                        title += "%s=%.3e, " % (key, results[key][counter])

                ax.set_title(title)
                ax.legend()
                counter += 1      
                ax.set_xlabel("Energy [eV]")
                ax.set_ylabel("Spectra [A/nm^2/eV]")
                ax.grid()
            plt.show()
    
    def testIVcurve(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg")
        xFN = np.linspace(.1, .5, 512)
        getelec.setField(1./xFN)
        getelec.run(["CurrentDensity"])
        currentDensities = getelec.getCurrentDensity()
        print(currentDensities)
        poly = np.polyfit(xFN, np.log(currentDensities), 2)

        plt.semilogy(xFN, currentDensities, ".")
        plt.semilogy(xFN, np.exp(np.polyval(poly, xFN)))
        plt.xlabel("1/Field (nm/V)")
        plt.ylabel("Current Density (A/nm^2)")
        plt.show()

        # expected_poly = [12.40708813, -75.37623784, -4.9849858]
        expected_poly = [ 7.234612, -73.16,  -5.083]

        np.testing.assert_allclose(poly, expected_poly, rtol=1e-1)  


    

if (__name__ == "__main__"):
    # Example usage:
    
    tst.main()