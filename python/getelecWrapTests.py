import getelec_wrap as gt

import numpy as np
import matplotlib.pyplot as plt

import unittest as tst

class TestGetelecInterface(tst.TestCase):
    def testDensityAndSpectra(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg")
        getelec.setRadius([5., 6., 7. , 8., 9.])
        getelec.run(calculate_spectra=True)
        densities = getelec.getCurrentDensity()
        spectra = getelec.getSpectra()

        for spline, density in zip(spectra, densities):
            xPlot = np.linspace(min(spline.x), max(spline.x), 512)
            yPlot = spline(xPlot)
            plt.plot(xPlot, yPlot)
            currentDensityFromSpectra = spline.integrate(min(spline.x), max(spline.x))
            self.assertAlmostEqual(np.log(currentDensityFromSpectra), np.log(density), delta=1e-3)
        plt.xlabel("Energy (eV)")
        plt.ylabel("spectra [A/nm^2/eV]")
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
        expected_poly = [7.2848  , -73.159786,  -5.200473]

        np.testing.assert_allclose(poly, expected_poly, rtol=1e-3)  


    

if (__name__ == "__main__"):
    # Example usage:
    
    tst.main()