import getelec_wrap as gt

import numpy as np
import matplotlib.pyplot as plt

import unittest as tst

class TestGetelecInterface(tst.TestCase):
    def testDensityAndSpectra(self):
        getelec = gt.GetelecInterface()
        getelec.setRadius([5., 6., 7. , 8., 9.])
        getelec.run(calculate_spectra=True)
        densities = getelec.getCurrentDensity()
        spectra = getelec.getSpectra()

        for spline, density in zip(spectra, densities):
            xPlot = np.linspace(min(spline.x), max(spline.x), 512)
            yPlot = spline(xPlot)
            currentDensityFromSpectra = np.trapz(yPlot, xPlot)
            self.assertAlmostEqual(np.log(currentDensityFromSpectra), np.log(density), delta=1e-3)
    
    def testIVcurve(self):
        getelec = gt.GetelecInterface()
        xFN = np.linspace(.1, .3, 512)
        getelec.setField(1./xFN)
        getelec.run()
        currentDensitiesLog = np.log(getelec.getCurrentDensity())
        poly = np.polyfit(xFN, currentDensitiesLog, 2)
        expected_poly = [12.40708813, -75.37623784, -4.9849858]
        np.testing.assert_allclose(poly, expected_poly, rtol=1e-3)        

    

if (__name__ == "__main__"):
    # Example usage:
    
    tst.main()