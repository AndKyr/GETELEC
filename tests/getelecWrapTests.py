
#! 
from python import getelec_wrap as gt
from python import dataFitterIV as fitter


import numpy as np
import matplotlib.pyplot as plt

import unittest as tst
from itertools import cycle

font = 15
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
mb.rcParams["lines.markersize"] = 10
figureSize = [16,10]
colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

class TestGetelecInterface(tst.TestCase):
    def testDensityAndSpectra(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg", numberOfThreads=1)

        Ntests = 8
        getelec.setRandomParams(Ntests)
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
                        ax.semilogy(results[key][counter][0], results[key][counter][1], label=key)
                    else:
                        title += "%s=%.3e, " % (key, results[key][counter])

                ax.set_title(title)
                ax.legend()
                counter += 1      
                ax.set_xlabel(r"Energy [eV]")
                ax.set_ylabel(r"Spectra [A/nm$^2$/eV]")
                ax.grid()
            plt.show()
    
    def testIVcurve(self):
        getelec = gt.GetelecInterface(configPath="getelec.cfg")
        xFN = np.linspace(.1, .5, 512)
        getelec.setField(1./xFN)
        getelec.run(["CurrentDensity"])
        currentDensities = getelec.getCurrentDensity()
        poly = np.polyfit(xFN, np.log(currentDensities), 2)

        plt.semilogy(xFN, currentDensities, ".")
        plt.semilogy(xFN, np.exp(np.polyval(poly, xFN)))
        plt.xlabel("1/Field (nm/V)")
        plt.ylabel("Current Density (A/nm2)")
        plt.show()
        expected_poly = [ 7.234612, -73.16,  -5.083]
        np.testing.assert_allclose(poly, expected_poly, rtol=1e-1)  
    
    def testFitter(self):

        voltageData = np.array([6666.66666667, 6528.49740933, 6395.93908629, 6268.65671642, \
        6146.34146341, 6028.70813397, 5915.49295775, 5806.4516129, 5701.35746606, 5600., 5502.18340611, 5407.72532189,\
        5316.4556962 , 5228.21576763, 5142.85714286, 5060.24096386,\
        4980.23715415, 4902.72373541, 4827.5862069 , 4754.71698113,\
        4684.01486989, 4615.38461538, 4548.73646209, 4483.98576512,\
        4421.05263158, 4359.8615917 , 4300.34129693, 4242.42424242,\
        4186.04651163, 4131.14754098, 4077.66990291, 4025.55910543,\
        3974.76340694, 3925.23364486, 3876.92307692, 3829.78723404,\
        3783.78378378, 3738.87240356, 3695.01466276, 3652.17391304,\
        3610.31518625, 3569.40509915, 3529.41176471, 3490.30470914,\
        3452.05479452, 3414.63414634, 3378.01608579, 3342.17506631,\
        3307.08661417, 3272.72727273, 3239.07455013, 3206.10687023,\
        3173.80352645, 3142.1446384 , 3111.11111111, 3080.68459658,\
        3050.84745763, 3021.58273381, 2992.87410926, 2964.70588235,\
        2937.06293706, 2909.93071594, 2883.29519451, 2857.14285714])
    
        currentData = np.array([5.65354016e-05, 4.28464283e-05, 3.24541251e-05, 2.45683389e-05, \
            1.85875533e-05, 1.40539879e-05, 1.06193664e-05, 8.01882052e-06,\
            6.05098032e-06, 4.56287757e-06, 3.43825950e-06, 2.58892706e-06,\
            1.94793392e-06, 1.46451644e-06, 1.10020423e-06, 8.25856631e-07,\
            6.19415049e-07, 4.64192603e-07, 3.47574447e-07, 2.60030670e-07,\
            1.94366991e-07, 1.45156296e-07, 1.08307623e-07, 8.07395611e-08,\
            6.01329743e-08, 4.47437594e-08, 3.32614033e-08, 2.47019853e-08,\
            1.83274191e-08, 1.35845100e-08, 1.00590041e-08, 7.44096972e-09,\
            5.49873640e-09, 4.05929377e-09, 2.99355710e-09, 2.20530977e-09,\
            1.62290148e-09, 1.19302882e-09, 8.76074733e-10, 6.42626597e-10,\
            4.70867981e-10, 3.44633799e-10, 2.51960496e-10, 1.83999483e-10,\
            1.34216889e-10, 9.77911792e-11, 7.11688901e-11, 5.17338317e-11,\
            3.75620267e-11, 2.72401384e-11, 1.97311118e-11, 1.42748585e-11,\
            1.03149269e-11, 7.44441693e-12, 5.36613434e-12, 3.86326821e-12,\
            2.77783297e-12, 1.99486455e-12, 1.43077440e-12, 1.02488643e-12,\
            7.33204336e-13, 5.23857148e-13, 3.73798141e-13, 2.66375053e-13])


        outData = fitter.performFullEmissionFitting(voltageData, currentData, radius=[1., 5., 50.])

        fittedCurrent = outData["fittedCurrents"]
        plt.semilogy(1./voltageData, currentData, '.', label="data" )
        plt.semilogy(outData["fieldConversionFactor"] / outData["plotFields"] , fittedCurrent, label = "beta=%.2g, err=%.2g"%(outData["fieldConversionFactor"], outData["fittingError"]))

        plt.show()

    def testTransmissionProbability(self):
        # getelec = gt.GetelecInterface(configPath="getelec.cfg")
        # energies = np.linspace(-2., 2., 64)
        # waveVectors = 12. * np.ones(64)
        # coeffs = getelec.calculateTransmissionCoefficients(energies, waveVectors)

        # plt.semilogy(energies, np.abs(coeffs)**2)
        # plt.show()


        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        waveVectors = np.linspace(1.e-5, 20., 32)
        energies = np.zeros(len(waveVectors))

        emitter = gt.GetelecInterface(configPath="getelec.cfg")
        coefficients = emitter.calculateTransmissionCoefficients(energies, waveVectors)
        fig, ax1 = plt.subplots(figsize = figureSize, tight_layout=True)
        ax2 = ax1.twinx()

        ax1.plot(waveVectors, np.abs(coefficients), label=r"$|T|$", color = colors[0])
        ax2.plot(waveVectors, np.angle(coefficients), label=r"$\arg(T)$", color = colors[1])

        ax1.set_xlabel(r"$k_z \textrm{[nm}^{-1}\textrm{]}$")
        ax1.set_ylabel(r"$|T|$")
        ax2.set_ylabel(r"$\arg(T)$")
        ax1.grid()
        ax1.legend(loc="upper left")
        ax2.legend(loc="upper right")


        plt.savefig("solverComparison.pdf")
        plt.show()



    

if (__name__ == "__main__"):
    # Example usage:
    
    tst.main()