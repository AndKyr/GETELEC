import numpy as np 
import py4vasp
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline  
from itertools import cycle
from getelec import Barrier, Globals, BandEmitter, ConductionBandEmitter

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


class dftBandEmitter(BandEmitter):
    
    def __init__(self, dftRunPath = ".") -> None:
        self.vaspData = py4vasp.Calculation.from_path(dftRunPath)
        self.kPointsRaw = self.vaspData.kpoint.to_dict()["coordinates"]
        self.bandsRaw = self.vaspData.band.to_dict()["bands"]
        self.reciprocalLatticeVectors = 2 * np.pi / np.diagonal(self.vaspData.structure.to_dict()["lattice_vectors"])
        self.weightsRaw = self.vaspData.kpoint.to_dict()["weights"]
        self.brillouinZoneVolume = np.prod(self.reciprocalLatticeVectors)
        #TODO: assert that the lattice is orthogonal
        self.NkPoints, self.nBands = np.shape(self.bandsRaw)
        self.getReshapedArrays()


    def getReshapedArrays(self):
        self.kxPoints = np.unique(self.kPointsRaw[:,0]) * self.reciprocalLatticeVectors[0]
        self.kyPoints = np.unique(self.kPointsRaw[:,1]) * self.reciprocalLatticeVectors[1]
        self.kzPoints = np.unique(self.kPointsRaw[:,2]) * self.reciprocalLatticeVectors[2]

        self.Nkx = len(self.kxPoints)
        self.Nky = len(self.kyPoints)
        self.Nkz = len(self.kzPoints)

        # assert self.Nkz == 2, "There must be exactly 2 kpoints in the z direction"
        
        # self.kPointsReshaped = np.reshape(self.kPoints, (len(self.kxPoints), len(self.kyPoints), 2,3))
        self.kPoints = np.reshape(self.kPointsRaw, (len(self.kzPoints), len(self.kyPoints), len(self.kxPoints), 3)) * self.reciprocalLatticeVectors
        self.bands = np.reshape(self.bandsRaw, (len(self.kzPoints), len(self.kyPoints), len(self.kxPoints), self.nBands))
        self.weights = np.copy(self.bands)
        for i in range(self.nBands):
            self.weights[:,:,:,i] = np.reshape(self.weightsRaw, (len(self.kzPoints), len(self.kyPoints), len(self.kxPoints)))

        self.zBandWidth = self.bands[1, :,:, :] - self.bands[0,:,:,:]
        # self.zBandWidth2 = self.bands[1, :,:, :] - self.bands[0,:,:,:]

        
    def plotKsurfaces(self):
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        for iBand in range(0, self.nBands, 5):
            surf = ax.plot_surface(self.kPoints[0, :, :, 0],self.kPoints[0, :, :, 1], self.bands[0,:,:,iBand], cmap=mb.cm.coolwarm)
        # fig.colorbar(surf)
        plt.show()  
        
    def get2DSplines(self):
        #TODO
        """get a list of 2D spline functions(kx, ky) for all bands"""
        
    def plotBands(self, mode = "totalEnergy", show = False):
        formats = ["-", ".", ":"]
        fig1, axesForDifferentX = plt.subplots(1,len(self.kxPoints))
        fig2, axesForDifferentY = plt.subplots(1, len(self.kyPoints))
        for axis in axesForDifferentY:
            axis.set_xlabel(r"$k_x [\textrm{\AA}]$")
        axesForDifferentY[0].set_ylabel(r"$E [\textrm{eV}]$")
        
        for axis in axesForDifferentX:
            axis.set_xlabel(r"$k_y [\textrm{\AA}]$")
        axesForDifferentX[0].set_ylabel(r"$E [\textrm{eV}]$")

        zRange = len(self.kzPoints)
        if (mode == "totalEnergy"):
            plotData = self.bands
            fig1.suptitle("Total Energy")
            fig2.suptitle("Total Energy")
        elif(mode == "parallelEnergy"):
            plotData = self.parallelEnergy
            fig1.suptitle("Parallel Energy")
            fig2.suptitle("Parallel Energy")
        elif (mode == "perpendicularEnergy"):
            plotData = self.perpendicularEnergy
            fig1.suptitle("Perpendicular Energy")
            fig2.suptitle("Perpendicular Energy")
        elif (mode == "energyDifference"):
            plotData = abs(np.array([self.zBandWidth]))
            zRange = 1
            fig1.suptitle("Energy difference in kz")
            fig2.suptitle("Energy difference on kz")
        else:
            assert False, "wrong mode. Should be 'totalEnergy', 'parallelEnergy', 'energyDifference' or 'perpendicularEnergy'"

        for iz in range(zRange):      
            for iy in range(len(self.kyPoints)):
                axis = axesForDifferentY[iy]
                axis.set_title(r"$k_y = %.2f \textrm{ \AA}^{-1}$"%self.kyPoints[iy])
                if (mode == "energyDifference"):
                    axis.set_yscale("log")
                # xPlot = np.linspace(0, max(self.kxPoints), 256)

                for iBand in range(self.nBands):
                    # spline = CubicSpline(self.kPoints[iz, iy, :, 0],  plotData[iz, iy, :, iBand], bc_type="clamped")
                    # derivative = spline.derivative()
                    color = next(colors)
                    axis.plot(self.kPoints[iz, iy, :, 0], plotData[iz, iy, :, iBand],formats[iz], color = color)
                    # axis.plot(xPlot, spline(xPlot), lines[iz], color = color)
                
        for iz in range(zRange):        
            for ix in range(len(self.kxPoints)):
                axis = axesForDifferentX[ix]
                axis.set_title(r"$k_x = %.2f \textrm{ \AA}^{-1}$"%self.kxPoints[ix])
                # xPlot = np.linspace(0., max(self.kyPoints), 256)
                if (mode == "energyDifference"):
                    axis.set_yscale("log")

                for iBand in range(self.nBands):
                    # spline = CubicSpline(self.kPoints[iz, :, ix, 1], plotData[iz, :, ix, iBand], bc_type="clamped")
                    # derivative = spline.derivative()

                    color = next(colors)
                    axis.plot(self.kPoints[iz, :, ix, 1], plotData[iz, :, ix, iBand], formats[iz], color = color)
                    
                    # axis.plot(xPlot, spline(xPlot), lines[iz], color = color)
        if (show):
            plt.show()
        
    def getBAndDerivativesGradient(self):
        self.partialDerivativesOfBands = np.array(np.gradient(self.bands, axis=(2,1,0)))
        self.partialDerivativesOfBands[0, :, :, 0, :] = 0.
        self.partialDerivativesOfBands[0, :, :, -1, :] = 0.
        self.partialDerivativesOfBands[1, :, 0, :, :] = 0.
        self.partialDerivativesOfBands[1, :, -1, :, :] = 0.


        self.partialDerivativesOfBands[0] /= (self.kxPoints[1] - self.kxPoints[0])
        self.partialDerivativesOfBands[1] /= (self.kyPoints[1] - self.kyPoints[0])
        self.partialDerivativesOfBands[2] /= (self.kzPoints[1] - self.kzPoints[0])

    def getCurrentOfStates(self, barrier: Barrier, workFunction = 5.25):
        self.currentDensityElements = barrier.transmissionCoefficient(workFunction - self.parallelEnergy[0, :,:,:]) * self.FermiDiracFunction(self.bands[0,:,:,:], self.kT) * abs(self.zBandWidth) * Globals.bandIntegrationPrefactor * self.weights[0,:,:,:] * self.reciprocalLatticeVectors[0] * self.reciprocalLatticeVectors[1]
        return np.sum(self.currentDensityElements)

    def getTotalEnergyDistribution(self):
        # flatEnergies = self.bands.flatten()
        # flatCurrents = self.currentIntegrand.flatten()
        # sortInds = np.argsort(flatEnergies)
        # self.sortedEnergies = flatEnergies[sortInds]
        # self.sortedCurrents = flatCurrents[sortInds]
        # self.cumulativeCurrent = np.cumsum(self.sortedCurrents)
        # self.totalEnergyDistribution = np.gradient(self.cumulativeCurrent) / np.gradient(self.sortedEnergies)
        self.totalEnergyDistribution, energyBinEdges = np.histogram(self.bands[0,:,:,:].flatten(), bins = 64, weights=self.currentDensityElements.flatten(), density=False, range=(-2., 0.5))
        self.energiesForDistribution = 0.5 * ( energyBinEdges[1:] + energyBinEdges[:-1])
        self.binWidth = energyBinEdges[1:] - energyBinEdges[:-1]
        self.totalEnergyDistribution /= self.binWidth

    def getDensityOfStates(self, bins = 64, range = (-1.5, 0.5)):
        self.densityOfStates, energyBinEdges = np.histogram(self.bands.flatten(), bins = bins, weights=self.weights.flatten(), density=False, range=range)
        self.dosEnergies = 0.5 * ( energyBinEdges[1:] + energyBinEdges[:-1])
        self.dosBinWidth = energyBinEdges[1:] - energyBinEdges[:-1]
        self.densityOfStates /= self.dosBinWidth
    


    def getBandDerivativesCspline(self):
        """Calculate an array of derivatives of energy with respect to kx,ky"""

        self.partialDerivativesOfBands = np.array(3*[np.copy(self.bands)])

        for iz in range(2):        
            for iy in range(len(self.kyPoints)):
                for iBand in range(self.nBands):
                    spline = CubicSpline(self.kPoints[iz, iy, :, 0], self.bands[iz, iy, :, iBand], bc_type="clamped")
                    derivative = spline.derivative()
                    self.partialDerivativesOfBands[0, iz,iy,:, iBand] = derivative(self.kPoints[iz, iy, :, 0])
                
        for iz in range(2):        
            for ix in range(len(self.kxPoints)):
                for iBand in range(self.nBands):
                    spline = CubicSpline(self.kPoints[iz, :, ix, 1], self.bands[iz, :, ix, iBand], bc_type="clamped")
                    derivative = spline.derivative()
                    self.partialDerivativesOfBands[1, iz, :, ix, iBand] = derivative(self.kPoints[iz, :, ix, 1])
        
        self.partialDerivativesOfBands[2, :, :, :, :] = (self.bands[1, :, :, :] - self.bands[0, :, :, :]) / (self.kzPoints[1] - self.kzPoints[0])

            
    def calculateParallelEnergy(self):
        # self.getBandDerivativesCspline()
        self.getBAndDerivativesGradient()
        halfMeOverHbarSquare = 0.06561710580766435 #eV^-1 Angstrom^-2
        self.perpendicularEnergy = halfMeOverHbarSquare * (self.partialDerivativesOfBands[0, :, :, :, :]**2 + self.partialDerivativesOfBands[1, :,  :, :, :]**2)
        self.parallelEnergy = self.bands - self.perpendicularEnergy


workFunction = 5.25
emitter = dftBandEmitter("/home/kyritsak/vasp_runs/W_surfaces/W110_slab/shortSlab")
emitter.kT = 0.025
emitter.calculateParallelEnergy()

barrier = Barrier(tabulationFolder="./tabulated/1D_1024", field=7.5)

currentDft = emitter.getCurrentOfStates(barrier, workFunction)
emitter.getTotalEnergyDistribution()

em2 = ConductionBandEmitter(barrier, workFunction=workFunction, kT = emitter.kT)
em2.calculateTotalEnergySpectrum()
emitter.getDensityOfStates(bins=128, range=(-1, 0.3))

# plt.semilogy(emitter.parallelEnergy[0].flatten(), emitter.currentDensityElements.flatten(), ".")
# plt.semilogy(emitter.parallelEnergy[0].flatten(), barrier.transmissionCoefficient(workFunction - emitter.parallelEnergy[0].flatten()), ".")

# plt.plot(emitter.dosEnergies, emitter.densityOfStates)
plt.plot(emitter.energiesForDistribution, emitter.totalEnergyDistribution)
plt.plot(em2.totalEnergySpectrumArrays()[0], em2.totalEnergySpectrumArrays()[1])
# plt.yscale("log")
plt.grid()

plt.show()

# plt.semilogy(emitter.bands[0].flatten(), abs(emitter.zBandWidth.flatten()), ".")
# plt.show()

print("total current density dft = ", currentDft, " getelec: ", em2.currentDensity())
# emitter.plotKsurfaces()
# emitter.plotBands(show=True, mode="totalEnergy")
# emitter.plotBands(show=True, mode="energyDifference")
# emitter.plotBands(show=True, mode="parallelEnergy")
# emitter.plotBands(show=True, mode="perpendicularEnergy")

# emitter.calculateZBandWidth()
# bandXY_0 = myBand[:15]
# bandXY_top = myBand[15:]

# print(len(bandXY_0))
# print(len(bandXY_top))

# kx = kPoints[:15,0]
# ky = kPoints[:15,1]


# fig = plt.figure()
# ax = fig.add_subplot(projection="3d")
# ax.scatter3D(kx,ky,bandXY_0, c=bandXY_0)
# # ax.plot3D(kx,ky,bandXY_top, ".")
plt.show()
