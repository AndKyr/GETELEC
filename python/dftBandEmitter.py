import numpy as np 
import py4vasp
import matplotlib.pyplot as plt
import scipy.interpolate as interp  
from itertools import cycle
from getelec_wrap import GetelecInterface

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


class Globals:
    """Keeps global constants and variables"""
    BoltzmannConstant:float = 8.617333262e-5
    SommerfeldConstant:float = 1.618311e-4 
    electronMass:float = 9.1093837e-31
    imageChargeConstant:float = .359991137
    gamowPrefactor:float = 10.24633444
    bandIntegrationPrefactor:float = 4.86826961e-2  #A/nm^2 (energy in eV, k in 1/Angstrom)

class dftBandEmitter:
    
    def __init__(self, dftRunPath = ".") -> None:
        self.vaspData = py4vasp.Calculation.from_path(dftRunPath)
        self.kPointsRaw = self.vaspData.kpoint.to_dict()["coordinates"]
        self.bandsRaw = self.vaspData.band.to_dict()["bands"]
        self.reciprocalLatticeVectors = 2 * np.pi / np.diagonal(self.vaspData.structure.to_dict()["lattice_vectors"])
        self.weightsRaw = self.vaspData.kpoint.to_dict()["weights"]
        self.brillouinZoneVolume = np.prod(self.reciprocalLatticeVectors)
        self.NkPoints, self.nBands = np.shape(self.bandsRaw)
        self.getReshapedArrays()


    def getReshapedArrays(self):
        # self.kPoints = self.kPointsRaw# * self.reciprocalLatticeVectors
        self.kPointsX = np.unique(self.kPointsRaw[:,0])
        self.kPointsY = np.unique(self.kPointsRaw[:,1])
        self.kPointsZ = np.unique(self.kPointsRaw[:,2])

        self.kPoints = np.flip(np.array(np.meshgrid(self.kPointsZ, self.kPointsY, self.kPointsX, indexing="ij")), axis = 0)
        self.Nkz = len(self.kPointsZ)

        self.NkXY = int(len(self.kPointsRaw) / len(self.kPointsZ))


        self.bands = np.ones(np.shape(self.kPoints[0]) + (self.nBands,))
        self.weights = np.copy(self.kPoints[0])

        if (self.kPointSymmetry() > 7.9):
            kPointsExtended = np.copy(self.kPointsRaw)
            kPointsExtended[:,0] = self.kPointsRaw[:,1]
            kPointsExtended[:,1] = self.kPointsRaw[:,0]
            self.kPointsRaw = np.concatenate((self.kPointsRaw, kPointsExtended))

        for i in range(len(self.kPointsZ)):
            zIndices  = np.nonzero(abs(self.kPointsRaw[:,2] - self.kPointsZ[i]) < 1.e-8)[0]
            for j in range(len(self.kPointsY)):
                zyIndices = zIndices[abs(self.kPointsRaw[zIndices,1] - self.kPointsY[j]) < 1.e-8]
                for k in range(len(self.kPointsX)):
                    kPointIndex = np.min(zyIndices[abs(self.kPointsRaw[zyIndices,0] - self.kPointsX[k]) < 1.e-8])
                    self.bands[i,j,k,:] = self.bandsRaw[kPointIndex % self.NkPoints, :]
                    self.weights[i,j,k] = self.weightsRaw[kPointIndex % self.NkPoints]


        # self.kPoints = np.reshape(self.kPointsRaw, (len(self.kPointsZ), self.NkXY, 3)) * self.reciprocalLatticeVectors
        # self.bands = np.reshape(self.bandsRaw, (len(self.kPointsZ),self.NkXY, self.nBands))
        # self.weights = np.reshape(self.weightsRaw, (len(self.kPointsZ),self.NkXY))


        self.zBandWidth = self.bands[-1, :,:, :] - self.bands[0,:,:,:]
        # self.zBandWidth2 = self.bands[1, :,:, :] - self.bands[0,:,:,:]

    def findNumberOfPlotsPerFigure(self):
        for nPlotCandidate in range(8,0,-2):
            if (self.nBands % nPlotCandidate == 0):
                nPlots = nPlotCandidate
                break
        
        return int(self.nBands / nPlots), nPlots
        
    def colorPlotBandsOnkPoints(self, mode = "totalEnergy", show = True):
        
        Nfigures, nPlots = self.findNumberOfPlotsPerFigure()

        for iFigure in range(Nfigures):
            fig, axes = plt.subplots(2, int(nPlots/2))
            for jPlot in range(nPlots):
                iBand = nPlots * iFigure + jPlot
                iRow = jPlot // 2
                iColumn = jPlot % 2

                if (mode == "totalEnergy"):
                    pointValues = self.bands[0, :, :, iBand].flatten()
                elif(mode == "normalEnergy"):
                    pointValues = self.normalEnergy[0, :, :, iBand].flatten()
                elif(mode == "parallelEnergy"):
                    pointValues = self.parallelEnergy[0, :, :, iBand].flatten()
                elif(mode == "deltaEz"):
                    pointValues = self.bands[-1, :, :, iBand].flatten() - self.bands[0, :, :, iBand].flatten()
                else:
                    assert False, "wrong mode"

                scatterPlot = axes[iColumn, iRow].scatter(self.kPoints[0, 0, :, :].flatten(), self.kPoints[1, 0, :, :].flatten(), c = pointValues , cmap = "jet")

                cbar = plt.colorbar(scatterPlot)
                cbar.ax.set_ylabel(r"%s [eV]"%mode)
                # axes[iz, jRow].set_aspect("equal")
                axes[iColumn, iRow].set_xlabel(r"$k_x [\textrm{\AA}^{-1}]$")
                axes[iColumn, iRow].set_ylabel(r"$k_y [\textrm{\AA}^{-1}]$")
                axes[iColumn, iRow].set_title("iBand=%d"%(iBand))
            
            if (show):
                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())
                plt.show()


    def plotBands(self, mode = "totalEnergy", show = False):
        formats = ["-", ".", ":", "--"]
        colors = plt.cm.jet(np.linspace(0,1,self.nBands))
        fig1, axesForDifferentX = plt.subplots(1,len(self.kPointsX))
        fig2, axesForDifferentY = plt.subplots(1, len(self.kPointsY))
        for axis in axesForDifferentY:
            axis.set_xlabel(r"$k_x [\textrm{\AA}^{-1}]$")
        axesForDifferentY[0].set_ylabel(r"$E [\textrm{eV}]$")
        
        for axis in axesForDifferentX:
            axis.set_xlabel(r"$k_y [\textrm{\AA}]$")
        axesForDifferentX[0].set_ylabel(r"$E [\textrm{eV}]$")

        zRange = len(self.kPointsZ)
        if (mode == "totalEnergy"):
            plotData = self.bands
            fig1.suptitle("Total Energy")
            fig2.suptitle("Total Energy")
        elif(mode == "normalEnergy"):
            plotData = self.normalEnergy
            fig1.suptitle("Normal Energy")
            fig2.suptitle("Nromal Energy")
        elif (mode == "parallelEnergy"):
            plotData = self.parallelEnergy
            fig1.suptitle("Parallel Energy")
            fig2.suptitle("Parallel Energy")
        elif (mode == "deltaEz"):
            plotData = abs(np.array([self.zBandWidth]))
            zRange = 1
            fig1.suptitle("Energy difference in kz")
            fig2.suptitle("Energy difference on kz")
        else:
            assert False, "wrong mode. Should be 'totalEnergy', 'normalEnergy', 'parallelEnergy' or 'deltaEz'"

        for iz in range(zRange):      
            for iy in range(len(self.kPointsY)):
                axis = axesForDifferentY[iy]
                axis.set_title(r"$k_y = %.2f \textrm{ \AA}^{-1}$"%self.kPointsY[iy])
                if (mode == "energyDifference"):
                    axis.set_yscale("log")
                # xPlot = np.linspace(0, max(self.kxPoints), 256)

                for iBand in range(self.nBands):
                    # spline = CubicSpline(self.kPoints[iz, iy, :, 0],  plotData[iz, iy, :, iBand], bc_type="clamped")
                    # derivative = spline.derivative()
                    # color = next(colors)
                    color = colors[iBand]
                    axis.plot(self.kPoints[0, iz, iy, :], plotData[iz, iy, :, iBand],formats[iz], color = color)
                    # axis.plot(xPlot, spline(xPlot), lines[iz], color = color)

        plt.colorbar(plt.cm.ScalarMappable(cmap="jet"), ax = axis)
        for iz in range(zRange):        
            for ix in range(len(self.kPointsX)):
                axis = axesForDifferentX[ix]
                axis.set_title(r"$k_x = %.2f \textrm{ \AA}^{-1}$"%self.kPointsX[ix])
                # xPlot = np.linspace(0., max(self.kyPoints), 256)
                if (mode == "energyDifference"):
                    axis.set_yscale("log")

                for iBand in range(self.nBands):
                    # spline = CubicSpline(self.kPoints[iz, :, ix, 1], plotData[iz, :, ix, iBand], bc_type="clamped")
                    # derivative = spline.derivative()

                    # color = next(colors)
                    color = colors[iBand]
                    axis.plot(self.kPoints[1, iz, :, ix], plotData[iz, :, ix, iBand], formats[iz], color = color)
                    
        plt.colorbar(plt.cm.ScalarMappable(cmap="jet"), ax = axis)
        if (show):
            plt.show()

    def plotBandsOnKz(self):
        nPlots = self.nBands  * len(self.kPointsX) * len(self.kPointsY)
        # colors = plt.cm.jet(np.linspace(0,1,nPlots))  
        # maxEnergyOnKz0 = np.max(self.bands[0,:,:,:].flatten())
        # minEnergyOnKz0 = np.min(self.bands[0,:,:,:].flatten())
        colorNumberArray = np.linspace(0, 1, 128)
        colors = plt.cm.jet(colorNumberArray)
        # fig2, axesForDifferentY = plt.subplots(1, len(self.kPointsY))
        fig, axes = plt.subplots(4,5, squeeze=True)
        fig.tight_layout()
        iPlot = 0
        for ix in range(len(self.kPointsX)):
            for iy in range(len(self.kPointsY)): 
                for iBand in range(self.nBands):
                    if (self.bands[0, iy, ix, iBand] < 0.1 and self.bands[0, iy, ix, iBand] > -1.):
                        axes.flatten()[iPlot % 20].plot(self.kPoints[2, :, iy, ix], self.bands[:, iy, ix, iBand] - self.bands[0, iy, ix, iBand], color = colors[iPlot % 128])
                        iPlot += 1
        
        # plt.colorbar(plt.cm.ScalarMappable(cmap="jet"), ax = ax)
        plt.show()
        
    def getBAndDerivativesGradient(self):
        self.partialDerivativesOfBands = np.array(np.gradient(self.bands, axis=(2,1,0)))
        # self.partialDerivativesOfBands[0, :, :, 0, :] = 0.
        # self.partialDerivativesOfBands[0, :, :, -1, :] = 0.
        # self.partialDerivativesOfBands[1, :, 0, :, :] = 0.
        # self.partialDerivativesOfBands[1, :, -1, :, :] = 0.


        self.partialDerivativesOfBands[0] /= (self.kPointsX[1] - self.kPointsX[0])
        self.partialDerivativesOfBands[1] /= (self.kPointsY[1] - self.kPointsY[0])
        self.partialDerivativesOfBands[2] /= (self.kPointsZ[1] - self.kPointsZ[0])

    def getCurrentOfStates(self, getelecInterface: GetelecInterface, workFunction = 5.25, kT = 0.025):
        weightsOnAllBands = np.copy(self.bands[0,:,:,:])
        for i in range(self.nBands):
            weightsOnAllBands[:,:,i] = self.weights[0,:,:]
        
        self.currentDensityElements = getelecInterface.calculateTransmissionForManyEnergies(self.normalEnergy) * abs(self.zBandWidth) * Globals.bandIntegrationPrefactor * weightsOnAllBands * self.reciprocalLatticeVectors[0] * self.reciprocalLatticeVectors[1]
        return np.sum(self.currentDensityElements)

    def getTotalEnergyDistribution(self, bins = 64, range = (-1.5, 0.5), mode = "grid"):

        if (mode == "kpoints"):
            energyPoints = self.bands[0,:,:,:].flatten()
            currentDensityPoints = self.currentDensityElements.flatten()
        elif(mode == "grid"):
            energyPoints = self.totalEnergyOnGrid.flatten()
            currentDensityPoints = self.currentDensityElementsOnGrid.flatten()

        self.totalEnergyDistribution, energyBinEdges = np.histogram(energyPoints, bins = bins, weights=currentDensityPoints, density=False, range=range)
        self.energiesForDistribution = 0.5 * ( energyBinEdges[1:] + energyBinEdges[:-1])
        self.binWidth = energyBinEdges[1:] - energyBinEdges[:-1]
        self.totalEnergyDistribution /= self.binWidth
    
    def getCumulativeEnergyDistribution(self):
        flatEnergies = self.bands[0].flatten()
        # flatCurrents = self.currentIntegrand.flatten()
        sortInds = np.argsort(flatEnergies)
        self.sortedEnergies = flatEnergies[sortInds]
        self.sortedCurrentElements = self.currentDensityElements.flatten()[sortInds]
        self.cumulativeCurrent = np.cumsum(self.sortedCurrentElements)

    def getDensityOfStates(self, bins = 64, range = (-1.5, 0.5)):
        self.densityOfStates, energyBinEdges = np.histogram(self.bands.flatten(), bins = bins, weights=self.weights.flatten(), density=False, range=range)
        self.dosEnergies = 0.5 * ( energyBinEdges[1:] + energyBinEdges[:-1])
        self.dosBinWidth = energyBinEdges[1:] - energyBinEdges[:-1]
        self.densityOfStates /= self.dosBinWidth

    def kPointSymmetry(self):
        multiplicities = self.weightsRaw / np.min(self.weightsRaw)
        return max(multiplicities)
            
    def calculateParallelEnergy(self):
        # self.getBandDerivativesCspline()
        self.parallelEnergy = np.copy(self.kPoints)
        hbarSquareOver2m = 3.80998211 # Å^2 eV
        
        self.parallelEnergy = np.copy(self.bands)
        for i in range(self.nBands):
            self.parallelEnergy[:,:,:, i] = hbarSquareOver2m * (self.kPoints[0]**2 + self.kPoints[1]**2)
        self.normalEnergy = self.bands - self.parallelEnergy



workFunction = 5.25
kT = 0.025
emitter = dftBandEmitter("/home/kyritsak/vasp_runs/manyKz")
# fig = emitter.vaspData.band.plot()
# fig.show()
# emitter.plotBands(show=True)
emitter.plotBandsOnKz()
exit()
emitter.calculateParallelEnergy()
emitter.colorPlotBandsOnkPoints(mode="deltaEz")
exit()
barrier = Barrier(tabulationFolder="./tabulated/1D_1024", field=7)

em2 = ConductionBandEmitter(barrier, workFunction=workFunction, kT = kT)
em2.calculateTotalEnergySpectrum()
# emitter.getBivariateSplinesForBands(verbose=True)
# emitter.getBivariateSplinesForDeltaEz(verbose=True)
# emitter.calculateParallelEnergyOnGrid()
# currentDft = emitter.getCurrentElementsOnGrid(barrier=barrier)
# emitter.plotBandContours(show=True).
# plt.plot(emitter.totalEnergyOnGrid.flatten(), emitter.currentDensityElementsOnGrid.flatten(), ".")

# for mode in ["totalEnergy", "parallelEnergy", "perpendicularEnergy", "deltaEz"]:
#     emitter.plotBandContoursForZeroKz(mode = mode, show=True)

emitter.calculateParallelEnergy()
# for mode in ["totalEnergy", "parallelEnergy", "perpendicularEnergy", "deltaEz"]:
#     emitter.colorPlotBandsOnkPoints(mode=mode)

currentDft = emitter.getCurrentOfStates(barrier, workFunction, kT=kT)
emitter.getTotalEnergyDistribution(bins=32, range=(-1.5, 0.3), mode="kpoints")
# emitter.getCumulativeEnergyDistribution()



# plt.plot(emitter.bands[0].flatten(), emitter.currentDensityElements.flatten(), ".")

# plt.plot(emitter.dosEnergies, emitter.densityOfStates)
plt.bar(emitter.energiesForDistribution, emitter.totalEnergyDistribution, width=0.95 * np.gradient(emitter.energiesForDistribution)[1])
# plt.plot(emitter.energiesForDistribution, emitter.totalEnergyDistribution)
# gtEnergies, gtCurrents = em2.totalEnergySpectrumArrays()


# plt.plot(emitter.sortedEnergies, emitter.cumulativeCurrent)
# plt.plot(gtEnergies, np.cumsum(gtCurrents) * np.gradient(gtEnergies))


plt.plot(em2.totalEnergySpectrumArrays()[0], em2.totalEnergySpectrumArrays()[1])
# plt.yscale("log")
plt.grid()


print("total current density dft = ", currentDft, " getelec: ", em2.currentDensity())

plt.show()
