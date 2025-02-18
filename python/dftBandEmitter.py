import numpy as np 
import py4vasp
import matplotlib.pyplot as plt
import scipy.interpolate as interp  
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
                elif(mode == "parallelEnergy"):
                    pointValues = self.parallelEnergy[0, :, :, iBand].flatten()
                elif(mode == "perpendicularEnergy"):
                    pointValues = self.perpendicularEnergy[0, :, :, iBand].flatten()
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

    def fitBivariateSpline(self, X, Y, Z, weights = None, eps = 1.e-10):
        
        rmsValueOfZ = np.sqrt(np.mean(Z**2))
        smoothnessTolerances = np.logspace(-4, 3, 16) * rmsValueOfZ
        for smoothnessTolerance in smoothnessTolerances:
            try:
                spline = interp.SmoothBivariateSpline(X, Y, Z, w=weights, s=smoothnessTolerance, eps=eps, kx=4, ky=4)
                return spline
            except:
                pass

    def getBivariateSplinesForDeltaEz(self, verbose = False):
        """get a list of 2D spline functions(kx, ky) for all bands"""


        self.bivariateSplinesForEz = []
                                

        if (self.kPointSymmetry() >= 4):
            kPointsX = np.concatenate((kPointsX, kPointsX, -kPointsX, -kPointsX, 2*np.max(kPointsX)-kPointsX, kPointsX,                    2*np.max(kPointsX)-kPointsX, 2*np.max(kPointsX)-kPointsX, -kPointsX))
            kPointsY = np.concatenate((kPointsY, -kPointsY, kPointsY, -kPointsY, kPointsY,                    2*np.max(kPointsY)-kPointsY, 2*np.max(kPointsY)-kPointsY, -kPointsY, 2*np.max(kPointsY)-kPointsY))

        for iBand in range(self.nBands):
            energyDifference = self.bands[-1, :, iBand] - self.bands[0, :, iBand]
            if (self.kPointSymmetry() >= 8):
                energyDifference = np.concatenate((energyDifference, energyDifference))
            if (self.kPointSymmetry() >= 4):
                energyDifference = np.concatenate(tuple([energyDifference] * 9))

            newSpline = self.fitBivariateSpline(kPointsX, kPointsY, energyDifference, eps=1.e-5)
            if (verbose):
                print("got deltaEz spline for iBand, iz = ", iBand, "with relative error: ", np.sqrt(newSpline.get_residual() / np.sum(energyDifference**2)))

            self.bivariateSplinesForEz.append(newSpline)
        
    def getBivariateSplinesForBands(self, verbose = False):
        """get a list of 2D spline functions(kx, ky) for all bands"""

        self.bivariateSplinesForBands = []
                                
        if (self.kPointSymmetry() >= 8):
            kPointsX = np.concatenate((self.kPoints[0, :, 0], self.kPoints[0, :, 1]))
            kPointsY = np.concatenate((self.kPoints[0, :, 1],  self.kPoints[0, :, 0]))
        if (self.kPointSymmetry() >= 4):
            kPointsX = np.concatenate((kPointsX, kPointsX, -kPointsX, -kPointsX, 2*np.max(kPointsX)-kPointsX, kPointsX,                    2*np.max(kPointsX)-kPointsX, 2*np.max(kPointsX)-kPointsX, -kPointsX))
            kPointsY = np.concatenate((kPointsY, -kPointsY, kPointsY, -kPointsY, kPointsY,                    2*np.max(kPointsY)-kPointsY, 2*np.max(kPointsY)-kPointsY, -kPointsY, 2*np.max(kPointsY)-kPointsY))
            weights = np.concatenate(tuple([self.weights.flatten()] * 9))

        for iBand in range(self.nBands):
            bandSplines = []

            for iz in range(self.Nkz):
                energy = self.bands[iz, :, iBand]

                if (self.kPointSymmetry() >= 8):
                    energy = np.concatenate((energy, energy))
                if (self.kPointSymmetry() >= 4):
                    energy = np.concatenate(tuple([energy] * 9))

                newSpline = self.fitBivariateSpline(kPointsX, kPointsY, energy, weights=weights, eps=1.e-5)
                bandSplines.append(newSpline)
                if (verbose):
                    print("got Spline for iBand, iz = ", iBand, iz, "with relative error: ", np.sqrt(newSpline.get_residual() / np.sum(energy**2)))


            self.bivariateSplinesForBands.append(bandSplines)


    def getIntegrationMesh2D(self, Npoints = 128):
        kXPoints = np.linspace(0,np.max(self.kPoints[0, :, 0]), Npoints)
        kYPoints = np.linspace(0,np.max(self.kPoints[0, :, 1]), Npoints)
        self.kXgrid, self.kYgrid = np.meshgrid(kXPoints, kYPoints)
        self.weightsOnGrid = np.copy(self.kXgrid)
        self.weightsOnGrid[:,:] = 4 * np.gradient(kXPoints)[1] * np.gradient(kYPoints)[1]
        self.weightsOnGrid[0, :] /= 2.
        self.weightsOnGrid[: ,  0] /= 2.
        self.weightsOnGrid[-1, :] /= 2.
        self.weightsOnGrid[:, -1] /= 2.
        # self.weightsOnGrid *= np.gradient(kXPoints)[1] * np.gradient(kYPoints)[1]

    
    def calculateParallelEnergyOnGrid(self, NgridPoints = 128):
        self.getIntegrationMesh2D(NgridPoints)
        self.perpendicularEnergyOnGrid = np.ones((self.nBands, NgridPoints, NgridPoints))
        self.parallelEnergyOnGrid = np.copy(self.perpendicularEnergyOnGrid)
        self.totalEnergyOnGrid = np.copy(self.perpendicularEnergyOnGrid)
        self.deltaEzOnGrid = np.copy(self.perpendicularEnergyOnGrid)

        for iBand in range(self.nBands):
            self.totalEnergyOnGrid[iBand, :, :] = self.bivariateSplinesForBands[iBand][0](self.kXgrid, self.kYgrid, grid = False)

            dE_dkx = self.bivariateSplinesForBands[iBand][0](self.kXgrid, self.kYgrid, 1, 0, grid = False)
            dE_dky = self.bivariateSplinesForBands[iBand][0](self.kXgrid, self.kYgrid, 1, 0, grid = False)
            halfMeOverHbarSquare = 0.06561710580766435 #eV^-1 Angstrom^-2
            self.perpendicularEnergyOnGrid[iBand, :, :] = halfMeOverHbarSquare * (dE_dkx**2 + dE_dky**2)
            self.parallelEnergyOnGrid = self.totalEnergyOnGrid - self.perpendicularEnergyOnGrid
            self.deltaEzOnGrid[iBand, : , :] = self.bivariateSplinesForEz[iBand](self.kXgrid, self.kYgrid, grid = False)


    def getCurrentElementsOnGrid(self, barrier: Barrier, workFunction = 5.25, kT = 0.025):
        self.currentDensityElementsOnGrid = barrier.transmissionCoefficient(workFunction - self.parallelEnergyOnGrid) * \
              self.FermiDiracFunction(self.totalEnergyOnGrid, kT) * abs(self.deltaEzOnGrid) * Globals.bandIntegrationPrefactor * self.weightsOnGrid
        return np.sum(self.currentDensityElementsOnGrid)


    def plotBandContoursForBothKz(self, mode = "totalEnergy", show = True):

        for iFigure in range(6,0,-1):
            if (self.nBands % iFigure == 0):
                NplotsInRow = iFigure
                break
        
        Nfigures = int(self.nBands / NplotsInRow)

        for iFigure in range(Nfigures):
            fig, axes = plt.subplots(2, NplotsInRow)
            for jRow in range(NplotsInRow):
                iBand = NplotsInRow * iFigure + jRow
                for iz in range(len(self.kPointsZ)):
                    splineValues = self.bivariateSplinesForBands[iBand][iz](self.kXgrid, self.kYgrid, grid = False)
                    minValue = min(np.min(splineValues), np.min(self.bands[iz, :, iBand]))
                    maxValue = max(np.max(splineValues), np.max(self.bands[iz, :, iBand]))
                    contour = axes[iz, jRow].contourf(self.kXgrid, self.kYgrid, splineValues, levels = 64, cmap = "jet", vmin = minValue, vmax = maxValue)
                    cbar = plt.colorbar(contour)
                    cbar.ax.set_ylabel(r"Energy [eV]")
                    axes[iz, jRow].scatter(self.kPoints[0, :, 0], self.kPoints[0, :, 1],c =  self.bands[iz, :, iBand], cmap = "jet", vmin = minValue, vmax = maxValue )
                    # axes[iz, jRow].set_aspect("equal")
                    axes[iz, jRow].set_xlabel(r"$k_x [\textrm{\AA}^{-1}]$")
                    axes[iz, jRow].set_ylabel(r"$k_x [\textrm{\AA}^{-1}]$")
                    axes[iz, jRow].set_title("iBand = %d, iz = %d"%(iBand, iz))
            
            if (show):
                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())
                plt.show()

    def plotBandContoursForZeroKz(self, mode = "totalEnergy", show = True):

        for nPlotCandidate in range(8,0,-2):
            if (self.nBands % nPlotCandidate == 0):
                nPlots = nPlotCandidate
                break
        
        Nfigures = int(self.nBands / nPlots)

        for iFigure in range(Nfigures):
            fig, axes = plt.subplots(2, int(nPlots/2))
            for jPlot in range(nPlots):
                iBand = nPlots * iFigure + jPlot
                iRow = jPlot // 2
                iColumn = jPlot % 2

                if (mode == "totalEnergy"):
                    splineValues = self.totalEnergyOnGrid[iBand]
                    pointValues = self.bands[0, :, iBand]
                    error = np.sqrt(self.bivariateSplinesForBands[iBand][0].get_residual())
                elif(mode == "parallelEnergy"):
                    splineValues = self.parallelEnergyOnGrid[iBand]
                    pointValues = splineValues
                    error = 0
                elif(mode == "perpendicularEnergy"):
                    splineValues = self.perpendicularEnergyOnGrid[iBand]
                    pointValues = splineValues
                    error = 0
                elif(mode == "deltaEz"):
                    splineValues = self.deltaEzOnGrid[iBand]
                    pointValues = self.bands[-1, :, iBand] - self.bands[0, :, iBand]
                    error = np.sqrt(np.mean(self.bivariateSplinesForEz[iBand].get_residual()))
                else:
                    assert False, "wrong mode"

                minValue = min(np.min(splineValues), np.min(pointValues))
                maxValue = max(np.max(splineValues), np.max(pointValues))
                contour = axes[iColumn, iRow].contourf(self.kXgrid, self.kYgrid, splineValues, levels = 64, cmap = "jet", vmin = minValue, vmax = maxValue)

                if (mode == "totalEnergy" or mode == "deltaEz"):
                    axes[iColumn, iRow].scatter(self.kPoints[0, :, 0], self.kPoints[0, :, 1],c = pointValues , cmap = "jet", vmin = minValue, vmax = maxValue )

                mappable = plt.cm.ScalarMappable(cmap="jet")
                mappable.set_clim(minValue, maxValue)

                cbar = plt.colorbar(mappable, ax = axes[iColumn, iRow])
                cbar.ax.set_ylabel(r"%s [eV]"%mode)
                # axes[iz, jRow].set_aspect("equal")
                axes[iColumn, iRow].set_xlabel(r"$k_x [\textrm{\AA}^{-1}]$")
                axes[iColumn, iRow].set_ylabel(r"$k_x [\textrm{\AA}^{-1}]$")
                axes[iColumn, iRow].set_title("iBand=%d, res=%g"%(iBand, error))
            
            if (show):
                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())
                plt.show()
        
    def plotBands(self, mode = "totalEnergy", show = False):
        formats = ["-", ".", ":"]
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
                    axis.plot(self.kPoints[iz, iy, :, 0], plotData[iz, iy, :, iBand],formats[iz], color = color)
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
                    axis.plot(self.kPoints[iz, :, ix, 1], plotData[iz, :, ix, iBand], formats[iz], color = color)
                    
        plt.colorbar(plt.cm.ScalarMappable(cmap="jet"), ax = axis)
        if (show):
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

    def getCurrentOfStates(self, barrier: Barrier, workFunction = 5.25, kT = 0.025):
        weightsOnAllBands = np.copy(self.bands[0,:,:,:])
        for i in range(self.nBands):
            weightsOnAllBands[:,:,i] = self.weights[0,:,:]
        self.currentDensityElements = barrier.transmissionCoefficient(workFunction - self.parallelEnergy[0, :,:,:]) * self.FermiDiracFunction(self.bands[0,:,:,:], kT) * abs(self.zBandWidth) * Globals.bandIntegrationPrefactor * weightsOnAllBands * self.reciprocalLatticeVectors[0] * self.reciprocalLatticeVectors[1]
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
    

    def getBandDerivativesCspline(self):
        """Calculate an array of derivatives of energy with respect to kx,ky"""

        self.partialDerivativesOfBands = np.array(3*[np.copy(self.bands)])

        for iz in range(2):        
            for iy in range(len(self.kPointsY)):
                for iBand in range(self.nBands):
                    spline = interp.CubicSpline(self.kPoints[iz, iy, :, 0], self.bands[iz, iy, :, iBand], bc_type="clamped")
                    derivative = spline.derivative()
                    self.partialDerivativesOfBands[0, iz,iy,:, iBand] = derivative(self.kPoints[iz, iy, :, 0])
                
        for iz in range(2):        
            for ix in range(len(self.kPointsX)):
                for iBand in range(self.nBands):
                    spline = interp.CubicSpline(self.kPoints[iz, :, ix, 1], self.bands[iz, :, ix, iBand], bc_type="clamped")
                    derivative = spline.derivative()
                    self.partialDerivativesOfBands[1, iz, :, ix, iBand] = derivative(self.kPoints[iz, :, ix, 1])
        
        self.partialDerivativesOfBands[2, :, :, :, :] = (self.bands[1, :, :, :] - self.bands[0, :, :, :]) / (self.kPointsZ[1] - self.kPointsZ[0])

            
    def calculateParallelEnergy(self):
        # self.getBandDerivativesCspline()
        self.getBAndDerivativesGradient()
        halfMeOverHbarSquare = 0.06561710580766435 #eV^-1 Angstrom^-2
        self.perpendicularEnergy = halfMeOverHbarSquare * (self.partialDerivativesOfBands[0, :, :, :, :]**2 + self.partialDerivativesOfBands[1, :,  :, :, :]**2)
        self.parallelEnergy = self.bands - self.perpendicularEnergy


workFunction = 5.25
kT = 0.025
emitter = dftBandEmitter("/home/kyritsak/vasp_runs/W_surfaces/W100_slab/smallSlab/16_16_2")

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
