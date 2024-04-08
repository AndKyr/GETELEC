import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os
from getelec import Globals
import scipy.optimize as opt
import scipy.integrate as ig
import scipy.special
import scipy.linalg as linalg


font = 15
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
figureSize = [16,10]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

class Schrodinger1DSolverFDM:
    def __init__(self, potential:np.ndarray = np.ones(254), dx:float = 1., energy:float = 0.5):

        assert max(potential[0], potential[-1]) < energy, "the energy has to be more than the potential on the left and right edges"
        self.Nx = len(potential)
        self.length = self.Nx * dx
        self.dx = dx
        self.hbarSqrOver2m = 3.80998212 # Å^2 (square Ångströms) (eV (electronvolt))
        self.MConstant = dx**2 / self.hbarSqrOver2m # Å^2 (square Ångströms) (eV (electronvolt))
        self.potential = potential
        self.energy = energy
        self.calculateWaveVectors()


    def calculateWaveVectors(self):
        self.kLeft = np.sqrt((self.energy - self.potential[0]) / self.hbarSqrOver2m)
        self.kRight = np.sqrt((self.energy - self.potential[-1]) / self.hbarSqrOver2m)


    def setSparseMatrixSystem(self):
        #sets the sparse matrix for solving the implicit FDM problem
        mainDiagonal = np.ones(self.Nx + 2, dtype=complex)
        mainDiagonal[2:-2] = 2. + self.MConstant * (self.potential[1:-1] - self.energy) # -ui+1 - ui-1 + 2 ui + dx^2 (2m/hbar^2) (Vi - E) ui = 0
        mainDiagonal[0] = -1. # u0 - r = 1
        mainDiagonal[1] = -1. # u1 - u0 + i k dx r = i k dx
        mainDiagonal[-2] = 1. # uN-1 - uN-2 - i k dx t = 0 
        mainDiagonal[-1] = -1. #  uN-1 - t = 0
        
        topDiagonal = -np.ones(self.Nx + 1, dtype=complex) # ui+1 + ui-1 - 2 ui - dx^2 (2m/hbar^2) (Vi - E) ui = 0
        topDiagonal[1] = 1. # u1 - u0 + i k dx r = i k dx
        topDiagonal[0] = 1. # u0 - r = 1
        topDiagonal[-1] = -1j * self.kRight * self.dx # # uN-1 - uN-2 - i k dx t = 0 

        bottomDiagonal = -np.ones(self.Nx + 1, dtype=complex) # ui+1 + ui-1 - 2 ui - dx^2 (2m/hbar^2) (Vi - E) ui = 0
        bottomDiagonal[0] = 1j * self.kLeft * self.dx # u1 - u0 + i k dx r = i k dx
        bottomDiagonal[-2] = -1. # uN-1 - uN-2 - i k dx t = 0 
        bottomDiagonal[-1] = 1.  #  uN-1 - t = 0
        
        self.sparseMatrix = sp.diags([bottomDiagonal, mainDiagonal, topDiagonal], offsets=[-1,0,1],format="lil")
        self.rightHandSide = np.zeros(self.Nx + 2, dtype=complex)
        self.rightHandSide[0] = 1.
        self.rightHandSide[1] = 1j * self.kLeft * self.dx

    def solveFirstTime(self):
        self.setSparseMatrixSystem()
        self.sparseMatrixCSC = sp.csc_matrix(self.sparseMatrix)
        self.solution = sp.linalg.spsolve(self.sparseMatrixCSC, self.rightHandSide)
    
    def solveWithNewEnergy(self, energy):
        self.energy = energy
        
        self.calculateWaveVectors()

        self.setSparseMatrixSystem()
        self.sparseMatrixCSC = sp.csc_matrix(self.sparseMatrix)

        self.solution = sp.linalg.spsolve(self.sparseMatrixCSC, self.rightHandSide)

        # self.solution, info = sp.linalg.bicg(self.sparseMatrixCSC, self.rightHandSide, x0=self.solution, tol = 1.e-16)
        # print(info)

    def getTransmissionProbability(self):
        return np.abs(self.solution[-1])**2
    
    def calculateTransmissionForEnergyRange(self, Emin, Emax, NoEnergies):
        self.energies = np.linspace(Emin, Emax, NoEnergies)
        self.transmissionCoefficients = np.copy(self.energies)
        self.solveFirstTime()
        for i in range(len(self.energies)):
            self.solveWithNewEnergy(self.energies[i])
            self.transmissionCoefficients[i] = self.getTransmissionProbability()
        return self.energies, self.transmissionCoefficients

    def plotBarrier(self):
        plt.figure(figsize=figureSize)
        fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
        z = np.linspace(0., (self.Nx - 1) * self.dx, self.Nx)
        ax1.plot(z, self.potential, label="U")
        ax1.plot(z, np.ones(self.Nx) * self.energy, label="E")
        ax1.grid()
        ax1.legend()
        ax1.set_ylabel("barrier [eV]")

        ax2.plot(z, np.abs(self.solution[1:-1]), label=r"|$\Psi$|", color = colors[2])
        ax2r = ax2.twinx()
        ax2r.plot(z, np.angle(self.solution[1:-1]), label=r"Phase", color=colors[1])
        ax2.legend()
        ax2.grid()
        ax2.set_ylabel("$\Psi$")
        ax2r.set_ylabel("Angle ($\Psi$)")
        ax2r.legend()
        ax2.set_xlabel("z [Angstrom]")
        plt.savefig("barrierTest.png")
     

class GamowCalculator:
    def __init__(self, field:float = 5., radius:float = 1000, gamma:float = 10., XCdataFile:str = "tabulated/XCdata_W110.npy") -> None:
        self.field = field
        self.radius = radius
        self.gamma = gamma
        self.readXCpotential(XCdataFile)

    def readXCpotential(self, XCTabulationFile:str = "") -> None:
        try:
            self.XCdata = np.load(XCTabulationFile, allow_pickle=True).item()
            self.zMinimum = self.XCdata["extensionStartPoint"] + 1.e-5
        except(IOError):
            print("Warning: XC data file not found. Reverting to simple image XC")
            self.XCdata = None
            self.zMinimum = 1.e-2

    def imagePotential(self, z:np.ndarray):
        return  Globals.imageChargeConstant * ( 1. / z - .5 / self.radius) 

    def electrostaticPotential(self, z:np.ndarray):
        return self.field * (self.radius * (self.gamma - 1.) * z + z**2) / (self.gamma * z + self.radius * (self.gamma - 1.))
    
    def XCpotential(self, z:np.ndarray):

        convertToFloat = False
        if not isinstance(z, np.ndarray):
            z = np.array([z])
            convertToFloat = True

        imagePotential = self.imagePotential(z)
        if (self.XCdata is None):
            outPut = imagePotential
        else:
            imagePotential[z <= 0] = 20.
            imagePotential[imagePotential > 20.] = 20.
            potentialPBE = np.polyval(self.XCdata["polynomial"], z)
            potentialPBE[z > self.XCdata["polynomialRange"][1]] = 0.

            lowRange = np.where(z < self.XCdata["polynomialRange"][0])
            potentialPBE[lowRange] = self.XCdata["extensionPrefactor"] / (z[lowRange] - self.XCdata["extensionStartPoint"])
            transitionFunction = .5 * scipy.special.erfc((z - self.XCdata["transitionPoint"]) / self.XCdata["transitionWidth"])
            outPut = transitionFunction * potentialPBE + (1. - transitionFunction) * imagePotential

        if convertToFloat:
            z = z[0]
            outPut = outPut[0]
        return outPut

    def plotBarrier(self):
        zPlot = np.linspace(self.zMinimum + .05, 3., 512)
        plt.plot(zPlot, self.barrierFunction(zPlot, 0.))
        plt.grid()
        plt.xlabel("z [nm]")
        plt.ylabel("barrier [eV]")
        plt.savefig("barrier.png")
        
    def findBarrierMax(self) -> float:
        optimization =  opt.minimize(self.negativeBarrierFunction,.5, bounds=[(self.zMinimum, 3.)])
        self.barrierMaximumLocation = optimization.x
        self.minBarrierDepth = self.negativeBarrierFunction(self.barrierMaximumLocation)[0] + 1.e-3
        return self.barrierMaximumLocation

    def barrierFunction(self, z, barrierDepth):
        return barrierDepth - self.electrostaticPotential(z) - self.XCpotential(z)
    
    def negativeBarrierFunction(self, z):
        return self.electrostaticPotential(z) + self.XCpotential(z)

    def findZeros(self, barrierDepth:float) -> None:
        rootResult = opt.root_scalar(self.barrierFunction, args=barrierDepth, bracket=[self.zMinimum, self.barrierMaximumLocation], xtol=1.e-8)
        self.leftRoot = rootResult.root
        rootResult = opt.root_scalar(self.barrierFunction, args=barrierDepth, bracket=[self.barrierMaximumLocation, 10.])
        self.rightRoot = rootResult.root

    def integrantWKB(self, z:float, barrierDepth:float) -> float:
        barrier = np.maximum(self.barrierFunction(z, barrierDepth), 0.)
        return barrier ** .5

    def calculateGamow(self, barrierDepth:float) -> float:
        self.findZeros(barrierDepth)
        try:
            return Globals.gamowPrefactor * ig.quad(self.integrantWKB, self.leftRoot, self.rightRoot, args=barrierDepth, epsabs=1.e-5, epsrel=1.e-5)[0]
        except(ig.IntegrationWarning):
            print("Warning: quad has a problem. Reverting to trapz")
            zTrapz = np.linspace(self.leftRoot, self.rightRoot, 256)
            return Globals.gamowPrefactor * np.trapz(self.integrantWKB(zTrapz, barrierDepth), zTrapz)

    def findMaxBarrierDepth(self, maximumGamow:float = 20.) -> float:
        self.maxBarrierDepth = self.minBarrierDepth
        for i in range(100):
            self.maxBarrierDepth *= 1.5
            Gnew = self.calculateGamow(self.maxBarrierDepth)
            if (Gnew > maximumGamow):
                break
        
        GminusGmax = lambda barrierDepth: self.calculateGamow(barrierDepth) - maximumGamow
        self.maxBarrierDepth = opt.root_scalar(GminusGmax, bracket=[self.maxBarrierDepth/1.5, self.maxBarrierDepth]).root
        return self.maxBarrierDepth

    def calculateGamowCurveWKB(self, numberOfPoints:int = 32) -> None:
        self.findBarrierMax()
        self.findMaxBarrierDepth()

        self.barrierDepthVector = np.linspace(self.minBarrierDepth, self.maxBarrierDepth, numberOfPoints)
        self.gamowVector = np.copy(self.barrierDepthVector)
        for i in range(numberOfPoints):
            self.gamowVector[i] = self.calculateGamow(self.barrierDepthVector[i])

        self.transmissionCoefficients = Globals.transmissionCoefficientForGamow(self.gamowVector)

    def calculateGamowCurveFDM(self, numberOfPoints:int = 32, potentialDepth:float = 10., Npoints:int = 510, maxDepth = 6., minDepth=-1.) -> None:
        self.findZeros(potentialDepth)

        z = np.linspace(self.leftRoot, self.rightRoot, Npoints)
        potential = self.barrierFunction(z, barrierDepth=0.)
        solver = Schrodinger1DSolverFDM(potential, dx=10. * (z[1] - z[0]), energy=0.)
        energies, self.transmissionCoefficients = solver.calculateTransmissionForEnergyRange(-maxDepth, -minDepth, numberOfPoints)
        self.barrierDepthVector = -energies

    def findMinimumZ(self, minimumPotential):
        func = lambda x : self.XCpotential(x) - minimumPotential
        self.zMinimum = opt.fsolve(func, 0.1)
        return self.zMinimum

class Schrodinger1DSolverIVP(GamowCalculator):
    def __init__(self, field: float = 5, radius: float = 1000, gamma: float = 10, XCdataFile: str = "tabulated/XCdata_W110.npy", barrierDepth=0.) -> None:
        super().__init__(field, radius, gamma, XCdataFile)
        self.barrierDepth = barrierDepth
        self.kConstant = 100./3.80998212 # 2m/hbar in Å^-2 eV^-1
        self.maxLength = 3.
        # self.zMinimum = 0.0001

    def differentialSystem(self, x, y):
        return [y[1], y[0] *  self.kConstant * (self.barrierFunction(x, self.barrierDepth))]
    
    def solveSystem(self):
        self.solution = ig.solve_ivp(self.differentialSystem, [self.maxLength, self.zMinimum], \
                                     [1., np.sqrt(0 * 1j + self.kConstant * self.barrierFunction(self.maxLength, self.barrierDepth))], rtol=1.e-6)#, jac=self.jacobian)

    def jacobian(self, x, y):
        return [[0., 1.], [self.kConstant * (self.barrierFunction(x, self.barrierDepth)), 0.]]
    
    def calculateTransmission(self):
        leftWaveNumber = np.sqrt(0 * 1j - self.kConstant * self.barrierFunction(self.maxLength, self.barrierDepth))
        matrix = np.array([[1., 1.], [1j * leftWaveNumber, - 1j * leftWaveNumber ]])
        sol = np.abs(linalg.solve(matrix, self.solution.y[:,-1]))
        transmission = (1 / sol[0])**2
        reflection = (sol[1] / sol[0])**2
        assert abs(reflection + transmission - 1) < 1.e-4, "reflection + transmission -1 = %g"%(reflection + transmission - 1.)
        return transmission
        
    
    def plotWaveFunction(self):
        plt.figure()
        fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
        z = np.linspace(self.zMinimum + 0.01, self.maxLength, 256)
        ax1.plot(z, self.barrierFunction(z, self.barrierDepth), label=r"U(z) - E")
        ax1.grid()
        ax1.legend()
        ax1.set_ylabel(r"barrier [eV]")

        ax2.plot(self.solution.t, np.real(self.solution.y[0]), label=r"$\Re(\Psi)$", color = colors[2])
        ax2.plot(self.solution.t, np.imag(self.solution.y[0]), label=r"$\Im(\Psi)$", color = colors[3])
        # ax2r = ax2.twinx()
        # ax2r.plot(z, np.angle(self.solution[1:-1]), label=r"Phase", color=colors[1])
        ax2.legend()
        ax2.grid()
        ax2.set_ylabel(r"$\Psi$")
        # ax2r.set_ylabel("Angle ($\Psi$)")
        # ax2r.legend()
        ax2.set_xlabel(r"z [$\textrm{\AA}$]")
        plt.savefig("barrierTest.png")

    def printFinalValue(self):
        print(self.solution.y[0])

    def calculateGamow(self, barrierDepth):
        self.barrierDepth = barrierDepth
        self.maxLength = self.lengthEstimation()
        self.solveSystem()
        transmissionCoefficient = self.calculateTransmission()
        return Globals.GamowForTransmissionCoefficient(transmissionCoefficient), transmissionCoefficient

    def lengthEstimation(self):
        requiredBarrierDepth = self.barrierDepth + self.field**(2./3) #how much the depth should be to satisfy JWKB approximation 2m(E-V) > (10 m hbar F)^2/3
        factor = self.radius * (self.gamma - 1.)
        poly = np.array([1., factor - requiredBarrierDepth * self.gamma / self.field, - requiredBarrierDepth * factor / self.field])
        roots = np.roots(poly)
        return max(roots)

    def calculateGamowForDepths(self, barrierDepths):
        gamowVector = np.copy(barrierDepths)
        for i in range(len(barrierDepths)):
            gamowVector[i] = self.calculateGamow(barrierDepths[i])[0]
        return gamowVector



    def calculateGamowCurve(self, numberOfPoints:int = 128) -> None:
        self.findBarrierMax()
        self.findMaxBarrierDepth()

        self.barrierDepthVector = np.linspace(self.minBarrierDepth, self.maxBarrierDepth, numberOfPoints)
        self.gamowVector = np.copy(self.barrierDepthVector)
        for i in range(numberOfPoints):
            self.gamowVector[i], self.transmissionCoefficients[i] = self.calculateGamow(self.barrierDepthVector[i])


         

class GamowTabulator(GamowCalculator):
    def __init__(self,  Nf = 256, Nr = 128, Ngamma = 1, Npoly = 4, XCdataFile = "") -> None:
        super().__init__(XCdataFile=XCdataFile)

        self.fieldRange = [0.5, 12.]
        self.minRadius = 0.5
        self.maxGamma = 1.e3

        #in this case the reasonable single Gamma value is 10
        if(Ngamma == 1):
            self.maxGamma = 10.

        self.numberOfFieldValues = Nf
        self.numberOfRadiusValues = Nr
        self.numberOfGammaValues = Ngamma
        self.numberOfPolynomialTerms = Npoly

        self.inverseFieldValues = np.linspace(1 / self.fieldRange[1], 1. / self.fieldRange[0], self.numberOfFieldValues)
        self.inverseRadiusValues = np.linspace(1.e-4, 1 / self.minRadius, self.numberOfRadiusValues)
        self.inverseGammaValues = np.linspace(1. / self.maxGamma, 1., self.numberOfGammaValues)
    
    def tabulateGamowTable(self, outputFolder = None):
        """Looks for the files where the precaculate barriers are stored. Then it uses interpolation methods to make the most accurate barrier for the given 
        input (electric field, tip radius and gamma exponent). Gtab is stores the polinomial that gives its shape to the barrier.
        """

        if outputFolder is None:
            outputFolder = Globals.tabulationPath

        self.gamowTable = np.ones([self.numberOfGammaValues, self.numberOfRadiusValues, self.numberOfFieldValues, self.numberOfPolynomialTerms + 2])

        for i in range(self.numberOfGammaValues):
            for j in range(self.numberOfRadiusValues):
                for k in range(self.numberOfFieldValues):
                    
                    self.field = 1 / self.inverseFieldValues[k]
                    self.radius = 1 / self.inverseRadiusValues[j]
                    self.gamma = 1 / self.inverseGammaValues[i]
                    
                    self.calculateGamowCurveWKB()
                    
                    try:
                        fittedPolynomialCoefficients = np.polyfit(self.barrierDepthVector, self.gamowVector, self.numberOfPolynomialTerms - 1)
                    except(np.RankWarning):
                        print("Rank Warning for F = %g, R = %g, gamma = %g"%(self.field, self.radius, self.gamma))
                        plt.plot(self.barrierDepthVector, self.gamowVector)
                        plt.show()
                    self.gamowTable[i, j, k, :] = np.append(fittedPolynomialCoefficients, [self.minBarrierDepth, self.maxBarrierDepth])    

        if (self.numberOfGammaValues == 1):
            self.gamowTable = np.reshape(self.gamowTable, (self.numberOfRadiusValues, self.numberOfFieldValues, self.numberOfPolynomialTerms + 2))
            if (self.numberOfRadiusValues == 1):
                self.gamowTable = np.reshape(self.gamowTable, (self.numberOfFieldValues, self.numberOfPolynomialTerms + 2))

        os.system("mkdir -p %s"%outputFolder)
        np.save(outputFolder + "/GamowTable.npy", self.gamowTable)

        np.save(outputFolder + "/tabLimits.npy", np.array([min(self.inverseFieldValues), max(self.inverseFieldValues), \
                min(self.inverseRadiusValues), max(self.inverseRadiusValues), min(self.inverseGammaValues), max(self.inverseGammaValues)]))
        


# calculator = GamowCalculator(XCdataFile="")

# solver = Schrodinger1DSolverIVP(barrierDepth=4.5, XCdataFile="")


calculator = GamowCalculator()

solver = Schrodinger1DSolverIVP()

# solver.solveSystem()

# solver.plotWaveFunction()
# solver.printFinalValue()
# solver.calculateTransmission()

# print(solver.barrierFunction(solver.lengthEstimation(), solver.barrierDepth))

# calculator.plotBarrier()
solver.findMinimumZ(7.)
# calculator.findMinimumZ(7.)

import time
start_time = time.time()
calculator.calculateGamowCurveWKB(32)
print("calculating WKB: %g seconds ---" %(time.time() - start_time))


gamowDerivative = (calculator.gamowVector[1] - calculator.gamowVector[0]) / (calculator.barrierDepthVector[1] - calculator.barrierDepthVector[0])
extraBarrierDepth = np.linspace(0., calculator.barrierDepthVector[0], 8)
extraGamow = gamowDerivative * (extraBarrierDepth - calculator.barrierDepthVector[0])
allGamow = np.concatenate((extraGamow, calculator.gamowVector))
# allTransmission = Globals.transmissionCoefficientForGamow(allGamow)
allDepths = np.concatenate((extraBarrierDepth, calculator.barrierDepthVector))

filterPoints = np.where(allGamow < 100.)

plt.figure()
plt.plot(allDepths[filterPoints], allGamow[filterPoints], ".-", label="WKB")


start_time = time.time()
D = solver.calculateGamowForDepths(allDepths[filterPoints])
print("calculating with IVP %g seconds ---" % (time.time() - start_time))


plt.plot(allDepths[filterPoints], D, ".-", label="RK")
plt.grid()
plt.legend()
plt.xlabel("barrier Depth [eV]")
plt.ylabel("Gamow factor")
plt.savefig("TransmissionComparison.png")
pltmod.pushFigureToServer(plt.gcf())
# plt.show()
