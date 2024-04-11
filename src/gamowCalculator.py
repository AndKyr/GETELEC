import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os
from getelec import Globals
import scipy.optimize as opt
import scipy.integrate as ig
import scipy.special
import scipy.linalg as linalg
import scipy.interpolate as interp
from collections.abc import Callable



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


class Schrodinger1DSolver:
    def __init__(self, potentialVector:np.ndarray = None, energy:float = 0., potentialFunction:Callable = None, xLimits:np.ndarray = np.array([0,1])):

        if potentialVector is None and potentialFunction is None:
            raise ValueError("The solver must be initialized with either a potentialFunction or a potentialVector")
        
        if potentialFunction is not None and xLimits is None:
            raise ValueError("If you initialize with potentialFunction, you have to give the limits of evaluation of the function")
        
        self.hbarSqrOver2m = 3.80998212e-2 # nm^2 (square Ångströms) (eV (electronvolt))
        self.potentialVector = potentialVector
        self.xLimits = np.copy(xLimits)
        self.potentialFunction = potentialFunction
        self.energy = energy
        
    def getPotentialFunctionForVector(self):
        if (self.potentialFunction is None):
            self.potentialFunction = interp.interp1d(self.xVector, self.potentialVector)

    def getPotentialVector(self, Npoints:int = 254, updatexVector:bool = True):
        if (self.potentialFunction is None):
            return
        
        if (updatexVector):
            self.xVector = np.linspace(self.xLimits[0], self.xLimits[1], Npoints)
            self.dx = self.xVector[1] - self.xVector[0]

        self.potentialVector = self.potentialFunction(self.xVector)

    def calculateWaveVectors(self):
        if (self.potentialVector is not None and len(self.potentialVector) > 0):
            leftEdgePotential = self.potentialVector[0]
            rightEdgePotential = self.potentialVector[-1]
        else:
            leftEdgePotential = self.potentialFunction(self.xLimits[0])
            rightEdgePotential = self.potentialFunction(self.xLimits[1])
        
        assert max(leftEdgePotential, rightEdgePotential) < self.energy, "the energy has to be more than the potential on the left and right edges"

        self.kLeft = np.sqrt((self.energy - leftEdgePotential) / self.hbarSqrOver2m)
        self.kRight = np.sqrt((self.energy - rightEdgePotential) / self.hbarSqrOver2m)

    def getSolutionVector(self):
        self.getPotentialVector()
        self.solutionVector = np.copy(self.potentialVector)

    def calculateTransmissionForEnergies(self, energies:np.ndarray)->np.ndarray:
        raise NotImplementedError("function not implemented in the base class")


    def calculateTransmissionForEnergy(self, energy:float):
        raise NotImplementedError("function not implemented in the base class")
    
    def plotWaveFunction(self):
        self.getSolutionVector()
        plt.figure()
        fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
        ax1.plot(self.xVector, self.potentialVector, label=r"U(z)")
        ax1.plot([self.xVector[0], self.xVector[-1]], [self.energy, self.energy], label=r"E")
        ax1.grid()
        ax1.legend()
        ax1.set_ylabel(r"barrier [eV]")

        ax2.plot(self.xVector, np.real(self.solutionVector), label=r"$\Re(\Psi)$", color = colors[2])
        ax2.plot(self.xVector, np.imag(self.solutionVector), label=r"$\Im(\Psi)$", color = colors[3])
        # ax2r = ax2.twinx()
        # ax2r.plot(z, np.angle(self.solution[1:-1]), label=r"Phase", color=colors[1])
        ax2.legend()
        ax2.grid()
        ax2.set_ylabel(r"$\Psi$")
        # ax2r.set_ylabel("Angle ($\Psi$)")
        # ax2r.legend()
        ax2.set_xlabel(r"z [$\textrm{\AA}$]")
        plt.savefig("barrierAndWaveFunctions.png")

            


class Schrodinger1DSolverFDM(Schrodinger1DSolver):
    def __init__(self, potentialVector:np.ndarray = None, energy:float = 0., potentialFunction:Callable = None, xLimits:np.ndarray = None, dx = 0.1, Npoints:int = 254):
        super().__init__(potentialVector=potentialVector, energy=energy, potentialFunction=potentialFunction, xLimits=xLimits)
        self.dx = dx
        self.getPotentialVector(Npoints=Npoints)
        
        self.Nx = Npoints
        self.MConstant = self.dx**2 / self.hbarSqrOver2m # Å^2 (square Ångströms) (eV (electronvolt))
        self.calculateWaveVectors()

    def setPotentialAndEnergy(self, potentialVector:np.ndarray, energy:float):
        self.potentialVector = potentialVector
        self.energy = energy

    def setSparseMatrixSystem(self):
        #sets the sparse matrix for solving the implicit FDM problem
        mainDiagonal = np.ones(self.Nx + 2, dtype=complex)
        mainDiagonal[2:-2] = 2. + self.MConstant * (self.potentialVector[1:-1] - self.energy) # -ui+1 - ui-1 + 2 ui + dx^2 (2m/hbar^2) (Vi - E) ui = 0
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


    def calculateTransmissionForEnergy(self, energy:float):
        self.energy = energy
        self.solveFirstTime()
        return self.getTransmissionProbability()
        
    def calculateTransmissionForEnergies(self, energies:np.ndarray)->np.ndarray:
        transmissionCoefficients = np.copy(energies)
        self.solveFirstTime()
        for i in range(len(energies)):
            self.solveWithNewEnergy(energies[i])
            transmissionCoefficients[i] = self.getTransmissionProbability()
        return transmissionCoefficients
    
    def getSolutionVector(self):
        super().getSolutionVector()
        assert len(self.solution) == len(self.potentialVector) + 2, "the solution vector length doesn't match with the potential vector"
        self.solutionVector = self.solution[1:-1]
    

    def plotBarrier(self):
        plt.figure(figsize=figureSize)
        fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
        z = np.linspace(0., (self.Nx - 1) * self.dx, self.Nx)
        ax1.plot(z, self.potentialVector, label="U")
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
     
class Schrodinger1DSolverIVP(Schrodinger1DSolver):
    def __init__(self, energy:float = 0., potentialFunction:Callable = None, xLimits:np.ndarray = None) -> None:
        super().__init__(potentialFunction=potentialFunction, xLimits=xLimits, energy=energy)
        self.kConstant = 1/self.hbarSqrOver2m # 2m/hbar in nm^-2 eV^-1
        self.calculateWaveVectors()

    def setEnergy(self, energy:float):
        self.energy = energy

    def differentialSystem(self, x:np.ndarray, y:np.ndarray) -> np.ndarray:
        return [y[1], y[0] *  self.kConstant * (self.potentialFunction(x) - self.energy)]
    
    def solveSystem(self) -> None:
        self.solution = ig.solve_ivp(self.differentialSystem, [self.xLimits[1], self.xLimits[0]], [1., 1j * self.kRight], rtol=1.e-8, vectorized=True)#, jac=self.jacobian)

    def jacobian(self, x, y):
        return [[0., 1.], [self.kConstant * (self.potentialFunction(x) - self.energy), 0.]]
    
    def getTransmissionProbability(self) -> float:
        matrix = np.array([[1., 1.], [1j * self.kRight, - 1j * self.kRight ]])
        sol = np.abs(linalg.solve(matrix, self.solution.y[:,-1]))
        transmission = (1 / sol[0])**2
        reflection = (sol[1] / sol[0])**2
        assert abs(reflection + transmission - 1) < 1.e-4, "reflection + transmission -1 = %g"%(reflection + transmission - 1.)
        return transmission
    
    def getSolutionVector(self):
        super().getSolutionVector()
        self.xVector = self.solution.t
        self.solutionVector = self.solution.y[0]


    def calculateTransmissionForEnergy(self, energy:float):
        self.energy = energy
        self.solveSystem()
        return self.getTransmissionProbability()

    def calculateTransmissionForEnergies(self, energies:np.ndarray)->np.ndarray:
        transmissionCoefficients = np.copy(energies)
        for i in range(len(energies)):
            self.energy = energies[i]
            self.solveSystem()
            transmissionCoefficients[i] = self.getTransmissionProbability()
        return transmissionCoefficients
    


class GamowCalculator:
    def __init__(self, field:float = 5., radius:float = 1000, gamma:float = 10., XCdataFile:str = "tabulated/XCdata_W110.npy", solverType:str = "WKB", minimumPotential:float = 7.) -> None:
        self.field = field
        self.radius = radius
        self.gamma = gamma
        self.readXCpotential(XCdataFile)
        self.potentialFunction = lambda x: - self.electrostaticPotential(x) - self.XCpotential(x)
        self.findMinimumZ(minimumPotential)
        self.maxLength = 10.

        if (solverType == "FDM"):
            self.solver = Schrodinger1DSolverFDM(potentialFunction=self.potentialFunction, xLimits=np.append(self.zMinimum, self.maxLength))
        elif(solverType == "IVP"):
            self.solver = Schrodinger1DSolverIVP(potentialFunction=self.potentialFunction, xLimits=np.append(self.zMinimum, self.maxLength))
        else:
            self.solver = None


    def setSolver(self, solver:Schrodinger1DSolver):
        self.solver = solver
        self.solver.potentialFunction = lambda x: - self.electrostaticPotential(x) - self.XCpotential(x)

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

    def barrierFunction(self, z:np.ndarray, barrierDepth:float) -> np.ndarray:
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

    def calculateGamowWKB(self, barrierDepth:float) -> float:
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
            if (self.barrierFunction(self.zMinimum, self.maxBarrierDepth * 1.5) > 0):
                return
            self.maxBarrierDepth *= 1.5
            Gnew = self.calculateGamowWKB(self.maxBarrierDepth)
            if (Gnew > maximumGamow):
                break
        
        GminusGmax = lambda barrierDepth: self.calculateGamowWKB(barrierDepth) - maximumGamow
        self.maxBarrierDepth = opt.root_scalar(GminusGmax, bracket=[self.maxBarrierDepth/1.5, self.maxBarrierDepth]).root
        return self.maxBarrierDepth

    def calculateGamowCurve(self, numberOfPoints:int = 32) -> None:
        self.findBarrierMax()
        self.findMaxBarrierDepth()

        self.barrierDepthVector = np.linspace(self.minBarrierDepth, self.maxBarrierDepth, numberOfPoints)
        self.gamowVector = np.copy(self.barrierDepthVector)
        if self.solver is None:  
            for i in range(numberOfPoints):
                self.gamowVector[i] = self.calculateGamowWKB(self.barrierDepthVector[i])
            self.transmissionCoefficients = Globals.transmissionCoefficientForGamow(self.gamowVector)
        else:
            self.maxLength = self.lengthEstimation(self.maxBarrierDepth)
            self.solver.xLimits = np.append(self.zMinimum, self.maxLength)
            self.transmissionCoefficients = self.solver.calculateTransmissionForEnergies(-self.barrierDepthVector)
            self.gamowVector = Globals.GamowForTransmissionCoefficient(self.transmissionCoefficients)


    def findMinimumZ(self, minimumPotential):
        func = lambda x : self.XCpotential(x) - minimumPotential
        self.zMinimum = opt.fsolve(func, 0.1)
        return self.zMinimum

    def setSolverpotential(self) -> None:
        if self.solver is None:
            return
        if (type(self.solver) is Schrodinger1DSolverFDM):
            self.solver.getPotentialVector(updatexVector=False)


    def calculateGamow(self, barrierDepth):
        self.maxLength = self.lengthEstimation()
        self.solver.xLimits = np.array([self.zMinimum, self.maxLength])

        if self.solver is not None:
            transmissionCoefficient = self.solver.calculateTransmissionForEnergy(-barrierDepth)
            return Globals.GamowForTransmissionCoefficient(transmissionCoefficient), transmissionCoefficient
        else:
            return self.calculateGamowWKB(barrierDepth)

    def lengthEstimation(self, barrierDepth):
        requiredBarrierDepth = barrierDepth + self.field**(2./3) #how much the depth should be to satisfy JWKB approximation 2m(E-V) > (10 m hbar F)^2/3
        factor = self.radius * (self.gamma - 1.)
        poly = np.array([1., factor - requiredBarrierDepth * self.gamma / self.field, - requiredBarrierDepth * factor / self.field])
        roots = np.roots(poly)
        return max(roots)


calculator = GamowCalculator()



import time
start_time = time.time()
calculator.calculateGamowCurve(32)
print("calculating WKB: %g seconds ---" %(time.time() - start_time))

gamowDerivative = (calculator.gamowVector[1] - calculator.gamowVector[0]) / (calculator.barrierDepthVector[1] - calculator.barrierDepthVector[0])
extraBarrierDepth = np.linspace(0., calculator.barrierDepthVector[0], 8)
extraGamow = gamowDerivative * (extraBarrierDepth - calculator.barrierDepthVector[0])
allGamow = np.concatenate((extraGamow, calculator.gamowVector))
# allTransmission = Globals.transmissionCoefficientForGamow(allGamow)
allDepths = np.concatenate((extraBarrierDepth, calculator.barrierDepthVector))
filterPoints = np.where(allGamow < 100.)

plt.plot(allDepths[filterPoints], allGamow[filterPoints], ".-", label="WKB")


calculator = GamowCalculator(solverType="FDM")

start_time = time.time()
calculator.calculateGamowCurve(32)

# calculator.solver.plotWaveFunction()
print("calculating FDM: %g seconds ---" %(time.time() - start_time))
print(calculator.solver.dx)
plt.plot(calculator.barrierDepthVector, calculator.gamowVector, ".-", label="FDM")

calculator = GamowCalculator(solverType="IVP")
start_time = time.time()
calculator.calculateGamowCurve(32)
print("calculating with IVP %g seconds ---" % (time.time() - start_time))
plt.plot(calculator.barrierDepthVector, calculator.gamowVector, ".-", label="IVP")


plt.grid()
plt.legend()
plt.xlabel("barrier Depth [eV]")
plt.ylabel("Gamow factor")
plt.savefig("TransmissionComparison.png")
# pltmod.pushFigureToServer(plt.gcf())
plt.show()
