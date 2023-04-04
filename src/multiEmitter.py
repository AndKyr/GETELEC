
import getelec as gt
import numpy as np

import matplotlib.pyplot as plt



def sampleVariable(distributionType:str = "normal", average:float= 0., std:float = 1., numberOfSamples:int = 32):
    if (distributionType == 'lognormal'):
        mu = np.log(average**2/np.sqrt(std**2 + average**2))
        sigma = np.sqrt(np.log(1+std**2 / average**2))
        return np.random.lognormal(mu, sigma, numberOfSamples)
    elif(distributionType == 'normal'):
        return np.random.normal(average, std, numberOfSamples)
    elif(distributionType == 'uniform'):
        low = average - np.sqrt(12) * std * 0.5
        high = average + np.sqrt(12) * std * 0.5
        return np.random.uniform(low, high, numberOfSamples)
    elif(distributionType == 'exponential'):
        return np.random.exponential(average, numberOfSamples)
    else:
        print ("wrong distribution for heights")
class MultiEmitter():
    """This class simulates the situation of many random emitters and their collective I-V curve"""
    
    def __init__(self, emitterModel = gt.GETELECModel(), numberOfEmitters = 128, averageHeight = 300, stdHeight = 30, averageRadius = 15, stdRadius = 1, minBetaValue = 1., maxBetaValue = 1000.):
        self.emitterModel = emitterModel
        self.numberOfEmitters = numberOfEmitters
        self.averageHeight = averageHeight
        self.stdHeight = stdHeight
        self.averageRadius = averageRadius
        self.stdRadius = stdRadius
        self.minBetaValue = minBetaValue
        self.maxBetaValue = maxBetaValue
        self.averageBeta = self.averageHeight / self.averageRadius
        self.stdBeta = np.sqrt(self.stdHeight**2 + (self.averageHeight * self.stdRadius / self.averageRadius**2)**2)

        self.currents = np.zeros(numberOfEmitters)
        self.betas = np.zeros(numberOfEmitters)
        self.radii = np.zeros(numberOfEmitters)

        
    def sampleEmitters(self, mode = "radiusAndHeight", heightDistribution = 'lognormal', radiiDistribution = 'lognormal', betaDistribution = "uniform"):
        
        if (mode == "radiusAndHeight"):
            self.heights = sampleVariable(distributionType=heightDistribution, average=self.averageHeight, std=self.stdHeight, numberOfSamples = self.numberOfEmitters)
            self.radii = sampleVariable(distributionType=radiiDistribution, average=self.averageRadius, std=self.stdRadius, numberOfSamples= self.numberOfEmitters) 
            self.betas = self.heights / self.radii
        elif(mode == "betaAndRadius"):
            self.betas = sampleVariable(betaDistribution, self.averageBeta, self.stdBeta, self.numberOfEmitters)
            self.radii = sampleVariable(distributionType=radiiDistribution, average=self.averageRadius, std=self.stdRadius, numberOfSamples= self.numberOfEmitters)
            self.heights = self.betas * self.radii
        else:
            assert False, "wrong sampling mode"
            
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.pi * self.radii**2

        self.betas[self.betas > self.maxBetaValue] = self.maxBetaValue
        self.betas[self.betas < self.minBetaValue] = self.minBetaValue
        
           
    def emit(self, Ffar):
        self.emitterModel.setParameters(field = Ffar * self.betas, radius = self.radii)
        self.emitterModel.calculateCurrentDensity()

        self.currents = np.array(self.emitterModel.getCurrentDensity()) * np.pi * self.radii**2            
        return sum(self.currents)   
        
        
    def plotHistograms(self):
        
        plt.figure()
        plt.hist(self.radii, 100)
        plt.xlabel("Radii")
        plt.ylabel("count")
        
        plt.figure()
        plt.hist(self.heights, 100)
        plt.xlabel("Heights")
        plt.ylabel("count")
        
        plt.figure()
        plt.hist(self.betas, 100)
        plt.xlabel(r"$\beta$")
        plt.ylabel("count")
        
        
        plt.figure()

    
    def print_data_out(self):
                   
        print ("heights = %.3g +- %.3g"%(np.mean(self.heights), np.std(self.heights)))
        print ("radii = %.3g +- %.3g"%(np.mean(self.radii), np.std(self.radii)))
        print ("betas = %.3g +- %.3g"%(np.mean(self.betas),np.std(self.betas)))
        print ("maxbeta = %g"%(np.max(self.betas)))
        

        
if (__name__ == "__main__"): #some testing operations
    farFields = np.linspace(0.035, 0.08, 128)
    totalCurrents = np.copy(farFields)
    for i in range(8):
        getelecModel = gt.GETELECModel(emitterType="metal", workFunction=4.5, temperature=300, gamma=10., effectiveMassConduction=1.)
        me = MultiEmitter(emitterModel= getelecModel, averageHeight=500, stdHeight=200., averageRadius=12., stdRadius=4., numberOfEmitters=256)
        me.sampleEmitters()

        me.print_data_out()
        
        for i in range(len(farFields)):
            totalCurrents[i] = me.emit(Ffar=farFields[i])

        
        plt.semilogy(1/farFields, totalCurrents)
    plt.savefig("multiemitterTest.png")




