
import getelec as gt
import numpy as np

import matplotlib.pyplot as plt

class MultiEmitter():
    
    def __init__(self, emitterModel = gt.GETELECModel(), numberOfEmitters = 128, averageHeight = 300, stdHeight = 30, averageRadius = 15, stdRadius = 1):
        self.emitterModel = emitterModel
        self.numberOfEmitters = numberOfEmitters
        self.averageHeight = averageHeight
        self.stdHeight = stdHeight
        self.averageRadius = averageRadius
        self.stdRadius = stdRadius

        self.currents = np.zeros(numberOfEmitters)
        self.betas = np.zeros(numberOfEmitters)
        self.radii = np.zeros(numberOfEmitters)
        
       
        
    def getEmitters(self, maxBetaValue = 1.e50, minBetaValue = 1., heightDistribution = 'lognormal', radiiDistribution = 'lognormal', appendMaxValues = True):
        
        if (heightDistribution == 'lognormal'):
            mu = np.log(self.averageHeight**2/np.sqrt(self.stdHeight**2 + self.averageHeight**2))
            sigma = np.sqrt(np.log(1+self.stdHeight**2/self.averageHeight**2))
            self.heights = np.random.lognormal(mu, sigma, self.numberOfEmitters)
        elif(heightDistribution == 'normal'):
            mu = self.averageHeight
            sigma = self.stdHeight
            self.heights = np.random.normal(mu, sigma, self.numberOfEmitters)
        elif(heightDistribution == 'uniform'):
            low = self.averageHeight - np.sqrt(12) * self.stdHeight * 0.5
            high = self.averageHeight + np.sqrt(12) * self.stdHeight * 0.5
            self.heights = np.random.uniform(low, high, self.numberOfEmitters)
        elif(heightDistribution == 'exponential'):
            self.heights = np.random.exponential(self.averageHeight, self.numberOfEmitters)
        else:
            print ("wrong distribution for heights")
        
        if (radiiDistribution == "lognormal"):
            mu = np.log(self.averageRadius**2/np.sqrt(self.stdRadius**2 + self.averageRadius**2))
            sigma = np.sqrt(np.log(1+self.stdRadius**2/self.averageRadius**2))
            self.radii = np.random.lognormal(mu, sigma, self.numberOfEmitters)
        elif(radiiDistribution == 'normal'):
            mu = self.averageRadius
            sigma = self.stdRadius
            self.radii = np.random.normal(mu, sigma, self.numberOfEmitters)
        elif(radiiDistribution == 'uniform'):
            low = self.averageRadius - np.sqrt(12) * self.stdRadius * 0.5
            high = self.averageRadius + np.sqrt(12) * self.stdRadius * 0.5
            self.radii = np.random.uniform(low, high, self.numberOfEmitters)
        elif(radiiDistribution == 'exponential'):
            curvatures = np.random.exponential(1./self.averageRadius, self.numberOfEmitters)
            self.radii = 1./curvatures
        else:
            print ("wrong distribution for radii"    )
            
        self.betas = self.heights / self.radii
        
        self.betas[self.betas > maxBetaValue] = maxBetaValue
        self.betas[self.betas < minBetaValue] = minBetaValue
        
        
    def getEmittersFromBetaDistribution(self, min_beta = 50, max_beta = 100, beta_dist = 'uniform', r_dist = 'lognormal'):

        if(beta_dist == 'normal'):
            mu = 0.5 * (min_beta + max_beta)
            sigma = (max_beta - min_beta) * 0.25
            self.betas = np.random.normal(mu, sigma, self.numberOfEmitters)
        elif(beta_dist == 'uniform'):
            low = self.averageHeight - np.sqrt(12) * self.stdHeight * 0.5
            high = self.averageHeight + np.sqrt(12) * self.stdHeight * 0.5
            self.betas = np.random.uniform(min_beta, max_beta, self.numberOfEmitters)
        else:
            print ("wrong distribution for betas")
            
        if (r_dist == "lognormal"):
            mu = np.log(self.averageRadius**2/np.sqrt(self.stdRadius**2 + self.averageRadius**2))
            sigma = np.sqrt(np.log(1+self.stdRadius**2/self.averageRadius**2))
            self.radii = np.random.lognormal(mu, sigma, self.numberOfEmitters)
        elif(r_dist == 'normal'):
            mu = self.averageRadius
            sigma = self.stdRadius
            self.radii = np.random.normal(mu, sigma, self.numberOfEmitters)
        elif(r_dist == 'uniform'):
            low = self.averageRadius - np.sqrt(12) * self.stdRadius * 0.5
            high = self.averageRadius + np.sqrt(12) * self.stdRadius * 0.5
            self.radii = np.random.uniform(low, high, self.numberOfEmitters)
        else:
            print ("wrong distribution for radii")
            
        self.heights = self.betas * self.radii
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.pi * self.radii**2
        
           
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
        me.getEmitters(minBetaValue=10., maxBetaValue=300.)

        me.print_data_out()
        
        for i in range(len(farFields)):
            totalCurrents[i] = me.emit(Ffar=farFields[i])

        
        plt.semilogy(1/farFields, totalCurrents)
    plt.show()




