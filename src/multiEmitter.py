
import getelec_mod as gt
import numpy as np

import matplotlib.pyplot as plt

class MultiEmitter():
    
    def __init__(self, em = gt.emission_create(), N = 1, mu_h = 100., std_h = 10., \
                    mu_r = 10., std_r = 1.):
        self.emitter = em
        self.emitter.approx = 1
        self.N_emitters = N
        self.mu_height = mu_h
        self.sigma_height = std_h
        self.mu_radii = mu_r
        self.sigma_radii = std_r
        
        
    def set_emitters(self, radii, heights = np.array([]), betas = np.array([])):
        self.radii = radii
        
        if(not heights.any()):
            self.betas = betas
            self.heights = betas * radii
        else:
            self.betas = heights / radii
            self.heights = heights
        
        self.N_emitters = len(self.betas)
        
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.copy(self.currents)
        
        
        
    def getEmitters(self, maxBetaValue = 1.e50, minBetaValue = 1., heightDistribution = 'normal', radiiDistribution = 'lognormal', appendMaxValues = True):
        
        if (heightDistribution == 'lognormal'):
            mu = np.log(self.mu_height**2/np.sqrt(self.sigma_height**2 + self.mu_height**2))
            sigma = np.sqrt(np.log(1+self.sigma_height**2/self.mu_height**2))
            heights = np.random.lognormal(mu, sigma, self.N_emitters)
        elif(heightDistribution == 'normal'):
            mu = self.mu_height
            sigma = self.sigma_height
            heights = np.random.normal(mu, sigma, self.N_emitters)
        elif(heightDistribution == 'uniform'):
            low = self.mu_height - np.sqrt(12) * self.sigma_height * 0.5
            high = self.mu_height + np.sqrt(12) * self.sigma_height * 0.5
            heights = np.random.uniform(low, high, self.N_emitters)
        elif(heightDistribution == 'exponential'):
            heights = np.random.exponential(self.mu_height, self.N_emitters)
        else:
            print ("wrong distribution for heights")
        
        if (radiiDistribution == "lognormal"):
            mu = np.log(self.mu_radii**2/np.sqrt(self.sigma_radii**2 + self.mu_radii**2))
            sigma = np.sqrt(np.log(1+self.sigma_radii**2/self.mu_radii**2))
            radii = np.random.lognormal(mu, sigma, self.N_emitters)
        elif(radiiDistribution == 'normal'):
            mu = self.mu_radii
            sigma = self.sigma_radii
            radii = np.random.normal(mu, sigma, self.N_emitters)
        elif(radiiDistribution == 'uniform'):
            low = self.mu_radii - np.sqrt(12) * self.sigma_radii * 0.5
            high = self.mu_radii + np.sqrt(12) * self.sigma_radii * 0.5
            radii = np.random.uniform(low, high, self.N_emitters)
        elif(radiiDistribution == 'exponential'):
            curvatures = np.random.exponential(1./self.mu_radii, self.N_emitters)
            radii = 1./curvatures
        else:
            print ("wrong distribution for radii"    )
            
        betas = heights / radii
        
        isgood = np.where(np.logical_and(betas < maxBetaValue, betas > minBetaValue))
        
        
        self.radii = radii[isgood]
        self.heights = heights[isgood]
        self.betas = betas[isgood]
        self.N_emitters = len(self.betas)
        
        if appendMaxValues:
            self.betas = np.append(self.betas, maxBetaValue)
            self.radii = np.append(self.radii, min(self.radii))
            self.heights = np.append(self.heights, self.betas[-1] * self.radii[-1])
            self.N_emitters += 1
            
        
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.copy(self.currents)
        
        
    def getEmittersFromBetaDistribution(self, min_beta = 50, max_beta = 100, beta_dist = 'uniform', r_dist = 'lognormal'):

        if(beta_dist == 'normal'):
            mu = 0.5 * (min_beta + max_beta)
            sigma = (max_beta - min_beta) * 0.25
            self.betas = np.random.normal(mu, sigma, self.N_emitters)
        elif(beta_dist == 'uniform'):
            low = self.mu_height - np.sqrt(12) * self.sigma_height * 0.5
            high = self.mu_height + np.sqrt(12) * self.sigma_height * 0.5
            self.betas = np.random.uniform(min_beta, max_beta, self.N_emitters)
        else:
            print ("wrong distribution for betas")
            
        if (r_dist == "lognormal"):
            mu = np.log(self.mu_radii**2/np.sqrt(self.sigma_radii**2 + self.mu_radii**2))
            sigma = np.sqrt(np.log(1+self.sigma_radii**2/self.mu_radii**2))
            self.radii = np.random.lognormal(mu, sigma, self.N_emitters)
        elif(r_dist == 'normal'):
            mu = self.mu_radii
            sigma = self.sigma_radii
            self.radii = np.random.normal(mu, sigma, self.N_emitters)
        elif(r_dist == 'uniform'):
            low = self.mu_radii - np.sqrt(12) * self.sigma_radii * 0.5
            high = self.mu_radii + np.sqrt(12) * self.sigma_radii * 0.5
            self.radii = np.random.uniform(low, high, self.N_emitters)
        else:
            print ("wrong distribution for radii")
            
        self.heights = self.betas * self.radii
        self.currents = np.copy(self.betas) * 0.
        self.areas = np.pi * self.radii**2
        
           
    def emit(self, Ffar):
        for i in range(len(self.betas)):
            self.emitter.F = self.betas[i] * Ffar
            self.emitter.R = self.radii[i]
            self.emitter.cur_dens()
            self.currents[i] = self.emitter.Jem * self.areas[i]
            
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
        
        # for i in range(len(betas)):
            # print heights[i], radii[i], betas[i]
            
        print ("heights = %.3g +- %.3g"%(np.mean(self.heights), np.std(self.heights)))
        print ("radii = %.3g +- %.3g", np.mean(self.radii), "std<radii> =  ", np.std(self.radii))
        print ("betas = %.3g +- %.3g"%(np.mean(self.betas),np.std(self.betas)))
        print ("maxbeta = %g"%(np.max(self.betas)))
        

        
if (__name__ == "__main__"): #some testing operations
    me = MultiEmitter(N=200)
    me.getEmittersFromBetaDistribution(min_beta=50, max_beta=110, beta_dist="normal", r_dist="lognormal")
    
    farFields = np.linspace(0.02, 0.08, 128)
    totalCurrents = np.copy(farFields)
    for i in range(len(farFields)):
        totalCurrents[i] = me.emit(Ffar=farFields[i])

    
    plt.semilogy(1/farFields, totalCurrents)
    plt.show()




