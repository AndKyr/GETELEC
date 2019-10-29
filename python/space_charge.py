"""This module facilitates space charge calculations for spherical, cylindrical
and planar diodes."""


import numpy as np
import scipy.optimize as opt
import scipy.integrate as ig
import sys
import os
sys.path.append("/home/kyritsak/Code/GETELEC/python")
import getelec_mod as gt


V0 = 5.e-10
tol = 1.e-10

    # ode system for planar geometry    
def planar(y, x, kJ):
    dy = y[1]
    ddy = (kJ / np.sqrt(y[0]))
    return np.array([dy, ddy])

#ode system for cylindrical geometry    
def cylindrical(y, x, kJ):
    dy = y[1]
    ddy = (kJ / y[0]**0.5 - y[1]) / x
    return np.array([dy, ddy])
    
#ode system for spherical geometry
def spherical(y, x, kJ):
    dy = y[1]
    ddy = (kJ / np.sqrt(y[0]))/x**2 - 2. * y[1] /x 
    return np.array([dy, ddy])


# systems for the "reduced" variable u = 2/3 V^1.5 / sqrt(kJ)
def spherical_u(y, x, Radius):
    dy = y[1]
    ddy = Radius**2 / (3 * y[0] * x**2) - 2 * y[1] / x - y[1]**2 / (3. * y[0])
    return np.array([dy, ddy])
    
def planar_u(y, x, Radius = None):
    dy = y[1]
    ddy = 1. / (1 * y[0]) - y[1]**2 / (1. * y[0])
    return np.array([dy, ddy])
    
    

    
""" This class accommodates all 1D space charge calculations"""        
class SpaceCharge():
    kappa = 1.904e5
    
        
    #initializes the class
    def __init__(self, Voltage = 100, Distance = 10, geo = "plane", 
                Npoints = 1024, emitter = None, Radius = 1.):
        if (emitter == None):
            print "SpaceCharge class created with default emitter"
            self.emitter = emission_create()
        else:
            self.emitter =  emitter
        
        self.voltage = float(Voltage)
        self.distance = float(Distance)
        self.geometry = geo
        self.radius = float(Radius)
                                                                                  
        if (self.geometry == "plane"):
            self.system = planar
            self.system_u = planar_u
            self.Radius = 1.
            self.beta = 1. / self.distance
        elif(self.geometry == "cylinder"):
            self.system = cylindrical
            self.beta = 1. / (self.radius * np.log(1 + self.distance / self.radius))
        elif (self.geometry == "sphere"):
            self.system = spherical
            self.system_u = spherical_u
            self.beta = (self.radius + self.distance) / (self.radius * self.distance)
        else:
            print "Error: Wrong geometry. Should be either, plane, cylinder or sphere"
        
        self.Np = Npoints
        
    def ksi(self, x):
        if (self.geometry == "plane"):
            return x
        elif(self.geometry == "cylinder"):
            return (self.radius * np.log(1 + x / self.radius))
        elif (self.geometry == "sphere"):
            return (self.radius * x) / (self.radius + x)
        else:
            print "Error: Wrong geometry. Should be either, plane, cylinder or sphere"
        
        
        
    # solves the ode initial value problem with cathode field F and finds potential at d
    def Vanode(self, F, J = None):
        ## calculate current density for given cathode field F
        if (J == None):
            self.emitter.F = F
            self.emitter.cur_dens()
            self.Jc = self.emitter.Jem
        else:
            self.Jc = J
        
        if (self.geometry == "plane"):
            factor = self.Jc * self.kappa
        elif (self.geometry == "sphere"):
            factor = self.kappa * self.Jc * self.radius**2
        elif (self.geometry == "cylinder"):
            factor = self.kappa * self.Jc * self.radius
        
        ## solve ode and calculate V(d)
        xs = np.linspace(self.radius, self.radius + self.distance, self.Np)
        y = ig.odeint(self.system, np.array([V0,F]), xs, (factor,), rtol = tol, tcrit = np.array([0]))
        return y[-1, 0]
        
        
        # solves the ode initial value problem with cathode field F and finds potential at d
    def xeff(self):
        ## solve ode and calculate V(d)
        xs = np.linspace(self.radius, self.radius + self.distance, self.Np)
        y = ig.odeint(self.system_u, np.array([V0,0]), xs, (self.radius,), rtol = tol, tcrit = np.array([0]))
        return xs, y[:,0], y[:,1]
        
    # uses the bisection method to calculate F(V, d)    
    def calcJF(self, V):
        rtol = 1.e-4
        Fmin = 0.
        
        Fmax = self.beta * V
        for i in range(50):
            Fi = .5 * (Fmax + Fmin)
            Vi = self.Vanode(Fi)
            if (Vi < (1 - rtol) * V):
                Fmin = Fi
            elif((Vi > (1 + rtol) * V)):
                Fmax = Fi
            else:
                break
            
            # print "Fi = %.3f, Vi = %.3f, Ji = %e"%(Fi, Vi, self.Jc)
        
        return Fi
        
    def Vanode_cl(self, J = None):
        ## calculate current density for given cathode field F
        if (J != None):
            self.Jc = J
        
        if (self.geometry == "plane"):
            factor = self.Jc * self.kappa
        elif (self.geometry == "sphere"):
            factor = self.kappa * self.Jc * self.radius**2
        elif (self.geometry == "cylinder"):
            factor = self.kappa * self.Jc * self.radius
        
        ## solve ode and calculate V(d)
        xs = np.linspace(self.radius, self.radius + self.distance, self.Np)
        y = ig.odeint(self.system, np.array([V0,0.]), xs, (factor,), rtol = tol, tcrit = np.array([0]))
        return y[-1, 0]
        
    def pot_solution(self, F, J = None):
        ## calculate current density for given cathode field F
        if (J == None):
            self.emitter.F = F
            self.emitter.cur_dens()
            self.Jc = self.emitter.Jem
        else:
            self.Jc = J
        
        if (self.geometry == "plane"):
            factor = self.Jc * self.kappa
        elif (self.geometry == "sphere"):
            factor = self.kappa * self.Jc * self.radius**2
        elif (self.geometry == "cylinder"):
            factor = self.kappa * self.Jc * self.radius
        
        ## solve ode and calculate V(d)
        xs = np.linspace(self.radius, self.radius + self.distance, self.Np)
        y = ig.odeint(self.system, np.array([V0,F]), xs, (factor,), rtol = tol, tcrit = np.array([0]))
        return xs, y[:,0], y[:,1]
        
        
        
            



    
