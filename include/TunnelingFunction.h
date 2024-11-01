#ifndef TUNNELINGFUNCTION_H_
#define TUNNELINGFUNCTION_H_

#include <math.h>
#include <vector>
#include "Utilities.h"
using namespace std;


class TunnelingFunctionBase{
private:
    double energy = 0.;

public:
    TunnelingFunctionBase(){}
    TunnelingFunctionBase(double E) : energy(E){}
    ~TunnelingFunctionBase(){}

    virtual double potentialFunction(double x){ return 0. * x;}

    virtual double potentialFunctionDerivative(double x){return 0. * x;}

    virtual double findLeftXLimit(double maxPotentialDepth){ return 0.;}

    virtual double findRightXLimit(double maxPotentialDepth){ return 0.;}

    double kappaSquared(double x){
        return CONSTANTS.kConstant * (energy - this->potentialFunction(x));
    }

    void setEnergy(double E){energy = E;}

    double getEnergy(){return energy;}
};

class ModifiedSNBarrier : public TunnelingFunctionBase{
private:
    double radius = 1.e3;
    double field = 5.;
    double gamma = 10.;

    double imagePotential(double z){
        return CONSTANTS.imageChargeConstant * ( 1. / z - .5 / radius);
    }

    double imagePotentialDerivative(double z){
        return -CONSTANTS.imageChargeConstant * (1./(z * z));
    }

    double electrostaticPotential(double z){
        return field * (radius * (gamma - 1.) * z + z*z) / (gamma * z + radius * (gamma - 1.));
    }

    double electrostaticPotentialDerivative(double z){
        return field * ((gamma-1)*(gamma-1) * radius*radius + 2*(gamma-1)*radius*z + gamma*z*z) / 
        (((gamma - 1)*radius + gamma*z) * ((gamma - 1)*radius + gamma*z));
    }

public:
    ModifiedSNBarrier(){}
    ModifiedSNBarrier(double f, double r, double g) : field(f), radius(r), gamma(g){}
    ~ModifiedSNBarrier(){}

    void setBarrierParameters(double f, double r, double g){
        field = f;
        radius = r;
        gamma = g;
    } 

    void setField(double f){
        field = f;
    }

    void setRadius(double R){
        radius = R;
    }

    void setGamma(double g){
        gamma = g;
    }

    double potentialFunction(double z) override{
        return -imagePotential(z) - electrostaticPotential(z);
    }

    double potentialFunctionDerivative(double z) override{
        return -imagePotentialDerivative(z) - electrostaticPotentialDerivative(z);
    }
    
    double findLeftXLimit(double maxPotentialDepth) override{
        return 1. / (maxPotentialDepth / CONSTANTS.imageChargeConstant + 0.5 / radius);
    }

    double findRightXLimit(double maxPotentialDepth) override{
        double b = field * radius * (gamma - 1.) - gamma * maxPotentialDepth;
        double c = - maxPotentialDepth * radius * (gamma - 1.);
        return 0.5 * (-b + sqrt(b*b - 4*field*c)) / field;
    }
};

#endif