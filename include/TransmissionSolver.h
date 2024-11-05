#ifndef TRANSMISSIONSOLVER_H_
#define TRANSMISSIONSOLVER_H_



#include <iostream>

#include "ODESolver.h"
#include "TunnelingFunction.h"
#include "Utilities.h"

using namespace std;


class TransmissionSolver : public ODESolver{
private:
    TunnelingFunction* tunnelingFunction;
     
    double kappaInitial;
    double kappaFinal;

    static int tunnelingDifferentialSystem(double x, const double y[], double f[], void *params);

    static int tunnelingSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

    int numberOfCalls = 0;

public:   

    TransmissionSolver(
            TunnelingFunction* tunnelFunctionPtr, 
            double relativeTolerance = 1.e-4,
            double absoluteTolerance = 1.e-4,
            const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
            int maxSteps = 4096,
            int stepExpectedForInitialStep = 64,
            double maxPotentialDepth = 10.
        );


    void setXlimits(double maxPotentialDepth){
        xInitial = tunnelingFunction->findRightXLimit(maxPotentialDepth);
        xFinal = tunnelingFunction->findLeftXLimit(maxPotentialDepth);
    }

    void updateKappaAtLimits(){
        kappaInitial = sqrt(tunnelingFunction->kappaSquared(xInitial));
        kappaFinal = sqrt(tunnelingFunction->kappaSquared(xFinal));
        initialValues = {0., kappaInitial, 0.};
    }

    void setEnergy(double E){
        tunnelingFunction->setEnergy(E);
        updateKappaAtLimits();
    }

    void resetNumberOfCalls(){numberOfCalls = 0;}

    double transmissionCoefficient() const;

    double calculateTransmissionCoefficientForEnergy(double energy){
        setEnergy(energy);
        solveNoSave();
        numberOfCalls++;
        return transmissionCoefficient();
    }

    int getNumberOfCalls(){return numberOfCalls;}

    void printXLimits(){ cout << "xInitial = " << xInitial << " xFinal = " << xFinal << endl;}
};

#endif /* TRANSMISSIONSOLVER_H_ */
