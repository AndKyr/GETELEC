#ifndef BANDEMITTER_H_
#define BANDEMITTER_H_

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

#include <vector>
#include <typeinfo>
#include <fstream>
#include <string>


#include "TransmissionSolver.h"

using namespace std;



class BandEmitter : public ODESolver{
private:


    /**
     * @brief A struct that keeps all the important parameters of the emitter. It is necessary so that a void pointer is passed to the
     * differentialSystem function as params, to avoid passing "this" object itself (ugly self-reference)
     */
    struct SystemParameters{
        double workFunction = 4.5;
        double kT = 0.025;
        double effectiveMass = 1.;
        double bandDepth = 10.; // the depth of the band as measured from the Fermi level
        TransmissionSolver* transmissionSolver;
    } systemParams;

    static int differentialSystem(double energy, const double y[], double f[], void *params){
        SystemParameters* sysParams = (SystemParameters*) params; //cast the void pointer as SystemParams

        //calculate transmission coefficient
        // TODO: this works only for effectiveMass = 1. FIX IT
        double D = sysParams->transmissionSolver->calculateTransmissionCoefficientForEnergy(-sysParams->workFunction + energy);

        f[0] = Utilities::fermiDiracFunction(energy, sysParams->kT) * (D - y[0] * exp(energy/sysParams->kT) / sysParams->kT);
        f[1] = y[0];
        return GSL_SUCCESS;

    }

    // static int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

public:

    void setParameters(double workFunction = 4.5, double kT = 0.025, double effectiveMass = 1., double bandDepth = 10.){
        systemParams.workFunction = workFunction;
        systemParams.effectiveMass = effectiveMass;
        systemParams.kT = kT;
        systemParams.bandDepth = bandDepth;
        systemParams.transmissionSolver->setXlimits(systemParams.workFunction + systemParams.bandDepth + 2.);
        xInitial = -bandDepth + 0.1;
        xFinal = 10. * systemParams.kT;
        initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;
    }
    

    BandEmitter(TransmissionSolver* solver, 
                int systemDimension = 2, 
                double rtol = 1.e-5,
                double atol = 1.e-12,
                const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
                int maxSteps = 4096,
                int stepExpectedForInitialStep = 64,
                double maxPotentialDepth = 10.
                )   :   ODESolver(vector<double>(2, 0.0), differentialSystem, 2, {-systemParams.bandDepth + 0.1, 10. * systemParams.kT},
                                rtol, atol, stepType, maxSteps, stepExpectedForInitialStep, NULL, &systemParams)
    {
        systemParams.transmissionSolver = solver;
        systemParams.transmissionSolver->setXlimits(systemParams.workFunction + systemParams.bandDepth);
        setParameters();
        setInitialValues({0.,0.});
    }

    int calculateCurrentDensityAndSpectra(double convergenceTolerance = 1.e-6){
        double x = xInitial;
        double dx = initialStep;
        int status;
        reinitialize();
        vector<double> previousSolution;
        for (size_t i = 0; i < maxAllowedSteps; i++){ //loop over blocks
            previousSolution = solutionVector;
            savedSolution.push_back(solutionVector);
            xSaved.push_back(x);
            status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
            bool hasConverged = abs(previousSolution[1] - solutionVector[1]) / abs(previousSolution[1]) < convergenceTolerance;
            if (x == xFinal || status != GSL_SUCCESS || hasConverged)    
                return status;         
        }
        return GSL_CONTINUE; 
    }
};  



#endif /* BANDEMITTER_H_ */
