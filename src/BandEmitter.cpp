#include "BandEmitter.h"
#include "Utilities.h"
#include <cmath>
#include <iostream>

int BandEmitter::differentialSystem(double energy, const double y[], double f[], void *params) {
    BandEmitter* emitter = (BandEmitter*) params; // Cast the void pointer as SystemParams
    double D = emitter->calculateIntegrand(energy);
    f[0] = Utilities::fermiDiracFunction(energy, emitter->kT) * 
           (D - y[0] * exp(energy / emitter->kT) / emitter->kT);
    f[1] = y[0];
    f[2] = energy * y[0];
    return GSL_SUCCESS;
}

int BandEmitter::differentialSystemLog(double energy, const double y[], double f[], void *params) {
    BandEmitter* emitter = (BandEmitter*) params; // Cast the void pointer as SystemParams
    double D = emitter->calculateIntegrand(energy); // Calculate transmission coefficient
    f[0] = Utilities::fermiDiracFunction(energy, emitter->kT) * 
           (D * exp(-y[0]) - exp(energy / emitter->kT) / emitter->kT);
    f[1] = exp(y[0] - y[1]);
    return GSL_SUCCESS;
}

double BandEmitter::normalEnergyDistribution(double energy, void* params) {
    BandEmitter* emitter = (BandEmitter*) params; // Cast the void pointer as SystemParams
    double result = emitter->calculateIntegrand(energy);
    return result * Utilities::logFermiDiracFunction(energy, emitter->kT);
}

void BandEmitter::updateBarrier() {
    double minNormalEnergy = xInitial;
    if (effectiveMass > 1.)
        minNormalEnergy += (1. - effectiveMass) * (xFinal - xInitial);
    
    // Set the limits for the transmission solver. The argument counts the barrier depth from Evacuum.
    transmissionSolver.setXlimits(workFunction - minNormalEnergy + 2.); 
    
    // Initialize the interpolator with the new limits. Slightly extends th interpolator limits to avoid edge effects.
    interpolator.initialize(minNormalEnergy - 0.1, xFinal + 0.1, 8);
    interpolator.refineToTolerance();
}

double BandEmitter::calculateIntegrand(double energy) {
    //TODO: be careful with the effective mass. abarX might get below the bandDepth causing problems. The interpolation range must be fixed.
    double result = interpolator.evaluate(energy);
    if (effectiveMass != 1.) {
        double aBarX = -bandDepth + (1. - effectiveMass) * (energy + bandDepth);
        result -= (1. - effectiveMass) * interpolator.evaluate(aBarX);
    }
    return result;
}

void BandEmitter::setParameters(double workFunction_, double kT_, double effectiveMass_, double bandDepth_){
    //set the new parameters
    workFunction = workFunction_;
    bandDepth = bandDepth_;
    effectiveMass = effectiveMass_;
    kT = kT_;
    interpolator.setParameters(kT, workFunction);

    //set ODE integration limits and affected parameters
    xInitial = -bandDepth;
    xFinal = workFunction + 10. * kT;
    maxStepSize = (xFinal - xInitial) / minAllowedSteps;
    initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;

    updateBarrier();
    setInitialValues({0., 0., 0.});
}

int BandEmitter::calculateCurrentDensityAndSpectra(double convergenceTolerance) {
    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();
    vector<double> previousSolution;

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousSolution = solutionVector;
        savedSolution.push_back(solutionVector);
        xSaved.push_back(x);
        dx = GSL_SIGN(dx) * min(abs(maxStepSize), abs(dx));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
        bool hasConverged = abs(previousSolution[1] - solutionVector[1]) / previousSolution[1] < convergenceTolerance;
        if (x >= xFinal || status != GSL_SUCCESS || hasConverged)    
            return status;         
    }
    return GSL_CONTINUE; 
}

int BandEmitter::calculateCurrentDensityAndNottingham(double convergenceTolerance) {
    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();
    vector<double> previousSolution;

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousSolution = solutionVector;
        dx = GSL_SIGN(dx) * min(abs(maxStepSize), abs(dx));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
        bool hasConverged = abs(previousSolution[1] - solutionVector[1]) / previousSolution[1] < convergenceTolerance;
        if (x >= xFinal || status != GSL_SUCCESS || hasConverged)    
            return status;         
    }
    return GSL_CONTINUE; 
}

double BandEmitter::calcualteCurrentDensity() {
    double result, error;
    if (!integrationWorkspace)
        integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);

    gsl_integration_qag(&integrationFunction, xInitial, xFinal, 0., relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result * CONSTANTS.SommerfeldConstant * kT;
}

void BandEmitter::writeSavedLogD(string filename) {
    ofstream outFile(filename, ios::out);        
    for (size_t i = 0; i < savedEnergies.size(); i++)
        outFile << savedEnergies[i] << " " << savedLogD[i] << endl;
}

void BandEmitter::writePlottingData(string filename) {
    ofstream outFile(filename, ios::out);        
    outFile << " E D_calc D_interp error NED lFD(E) tolerance" << endl;
    for (double x = xInitial; x < xFinal; x += 0.001) {
        double D = transmissionSolver.calculateTransmissionCoefficientForEnergy(x - workFunction);
        double lFD = Utilities::logFermiDiracFunction(x, kT);
        double err = interpolator.calculateError(x, log(D));
        double tol = interpolator.calculateTolerance(x, log(D));
        outFile << x << " " << D << " " << interpolator.evaluate(x) << " " << err << " " << lFD * D << " " << lFD << " " << tol << endl;
    }
    interpolator.writeSplineNodes();
}
