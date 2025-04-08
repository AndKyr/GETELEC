#include "BandEmitter.h"
#include "Utilities.h"
#include <cmath>
#include <iostream>

namespace getelec{

int BandEmitter::differentialSystem(double energy, const double y[], double f[], void *params) {
    BandEmitter* emitter = (BandEmitter*) params; // Cast the void pointer as "this"
    double D = emitter->gPrimeFunction(energy);
    assert(isfinite(D) && "D is not finite in TED ODE system");
    f[0] = Utilities::fermiDiracFunction(energy, emitter->kT) * 
           (D - y[0] * exp(energy / emitter->kT) / emitter->kT);
    f[1] = y[0];
    f[2] = energy * y[0];
    return GSL_SUCCESS;
}

double BandEmitter::normalEnergyDistribution(double energy, void* params) {
    BandEmitter* emitter = (BandEmitter*) params; // Cast the void pointer as "this"
    double result = emitter->gPrimeFunction(energy);
    return result * Utilities::logFermiDiracFunction(energy, emitter->kT);
}

double BandEmitter::gPrimeFunction(double energy) {
    //TODO: be careful with the effective mass. abarX might get below the bandDepth causing problems. The interpolation range must be fixed.
    double waveVector = sqrt(energy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
    double result = interpolator.getTransmissionProbability(energy - workFunction, waveVector);
    if (effectiveMass != 1.) {
        double aBarX = -bandDepth + (1. - effectiveMass) * (energy + bandDepth);
        result -= (1. - effectiveMass) * interpolator.getTransmissionProbability(aBarX - workFunction, waveVector);
    }
    assert(isfinite(result) && "Transmission coefficient is not finite");
    return result;
}

void BandEmitter::setParameters(double workFunction_, double kT_, double effectiveMass_, double bandDepth_) {
    //set the new parameters
    workFunction = workFunction_;
    bandDepth = bandDepth_;
    effectiveMass = effectiveMass_;
    kT = kT_;

    interpolator.setParameters(kT, workFunction);
    interpolator.smartInitialSampling();
    interpolator.refineSamplingToTolerance();

    //set ODE integration limits and affected parameters
    xInitial = max(interpolator.getMinimumSampleEnergy() + workFunction, -bandDepth + 0.001);
    xFinal = interpolator.getMaximumSampleEnergy() + workFunction;
    maxStepSize = (xFinal - xInitial) / minAllowedSteps;
    initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;

    setInitialValues({0., 0., 0.});
}

int BandEmitter::integrateTotalEnergyDistributionODEAndSaveSpectra(double convergenceTolerance, bool makeSpectralSpline) {
    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();
    savedSpectra.clear();
    savedSpectraDerivative.clear();
    double previousCurrentDensity;

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousCurrentDensity = solutionVector[1];
        savedSpectra.push_back(solutionVector[0] * CONSTANTS.SommerfeldConstant);
        savedSpectraDerivative.push_back(getSolutionDerivative(x)[0] * CONSTANTS.SommerfeldConstant);
        xSaved.push_back(x);
        dx = GSL_SIGN(dx) * min(abs(maxStepSize), abs(dx));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
        bool hasConverged = abs(previousCurrentDensity - solutionVector[1]) / previousCurrentDensity < convergenceTolerance;
        if (x >= xFinal || status != GSL_SUCCESS || hasConverged)    
            return status;         
    }

    if (makeSpectralSpline) updateSpectraSpline();
    return GSL_CONTINUE; 
}

void BandEmitter::updateSpectraSpline(){
    spectraSpline.initialize(xSaved, savedSpectra, savedSpectraDerivative);
}

int BandEmitter::integrateTotalEnergyDistributionODE(double convergenceTolerance) {
    if (convergenceTolerance <= 0.) convergenceTolerance = relativeTolerance;
    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();
    double previousCurrentDensity;

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousCurrentDensity = solutionVector[1];
        dx = GSL_SIGN(dx) * min(abs(maxStepSize), abs(dx));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
        bool hasConverged = abs(previousCurrentDensity - solutionVector[1]) / previousCurrentDensity < convergenceTolerance;
        if (x >= xFinal || status != GSL_SUCCESS || hasConverged)    
            return status;         
    }
    return GSL_CONTINUE; 
}

double BandEmitter::integrateNormalEnergyDistribution() {
    double result, error;
    if (!integrationWorkspace)
        integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);

    gsl_function integrationFunction = {&normalEnergyDistribution, this};
    gsl_integration_qag(&integrationFunction, xInitial, xFinal, 0., relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result * CONSTANTS.SommerfeldConstant * kT;
}


void BandEmitter::writePlottingData(string filename) {
    ofstream outFile(filename, ios::out);        
    outFile << " E D_calc D_interp NED TED lFD(E)" << endl;
    for (double x = xInitial; x < xFinal; x += 0.001) {
        double D = transmissionSolver.calculateTransmissionProbability(x - workFunction);
        double lFD = Utilities::logFermiDiracFunction(x, kT);
        outFile << x << " " << D << " " << interpolator.evaluate(x) << " " << interpolator.emissionCurrentEstimate(x) << " " << spectraForEnergy(x) << " " << lFD << " " << endl;
    }
    interpolator.writeSplineNodes();
}


double BandEmitter::calculateParallelEnergyDistribution(double parallelEnergy) {
    double result, error;
    if (!integrationWorkspace)
        integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);
    
    gsl_function integrationFunction = {&normalEnergyDistribution, this};

    gsl_integration_qag(&integrationFunction, xInitial, xFinal, 0., relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result * CONSTANTS.SommerfeldConstant * kT;

}
}