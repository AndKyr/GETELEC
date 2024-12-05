#include "BandEmitter.h"
#include "Utilities.h"
#include <cmath>
#include <iostream>

int BandEmitter::differentialSystem(double energy, const double y[], double f[], void *params) {
    BandEmitter* sysParams = (BandEmitter*) params; // Cast the void pointer as SystemParams
    double D = sysParams->calculateIntegrand(energy);
    f[0] = Utilities::fermiDiracFunction(energy, sysParams->kT) * 
           (D - y[0] * exp(energy / sysParams->kT) / sysParams->kT);
    f[1] = y[0];
    f[2] = energy * y[0];
    return GSL_SUCCESS;
}

int BandEmitter::differentialSystemLog(double energy, const double y[], double f[], void *params) {
    BandEmitter* sysParams = (BandEmitter*) params; // Cast the void pointer as SystemParams
    double D = sysParams->calculateIntegrand(energy); // Calculate transmission coefficient
    f[0] = Utilities::fermiDiracFunction(energy, sysParams->kT) * 
           (D * exp(-y[0]) - exp(energy / sysParams->kT) / sysParams->kT);
    f[1] = exp(y[0] - y[1]);
    return GSL_SUCCESS;
}

double BandEmitter::normalEnergyDistribution(double energy, void* params) {
    BandEmitter* sysParams = (BandEmitter*) params; // Cast the void pointer as SystemParams
    double result = log(sysParams->calculateIntegrand(energy));
    sysParams->savedEnergies.push_back(energy);
    sysParams->savedLogD.push_back(result);
    return result;
}

void BandEmitter::updateBarrier() {
    transmissionSolver.setXlimits(workFunction + bandDepth + 1.);
    interpolator.initialize(-bandDepth, 10 * kT + workFunction, 8);
    interpolator.refineToTolerance();
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

double BandEmitter::calcualteCurrentDensity() {
    double result, error;
    if (!integrationWorkspace)
        integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);

    gsl_integration_qag(&integrationFunction, xInitial, xFinal, 0., relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result;
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
