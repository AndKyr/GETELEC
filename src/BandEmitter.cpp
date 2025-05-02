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

double BandEmitter::gPrimeFunction(double energy) {
    //TODO: be careful with the effective mass. abarX might get below the bandDepth causing problems. The interpolation range must be fixed.
    double waveVector = sqrt(energy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
    double result = interpolator.getTransmissionProbability(energy - workFunction, waveVector);
    assert(isfinite(result) && "Transmission coefficient is not finite");

    if (effectiveMass != 1.) {
        double reducedEnergy = -bandDepth + (1. - effectiveMass) * (energy + bandDepth);
        double reducedTransmission = interpolator.getTransmissionProbability(reducedEnergy - workFunction, waveVector);
        result -= (1. - effectiveMass) * reducedTransmission;

        // Debug code
        if (! isfinite(result))
            interpolator.getTransmissionProbability(reducedEnergy - workFunction, waveVector);
        
        // end of debug code


        assert(isfinite(result) && "Transmission coefficient is not finite");
    }
    return result;
}

void BandEmitter::setParameters(double workFunction_, double kT_, double effectiveMass_, double bandDepth_, bool doQuadrature) {
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

    if (!integrationWorkspace && doQuadrature)
        integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);

    setInitialValues({0., 0., 0.});
}

int BandEmitter::integrateTotalEnergyDistributionODEAndSaveSpectra(double convergenceTolerance, bool makeSpectralSpline) {
    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();

    totalEnergyDistribution.first.reserve(maxAllowedSteps);
    totalEnergyDistribution.second.reserve(maxAllowedSteps);
    totalEnergyDistributionDerivatives.reserve(maxAllowedSteps);
    totalEnergyDistributionDerivatives.clear();
    totalEnergyDistribution.first.clear();
    totalEnergyDistribution.second.clear();

    
    double previousCurrentDensity;
    if(convergenceTolerance <= 0.) convergenceTolerance = relativeTolerance;

    totalEnergyDistribution.second.push_back(solutionVector[0] * CONSTANTS.SommerfeldConstant);
    totalEnergyDistributionDerivatives.push_back(getSolutionDerivative(x)[0] * CONSTANTS.SommerfeldConstant);
    totalEnergyDistribution.first.push_back(x);

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousCurrentDensity = solutionVector[1];
        dx = GSL_SIGN(dx) * min(abs(maxStepSize), abs(dx));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());

        totalEnergyDistribution.second.push_back(solutionVector[0] * CONSTANTS.SommerfeldConstant);
        totalEnergyDistributionDerivatives.push_back(getSolutionDerivative(x)[0] * CONSTANTS.SommerfeldConstant);
        totalEnergyDistribution.first.push_back(x);

        bool hasConverged = abs(previousCurrentDensity - solutionVector[1]) / previousCurrentDensity < convergenceTolerance;
        if (x >= xFinal || status != GSL_SUCCESS || hasConverged){
            if (makeSpectralSpline) totalEnergyDistributionSpline.initialize(totalEnergyDistribution.first, totalEnergyDistribution.second, totalEnergyDistributionDerivatives);   
            return status;         
        }
    }

    if (makeSpectralSpline) totalEnergyDistributionSpline.initialize(totalEnergyDistribution.first, totalEnergyDistribution.second, totalEnergyDistributionDerivatives);
    return GSL_CONTINUE; 
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

double BandEmitter::currentDensityIntegrateNormal(bool saveNormalEnergyDistribution) {
    if (!externalIntegrationWorkSpace)
        externalIntegrationWorkSpace = gsl_integration_workspace_alloc(maxAllowedSteps);
    assert(externalIntegrationWorkSpace && "externalIntegrationWorkSpace is not initialized");

    if (saveNormalEnergyDistribution){
        normalEnergyDistribution.first.reserve(maxAllowedSteps);
        normalEnergyDistribution.first.clear();
        normalEnergyDistribution.second.reserve(maxAllowedSteps);
        normalEnergyDistribution.second.clear();
        saveSpectra = true;
    } else
        saveSpectra = false;

    auto integrationLambda = [](double normalEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        double result = thisObj->normalEnergyDistributionForEnergy(normalEnergy);
        if (thisObj->saveSpectra){
            thisObj->normalEnergyDistribution.first.push_back(normalEnergy);
            thisObj->normalEnergyDistribution.second.push_back(result);
        }
        return result;
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, xInitial, xFinal, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    
    if (saveNormalEnergyDistribution)
        Utilities::sortCoordinates(normalEnergyDistribution.first, normalEnergyDistribution.second);
                            
    return result * kT * CONSTANTS.SommerfeldConstant;
}

void BandEmitter::writePlottingData(string filename) {
    ofstream outFile(filename, ios::out);        
    outFile << "#E D_calc(En=E) D_interp(En=E) NED(En=E) TED(Et=E)" << endl;
    for (double energy = xInitial; energy < xFinal; energy += 0.001) {
        double waveVector = sqrt(energy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
        double D_calculated = transmissionSolver.calculateTransmissionProbability(energy - workFunction, waveVector);
        double D_interpolated = interpolator.getTransmissionProbability(energy - workFunction, waveVector);
        double NED = interpolator.normalEnergyDistributionEstimate(energy - workFunction, waveVector);
        double TED = getTotalEnergyDistributionForEnergy(energy);

        outFile << energy << " " << D_calculated << " " << D_interpolated << " " << NED << " " << TED << " " << endl;
    }
    interpolator.writeSplineNodes();
}

double BandEmitter::doubleIntegrandTotalParallel(double totalEnergy, double parallelEnergy) const{
    double waveVectorZ = sqrt((totalEnergy+bandDepth) * effectiveMass - parallelEnergy) * CONSTANTS.sqrt2mOverHbar;
    assert(waveVectorZ > 0 && "waveVectorZ is not positive");

    double normalEnergy = totalEnergy - parallelEnergy - workFunction;
    
    // Check if normalEnergy is within the valid range. The interpolator should have informed and valid energy ranges (even for effectiveMass!=1)
    if (normalEnergy < interpolator.getMinimumSampleEnergy() || normalEnergy > interpolator.getMaximumSampleEnergy())
        return 0.;
    
    double transmissionProbability = interpolator.getTransmissionProbability(normalEnergy, waveVectorZ);
    return Utilities::fermiDiracFunction(totalEnergy, kT) * transmissionProbability;
}

double BandEmitter::doubleIntegrandParallelNormal(double parallelEnergy, double normalEnergy) const {
    double totalEnergy = normalEnergy + parallelEnergy;
    double waveVectorZ = sqrt((totalEnergy+bandDepth) * effectiveMass - parallelEnergy) * CONSTANTS.sqrt2mOverHbar;
    assert(waveVectorZ > 0 && "waveVectorZ is not positive");

    double normalEnegyForInterpolator = normalEnergy - workFunction;
    
    // Check if normalEnergy is within the valid range. The interpolator should have informed and valid energy ranges (even for effectiveMass!=1)
    if (normalEnegyForInterpolator < interpolator.getMinimumSampleEnergy() || normalEnegyForInterpolator > interpolator.getMaximumSampleEnergy())
        return 0.;
    
    double transmissionProbability = interpolator.getTransmissionProbability(normalEnegyForInterpolator, waveVectorZ);
    return Utilities::fermiDiracFunction(totalEnergy, kT) * transmissionProbability;
}

double BandEmitter::totalEnergyDistributionIntegrateParallel(double totalEnergy){
    assert(integrationWorkspace && "integrationWorkspace is not initialized");


    helperEnergy = totalEnergy; // store the totalEnergy in the class variable to be used in the lambda function
    auto integrationLambda = [](double parallelEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        return thisObj->doubleIntegrandTotalParallel(thisObj->helperEnergy, parallelEnergy);
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};

    double maxParallelEnergy = effectiveMass * (totalEnergy + bandDepth);
    assert(maxParallelEnergy >= 0 && "maxParallelEnergy is not positive");

    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, 0., maxParallelEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result * CONSTANTS.SommerfeldConstant;
}

double BandEmitter::currentDensityIntegrateTotalParallel(bool saveTotalEnergyDistribution){

    if (!externalIntegrationWorkSpace)
        externalIntegrationWorkSpace = gsl_integration_workspace_alloc(maxAllowedSteps);
    assert(externalIntegrationWorkSpace && "externalIntegrationWorkSpace is not initialized");

    if (saveTotalEnergyDistribution){
        totalEnergyDistribution.first.reserve(maxAllowedSteps);
        totalEnergyDistribution.first.clear();
        totalEnergyDistribution.second.reserve(maxAllowedSteps);
        totalEnergyDistribution.second.clear();
        saveSpectra = true;
    } else{
        saveSpectra = false;
    }

    auto integrationLambda = [](double totalEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        double result = thisObj->totalEnergyDistributionIntegrateParallel(totalEnergy);
        if (thisObj->saveSpectra){
            thisObj->totalEnergyDistribution.first.push_back(totalEnergy);
            thisObj->totalEnergyDistribution.second.push_back(result);
        }
        return result;
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, xInitial, xFinal, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    if (saveSpectra)
        Utilities::sortCoordinates(totalEnergyDistribution.first, totalEnergyDistribution.second);
    
    
    return result;
}

double BandEmitter::nottinghamIntegrateTotalPrallel(){
    if (!externalIntegrationWorkSpace)
        externalIntegrationWorkSpace = gsl_integration_workspace_alloc(maxAllowedSteps);
    assert(externalIntegrationWorkSpace && "externalIntegrationWorkSpace is not initialized");

    auto integrationLambda = [](double totalEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        return totalEnergy * thisObj->totalEnergyDistributionIntegrateParallel(totalEnergy);
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, xInitial, xFinal, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    return result;
}

double BandEmitter::parallelEnergyDistributionForEnergy(double parallelEnergy){
    assert(integrationWorkspace && "integrationWorkspace is not initialized");


    helperEnergy = parallelEnergy; // store the totalEnergy in the class variable to be used in the lambda function
    auto integrationLambda = [](double totalEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        return thisObj->doubleIntegrandTotalParallel(totalEnergy, thisObj->helperEnergy);
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this}; //convert to gsl_function and pass this pointer

    double minTotalEnergy = -bandDepth + parallelEnergy / effectiveMass;

    if (minTotalEnergy > xFinal)
        return 0.;

    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minTotalEnergy, xFinal, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result * CONSTANTS.SommerfeldConstant;
}

double BandEmitter::currentDensityIntegrateParallelTotal(bool saveParallelEnergyDistribution){
    if (!externalIntegrationWorkSpace)
        externalIntegrationWorkSpace = gsl_integration_workspace_alloc(maxAllowedSteps);
    assert(externalIntegrationWorkSpace && "externalIntegrationWorkSpace is not initialized");

    if (saveParallelEnergyDistribution){
        parallelEnergyDistribution.first.reserve(maxAllowedSteps);
        parallelEnergyDistribution.first.clear();
        parallelEnergyDistribution.second.reserve(maxAllowedSteps);
        parallelEnergyDistribution.second.clear();
        saveSpectra = true;
    } else 
        saveSpectra = false;

    auto integrationLambda = [](double parallelEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        double result = thisObj->parallelEnergyDistributionForEnergy(parallelEnergy);
        if (thisObj->saveSpectra){
            thisObj->parallelEnergyDistribution.first.push_back(parallelEnergy);
            thisObj->parallelEnergyDistribution.second.push_back(result);
        }
        return result;
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, 0., xFinal - xInitial, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    
    if (saveParallelEnergyDistribution) 
        Utilities::sortCoordinates(parallelEnergyDistribution.first, parallelEnergyDistribution.second);
    return result;
}

double BandEmitter::normalEnergyDistributionIntegratePrallel(double normalEnergy){
    assert(integrationWorkspace && "integrationWorkspace is not initialized");


    helperEnergy = normalEnergy; // store the totalEnergy in the class variable to be used in the lambda function
    auto integrationLambda = [](double parallelEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        return thisObj->doubleIntegrandParallelNormal(parallelEnergy, thisObj->helperEnergy);
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};

    double maxParalleEnergy = min(xFinal - normalEnergy, (normalEnergy + bandDepth) * effectiveMass / (1. - effectiveMass));
    double minParallelEnergy = max(xInitial - normalEnergy, 0.);

    if (minParallelEnergy > maxParalleEnergy)
        return 0.;

    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minParallelEnergy, maxParalleEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
    return result * CONSTANTS.SommerfeldConstant;
}

double BandEmitter::currentDensityIntegrateNormalParallel(bool saveNormalEnergyDistribution){
    if (!externalIntegrationWorkSpace)
        externalIntegrationWorkSpace = gsl_integration_workspace_alloc(maxAllowedSteps);
    assert(externalIntegrationWorkSpace && "externalIntegrationWorkSpace is not initialized");

    if (saveNormalEnergyDistribution){
        normalEnergyDistribution.first.reserve(maxAllowedSteps);
        normalEnergyDistribution.first.clear();
        normalEnergyDistribution.second.reserve(maxAllowedSteps);
        normalEnergyDistribution.second.clear();
        saveSpectra = true;
    } else
        saveSpectra = false;

    auto integrationLambda = [](double normalEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        double result = thisObj->normalEnergyDistributionIntegratePrallel(normalEnergy);
        if (thisObj->saveSpectra){
            thisObj->normalEnergyDistribution.first.push_back(normalEnergy);
            thisObj->normalEnergyDistribution.second.push_back(result);
        }
        return result;
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, xInitial, xFinal, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    
    if (saveNormalEnergyDistribution) 
        Utilities::sortCoordinates(normalEnergyDistribution.first, normalEnergyDistribution.second);
    
    
    return result;
}

}