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

void BandEmitter::setIntegrationLimits(){
    if (effectiveMass > 0.){ // conduction band curving upwards
        //pick the max of the band bottom and the min energy that the transmission coefficient is non-negligible
        minTotalEnergy = max(interpolator.getMinimumSampleEnergy() + workFunction, -bandDepth + 0.001); 
        maxTotalEnergy = interpolator.getMaximumSampleEnergy() + workFunction; //set the max to the maximum sample
    } else { //valence band (curving downwards)
        minTotalEnergy = interpolator.getMinimumSampleEnergy() + workFunction; //go to the minimum that has non-negligible transmission probability
        maxTotalEnergy = -bandDepth; //set the maximum to the band top
    }
    //set the inhereted parameters to follow for clarity
    xInitial = minTotalEnergy; 
    xFinal = maxTotalEnergy;    
    
    if (effectiveMass <= 1.)
        minNormalEnergy = minTotalEnergy;// -bandDepth + (1. - effectiveMass) * (minTotalEnergy + bandDepth);
    else
        minNormalEnergy = max(-bandDepth + (1. - effectiveMass) * (maxTotalEnergy + bandDepth), interpolator.getMinimumSampleEnergy() + workFunction);
    maxNormalEnergy = maxTotalEnergy;

    minParallelEnergy = 0.; //min parallel energy is always 0
    if (effectiveMass > 0.)
        maxParallelEnergy = min(effectiveMass * (maxTotalEnergy + bandDepth), maxNormalEnergy - minNormalEnergy);
    else
        maxParallelEnergy = effectiveMass * (minTotalEnergy + bandDepth);
}

double BandEmitter::gPrimeFunction(double normalEnergy) {
    assert(abs(effectiveMass - 1.) < numeric_limits<double>::epsilon() && "the g'(E) function is valid only for m* = 1");
    double waveVector = sqrt(normalEnergy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
    double result = interpolator.getTransmissionProbability(normalEnergy - workFunction, waveVector);
    assert(isfinite(result) && "Transmission coefficient is not finite");

    // if (effectiveMass != 1.) {
    //     double reducedEnergy = -bandDepth + (1. - effectiveMass) * (normalEnergy + bandDepth);
    //     double reducedTransmission = interpolator.getTransmissionProbability(reducedEnergy - workFunction, waveVector);
    //     result -= (1. - effectiveMass) * reducedTransmission;
    //     assert(isfinite(result) && "Transmission coefficient is not finite");
    // }
    return result;
}

void BandEmitter::setParameters(double workFunction_, double kT_, double effectiveMass_, double bandDepth_, bool doQuadrature) {
    //set the new parameters
    workFunction = workFunction_;
    bandDepth = bandDepth_;
    effectiveMass = effectiveMass_;
    kT = kT_;
    assert(bandDepth > 0. && kT > 0. && workFunction > 0. && abs(effectiveMass) > 1.e-5 && "Invalid BandEmitter input parameters");

    transmissionSolver.setWriteFlag(writePlottingFiles);
    interpolator.setParameters(kT, workFunction);
    interpolator.smartInitialSampling();
    interpolator.refineSamplingToTolerance();

    if (writePlottingFiles >= 0){
        interpolator.writeSplineNodes("splineNodes_i_" + to_string(writePlottingFiles) + ".dat");
        interpolator.writeSplineSolution("splineSolution_i_" + to_string(writePlottingFiles) + ".dat");
    }

    setIntegrationLimits();

    maxStepSize = (xFinal - xInitial) / minAllowedSteps;
    initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;

    if (!integrationWorkspace && doQuadrature)
        integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);

    setInitialValues({0., 0., 0.});
}

int BandEmitter::integrateTotalEnergyDistributionODEAndSaveSpectra(double convergenceTolerance, bool makeSpectralSpline) {
    double energy = minTotalEnergy;
    double energyStep = initialStep;
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
    totalEnergyDistributionDerivatives.push_back(getSolutionDerivative(energy)[0] * CONSTANTS.SommerfeldConstant);
    totalEnergyDistribution.first.push_back(energy);

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousCurrentDensity = solutionVector[1];
        energyStep = GSL_SIGN(energyStep) * min(abs(maxStepSize), abs(energyStep));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &energy, maxTotalEnergy, &energyStep, solutionVector.data());

        totalEnergyDistribution.second.push_back(solutionVector[0] * CONSTANTS.SommerfeldConstant);
        totalEnergyDistributionDerivatives.push_back(getSolutionDerivative(energy)[0] * CONSTANTS.SommerfeldConstant);
        totalEnergyDistribution.first.push_back(energy);

        bool hasConverged = abs(previousCurrentDensity - solutionVector[1]) / previousCurrentDensity < convergenceTolerance;
        if (energy >= maxTotalEnergy || status != GSL_SUCCESS || hasConverged){
            if (makeSpectralSpline) totalEnergyDistributionSpline.initialize(totalEnergyDistribution.first, totalEnergyDistribution.second, totalEnergyDistributionDerivatives);   
            return status;         
        }
    }

    if (makeSpectralSpline) totalEnergyDistributionSpline.initialize(totalEnergyDistribution.first, totalEnergyDistribution.second, totalEnergyDistributionDerivatives);
    return GSL_CONTINUE; 
}

int BandEmitter::integrateTotalEnergyDistributionODE(double convergenceTolerance) {
    if (convergenceTolerance <= 0.) convergenceTolerance = relativeTolerance;
    double energy = minTotalEnergy;
    double energyStep = initialStep;
    int status;
    reinitialize();
    double previousCurrentDensity;

    for (size_t i = 0; i < maxAllowedSteps; i++) { // Loop over blocks
        previousCurrentDensity = solutionVector[1];
        energyStep = GSL_SIGN(energyStep) * min(abs(maxStepSize), abs(energyStep));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &energy, maxTotalEnergy, &energyStep, solutionVector.data());
        bool hasConverged = abs(previousCurrentDensity - solutionVector[1]) / previousCurrentDensity < convergenceTolerance;
        if (energy >= maxTotalEnergy || status != GSL_SUCCESS || hasConverged)    
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
            thisObj->normalEnergyDistribution.second.push_back(result * CONSTANTS.SommerfeldConstant) ;
        }
        return result;
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minNormalEnergy, maxNormalEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    
    if (saveNormalEnergyDistribution)
        Utilities::sortCoordinates(normalEnergyDistribution.first, normalEnergyDistribution.second);
                            
    return result * kT * CONSTANTS.SommerfeldConstant;
}

void BandEmitter::writePlottingData() {
    interpolator.writeSplineSolution("splineSolution_i_" + to_string(writePlottingFiles) + ".dat", 256, true);
    interpolator.writeSplineNodes("splineNodes_" + to_string(writePlottingFiles) + ".dat");
}

double BandEmitter::doubleIntegrandTotalParallel(double totalEnergy, double parallelEnergy) const{
    double waveVectorZ = getWaveVectorZ(totalEnergy, parallelEnergy);
    
    double normalEnergy = totalEnergy - parallelEnergy - workFunction;
    // Check if normalEnergy is within the valid range. The interpolator should have informed and valid energy ranges (even for effectiveMass!=1)
    assert(normalEnergy >= interpolator.getMinimumSampleEnergy() && normalEnergy <= interpolator.getMaximumSampleEnergy() && "integration limits out of bounds");
    
    double transmissionProbability = interpolator.getTransmissionProbability(normalEnergy, waveVectorZ);
    return Utilities::fermiDiracFunction(totalEnergy, kT) * transmissionProbability;
}

double BandEmitter::doubleIntegrandParallelNormal(double parallelEnergy, double normalEnergy) const {
    double totalEnergy = normalEnergy + parallelEnergy;
    double waveVectorZ = getWaveVectorZ(totalEnergy, parallelEnergy);
    
    double normalEnegyForInterpolator = normalEnergy - workFunction;
    // Check if normalEnergy is within the valid range. The interpolator should have informed and valid energy ranges (even for effectiveMass!=1)
    assert(normalEnegyForInterpolator > interpolator.getMinimumSampleEnergy() && normalEnegyForInterpolator < interpolator.getMaximumSampleEnergy() && "integration limits out of bounds");
    
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

    double maxParallelEnergy = min(effectiveMass * (totalEnergy + bandDepth), totalEnergy - workFunction - interpolator.getMinimumSampleEnergy());
    assert(maxParallelEnergy >= -numeric_limits<double>::epsilon() && "maxParallelEnergy is not positive");
    if (maxParallelEnergy < absoluteTolerance) return 0.;

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
    gsl_integration_qag(&gslIntegrationFunction, minTotalEnergy, maxTotalEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
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
    gsl_integration_qag(&gslIntegrationFunction, minTotalEnergy, maxTotalEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    return result;
}

double BandEmitter::parallelEnergyDistributionIntegrateTotal(double parallelEnergy){
    assert(integrationWorkspace && "integrationWorkspace is not initialized");

    helperEnergy = parallelEnergy; // store the totalEnergy in the class variable to be used in the lambda function
    auto integrationLambda = [](double totalEnergy, void* params) {
        BandEmitter* thisObj =  static_cast<BandEmitter*>(params);
        return thisObj->doubleIntegrandTotalParallel(totalEnergy, thisObj->helperEnergy);
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this}; //convert to gsl_function and pass this pointer


    double minTotalEnergyLocal = max(-bandDepth + parallelEnergy / effectiveMass, workFunction + interpolator.getMinimumSampleEnergy() + parallelEnergy);

    assert(minTotalEnergyLocal <= maxTotalEnergy && "Total energy integration limits are wrong");

    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minTotalEnergyLocal, maxTotalEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
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
        double result = thisObj->parallelEnergyDistributionIntegrateTotal(parallelEnergy);
        if (thisObj->saveSpectra){
            thisObj->parallelEnergyDistribution.first.push_back(parallelEnergy);
            thisObj->parallelEnergyDistribution.second.push_back(result);
        }
        return result;
    };

    double(*integrationFunctionPointer)(double, void*) = integrationLambda; // convert the lambda into raw function pointer
    gsl_function gslIntegrationFunction = {integrationFunctionPointer, this};
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minParallelEnergy, maxParallelEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
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

    double maxParallelEnergyLocal;// = effectiveMass == 1. ? maxTotalEnergy - normalEnergy : min(maxTotalEnergy - normalEnergy, (normalEnergy + bandDepth) * effectiveMass / abs(1. - effectiveMass));
    double minParallelEnergyLocal = 0.;//  max(minTotalEnergy - normalEnergy, 0.);


    if (effectiveMass >= 1.){
        if (normalEnergy < -bandDepth && effectiveMass != 1.) minParallelEnergyLocal = -(normalEnergy + bandDepth) * effectiveMass / (effectiveMass - 1.);
        maxParallelEnergyLocal = maxTotalEnergy - normalEnergy;
    } else{
        maxParallelEnergyLocal = min((normalEnergy + bandDepth) * effectiveMass / (1. - effectiveMass), maxTotalEnergy - normalEnergy);
    }

    assert(minParallelEnergyLocal <= maxParallelEnergyLocal && "Parallel energy integration limits are wrong");

    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minParallelEnergyLocal, maxParallelEnergyLocal, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
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

    // double minNormalEnergy = min(xInitial, (1. - effectiveMass) * (maxTotalEnergy + bandDepth)); // capture the case that effectiveMass > 0
    // minNormalEnergy = max(minNormalEnergy, interpolator.getMinimumSampleEnergy() + workFunction);
    
    double result, error;
    gsl_integration_qag(&gslIntegrationFunction, minNormalEnergy, maxNormalEnergy, absoluteTolerance, relativeTolerance, maxAllowedSteps, 
                        GSL_INTEG_GAUSS31, externalIntegrationWorkSpace, &result, &error);
    
    if (saveNormalEnergyDistribution) 
        Utilities::sortCoordinates(normalEnergyDistribution.first, normalEnergyDistribution.second);
    
    
    return result;
}

}