#include "Getelec.h"

namespace getelec{

void Getelec::runIteration(size_t i, CalculationFlags flags) {
    setParamsForIteration(i);
    
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth); 

    if (abs(params.effectiveMass - 1.) > config.bandEmitterParams.effectiveMassTolerance){
        runIterationEffectiveMassNonUnity(i, flags);
    } else {
        runIterationEffectiveMassUnity(i, flags);
    }

}

void Getelec::runIterationEffectiveMassUnity(size_t i, CalculationFlags flags){
    auto& emitter = threadLocalEmitter.local();

    if (flags == CalculationFlags::None){
        calculationStatusFlags = CalculationFlags::None;
        return;
    }
    
    if (flags == CalculationFlags::CurrentDensity){ // if we ask only for the current density
        currentDensityVector[i] = emitter.currentDensityIntegrateNormal(false);
        calculationStatusFlags = CalculationFlags::CurrentDensity;
        return;
    }

    if (flags & CalculationFlags::NormalEnergyDistribution){ // the normal energy distribution is requested
        currentDensityVector[i] = emitter.currentDensityIntegrateNormal(true);
        normalEnergyDistributions[i] = emitter.getNormalEnergyDistribution();
        calculationStatusFlags = CalculationFlags::NormalEnergyDistribution | CalculationFlags::CurrentDensity;
    }

    if (flags & CalculationFlags::ParallelEnergyDistribution){ // if the parallel energy distribution is asked
        currentDensityVector[i] = emitter.currentDensityIntegrateParallelTotal(true);
        parallelEnergyDistributions[i] = emitter.getParallelEnergyDistribution();
        calculationStatusFlags |= CalculationFlags::ParallelEnergyDistribution | CalculationFlags::CurrentDensity;
    }


    if (flags & CalculationFlags::TotalEnergyDistribution){ // the TED has been requested
        emitter.integrateTotalEnergyDistributionODEAndSaveSpectra();
        totalEnergyDistributions[i] = emitter.getTotalEnergyDistribution();
        totalEnergyDistributionsDerivatives[i] = emitter.getTotalEnergyDistributionDerivatives();
        calculationStatusFlags |= CalculationFlags::TotalEnergyDistribution | CalculationFlags::TotalEnergyDistributionDerivatives;
    } else if (flags & (CalculationFlags::CurrentDensity | CalculationFlags::NottinghamHeat)){ // the TED has not been asked, but current density or Nottingham have
        emitter.integrateTotalEnergyDistributionODE(true);
    }
    
    if (flags & CalculationFlags::NottinghamHeat){ // if Nottingham has been asked
        nottinghamHeatVector[i] = emitter.getNottinghamHeatODE() / CONSTANTS.electronCharge;
        calculationStatusFlags |= CalculationFlags::NottinghamHeat;
    }
    
    if (flags & CalculationFlags::CurrentDensity){
        currentDensityVector[i] = emitter.getCurrentDensityODE();
        calculationStatusFlags |= CalculationFlags::CurrentDensity;
    }
}

void Getelec::runIterationEffectiveMassNonUnity(size_t i, CalculationFlags flags){
    auto& emitter = threadLocalEmitter.local();

    if (flags == CalculationFlags::None){
        calculationStatusFlags = CalculationFlags::None;
        return;
    }

    if (flags & CalculationFlags::NormalEnergyDistribution){ // the normal energy distribution is requested
        currentDensityVector[i] = emitter.currentDensityIntegrateNormal(true);
        normalEnergyDistributions[i] = emitter.getNormalEnergyDistribution();
        calculationStatusFlags = CalculationFlags::NormalEnergyDistribution | CalculationFlags::CurrentDensity;
    }

    if (flags & CalculationFlags::ParallelEnergyDistribution){ // if the parallel energy distribution is asked
        currentDensityVector[i] = emitter.currentDensityIntegrateParallelTotal(true);
        parallelEnergyDistributions[i] = emitter.getParallelEnergyDistribution();
        calculationStatusFlags |= CalculationFlags::ParallelEnergyDistribution | CalculationFlags::CurrentDensity;
    }

    if (flags & CalculationFlags::TotalEnergyDistribution){ // the TED has been requested
        currentDensityVector[i] = emitter.currentDensityIntegrateTotalParallel(true);
        totalEnergyDistributions[i] = emitter.getTotalEnergyDistribution();
        calculationStatusFlags |= CalculationFlags::TotalEnergyDistribution | CalculationFlags::CurrentDensity;
    }
    
    if (flags & CalculationFlags::NottinghamHeat){ // if Nottingham has been asked
        nottinghamHeatVector[i] = emitter.nottinghamIntegrateTotalPrallel() / CONSTANTS.electronCharge;
        calculationStatusFlags |= CalculationFlags::NottinghamHeat;
    }
    
    if (flags & CalculationFlags::CurrentDensity && !(calculationStatusFlags & CalculationFlags::CurrentDensity)){ //Current density has been asked and not already calculated
        currentDensityVector[i] = emitter.currentDensityIntegrateTotalParallel(false);
        calculationStatusFlags |= CalculationFlags::CurrentDensity;
    }
}

void Getelec::getBarrierValues(const double* x, double* potential, size_t size, size_t paramsIndex) {
    setParamsForIteration(paramsIndex);
    auto& params = threadLocalParams.local(); 
    auto& barrier = threadLocalBarrier.local();
    barrier->setBarrierParameters(params.field, params.radius, params.gamma);

    for (size_t j = 0; j < size; ++j) {
        potential[j] = barrier->potentialFunction(x[j]);
    }
}

pair<double, double> Getelec::getBarrierIntegrationLimits(size_t paramIndex){
    setParamsForIteration(paramIndex);
    auto& params = threadLocalParams.local(); 
    auto& barrier = threadLocalBarrier.local();
    auto& solver = threadLocalSolver.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth);

    return {solver.getXFinal(), solver.getXInitial()};  
}


size_t Getelec::run(CalculationFlags flags) {
    // Get how many iterations to run
    size_t maxIterations = getMaxIterations();

    // Reset the flags to nothing calculated
    calculationStatusFlags = CalculationFlags::None;

    //Do the necessary array resizings
    if (flags & CalculationFlags::CurrentDensity)
        currentDensityVector.resize(maxIterations);
        
    if (flags & CalculationFlags::NottinghamHeat)
        nottinghamHeatVector.resize(maxIterations);

    if (flags & CalculationFlags::TotalEnergyDistribution)
        totalEnergyDistributions.resize(maxIterations);
    
    if (flags & CalculationFlags::NormalEnergyDistribution)
        normalEnergyDistributions.resize(maxIterations);
    
    if (flags & CalculationFlags::ParallelEnergyDistribution)
        parallelEnergyDistributions.resize(maxIterations);

    //run the iterations in parallel
    tbb::parallel_for(size_t(0), maxIterations, [this, flags](size_t i) { runIteration(i, flags); });
    return maxIterations;
}

double Getelec::calculateTransmissionCoefficientForEnergy(double energy, size_t paramsIndex) {
    setParamsForIteration(paramsIndex);   
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth); 
    return emitter.calculateTransmissionCoefficientForEnergy(energy);
}

std::vector<double> Getelec::calculateTransmissionCoefficientForEnergies(const std::vector<double>& energies, size_t paramsIndex) {
    setParamsForIteration(paramsIndex);   
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth); 

    tbb::concurrent_vector<double> transmissionCoefficients(energies.size());
    
    tbb::parallel_for(size_t(0), energies.size(), [&transmissionCoefficients, &energies, &emitter, this](size_t i) { 
        transmissionCoefficients[i] = emitter.calculateTransmissionCoefficientForEnergy(energies[i]);
    });
    return std::vector<double>(transmissionCoefficients.begin(), transmissionCoefficients.end());
}

std::vector<double> Getelec::calculateTransmissionCoefficientForManyEnergies(const std::vector<double>& energies, size_t paramsIndex) {
    setParamsForIteration(paramsIndex);   
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth); 

    tbb::concurrent_vector<double> transmissionCoefficients(energies.size());
    
    tbb::parallel_for(size_t(0), energies.size(), [&transmissionCoefficients, &energies, &emitter, this](size_t i) { 
        transmissionCoefficients[i] = emitter.interpolateTransmissionProbability(energies[i]);
    });
    return std::vector<double>(transmissionCoefficients.begin(), transmissionCoefficients.end());
}

size_t Getelec::getMaxIterations() {
    const std::vector<const std::vector<double>*> allInputVectors = {&fieldsVector, &radiiVector, &gammasVector, &kTVector, &workFunctionVector, &bandDepthVector, &effectiveMassVector};
    int maxSize = 0;
    for (auto* inputVector : allInputVectors)
        if (inputVector->size() > maxSize)
            maxSize = inputVector->size();

    for (auto* inputVector : allInputVectors)
        if (inputVector->size() != maxSize && inputVector->size() > 1)
            throw std::runtime_error("All input vectors must have the same size or size 1.");
    
    return maxSize;
}

void Getelec::setParamsForIteration(size_t i) {
    ParamsForIteration& params = threadLocalParams.local();
    
    if (i < workFunctionVector.size())
        params.workFunction = workFunctionVector[i];
    else if (workFunctionVector.size() == 1)
        params.workFunction = workFunctionVector[0];
    else
        throw std::runtime_error("Work function vector must have size 1 or the same size as the other input vectors.");

    if (i < kTVector.size())
        params.kT = kTVector[i];
    else if (kTVector.size() == 1)
        params.kT = kTVector[0];
    else
        throw std::runtime_error("kT vector must have size 1 or the same size as the other input vectors.");

    if (i < effectiveMassVector.size())
        params.effectiveMass = effectiveMassVector[i];
    else if (effectiveMassVector.size() == 1)
        params.effectiveMass = effectiveMassVector[0];
    else
        throw std::runtime_error("Effective mass vector must have size 1 or the same size as the other input vectors.");
    
    if (i < bandDepthVector.size())
        params.bandDepth = bandDepthVector[i];
    else if (bandDepthVector.size() == 1)
        params.bandDepth = bandDepthVector[0];
    else
        throw std::runtime_error("Band depth vector must have size 1 or the same size as the other input vectors.");

    if (i < fieldsVector.size())
        params.field = fieldsVector[i];
    else if (fieldsVector.size() == 1)
        params.field = fieldsVector[0];
    else
        throw std::runtime_error("Field vector must have size 1 or the same size as the other input vectors.");
    
    if (i < radiiVector.size())
        params.radius = radiiVector[i];
    else if (radiiVector.size() == 1)
        params.radius = radiiVector[0];
    else
        throw std::runtime_error("Radius vector must have size 1 or the same size as the other input vectors.");
    
    if (i < gammasVector.size())
        params.gamma = gammasVector[i];
    else if (gammasVector.size() == 1)
        params.gamma = gammasVector[0];
    else
        throw std::runtime_error("Gamma vector must have size 1 or the same size as the other input vectors.");
}

} // namespace getelec