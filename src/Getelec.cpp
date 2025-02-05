#include "Getelec.h"

void Getelec::runIteration(size_t i, bool calculateSpectra) {
    setParamsForIteration(i);
    
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth); 

    if (calculateSpectra)
        emitter.calculateCurrentDensityAndSpectra();
    else
        emitter.calculateCurrentDensityAndNottingham(); 
    
    currentDensityVector[i] = emitter.getCurrentDensity();
    nottinghamHeatVector[i] = emitter.getNottinghamHeat();
    if (calculateSpectra)
        spectra[i] = emitter.getSpectra();
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
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, false);
    emitter.setTransmissionSolver();

    return {solver.getXFinal(), solver.getXInitial()};  
}


size_t Getelec::run(bool calculateSpectra) {
    size_t maxIterations = getMaxIterations();
    currentDensityVector.resize(maxIterations);
    nottinghamHeatVector.resize(maxIterations);
    if (calculateSpectra)
        spectra.resize(maxIterations);
    tbb::parallel_for(size_t(0), maxIterations, [this, calculateSpectra](size_t i) { runIteration(i, calculateSpectra); });
    return maxIterations;
}

double Getelec::calculateTransmissionCoefficientForEnergy(double energy, size_t paramsIndex) {
    setParamsForIteration(paramsIndex);   
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, false); 
    emitter.setTransmissionSolver();
    return emitter.calculateTransmissionCoefficientForEnergy(energy);
}

std::vector<double> Getelec::calculateTransmissionCoefficientForEnergies(const std::vector<double>& energies, size_t paramsIndex) {
    setParamsForIteration(paramsIndex);   
    auto& params = threadLocalParams.local();
    auto& barrier = threadLocalBarrier.local();
    auto& emitter = threadLocalEmitter.local();

    barrier->setBarrierParameters(params.field, params.radius, params.gamma);
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, false); 
    emitter.setTransmissionSolver();

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
    emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, true); 

    tbb::concurrent_vector<double> transmissionCoefficients(energies.size());
    
    tbb::parallel_for(size_t(0), energies.size(), [&transmissionCoefficients, &energies, &emitter, this](size_t i) { 
        transmissionCoefficients[i] = emitter.interpolateTransmissionCoefficientForEnergy(energies[i]);
    });
    return std::vector<double>(transmissionCoefficients.begin(), transmissionCoefficients.end());
}

size_t Getelec::getMaxIterations() {
    const std::vector<const std::vector<double>*> allInputVectors = {&fieldsVector, &radiiVector, &gammasVector, &kTVector, &workFunctionVector, &bandDepthVector, &effectiveMassVector};
    int maxSize = 0;
    for (auto inputVector : allInputVectors)
        if (inputVector && inputVector->size() > maxSize)
            maxSize = inputVector->size();
    return maxSize;
}

void Getelec::setParamsForIteration(size_t i) {
    ParamsForIteration& params = threadLocalParams.local();
    
    if (i < workFunctionVector.size())
        params.workFunction = workFunctionVector[i];

    if (i < kTVector.size())
        params.kT = kTVector[i];

    if (i < effectiveMassVector.size())
        params.effectiveMass = effectiveMassVector[i];
    
    if (i < bandDepthVector.size())
        params.bandDepth = bandDepthVector[i];

    if (i < fieldsVector.size())
        params.field = fieldsVector[i];
    
    if (i < radiiVector.size())
        params.radius = radiiVector[i];
    
    if (i < gammasVector.size())
        params.gamma = gammasVector[i];
}