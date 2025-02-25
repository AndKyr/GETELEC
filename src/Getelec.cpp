#include "Getelec.h"

namespace getelec{

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