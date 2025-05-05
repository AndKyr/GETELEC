#include "Getelec.h"
#include <tbb/global_control.h>

using namespace std;

extern "C" {

// Wrapper to create a new Getelec object
getelec::Getelec* Getelec_new() {
    return new getelec::Getelec();
}

getelec::Getelec* Getelec_new_with_config(const char* configFileName, const char* barrierType, int numberOfThreads, int seed) {
    if (numberOfThreads > 0)
        static tbb::global_control tbbGlobalControl(tbb::global_control::max_allowed_parallelism, numberOfThreads);

    return new getelec::Getelec(configFileName, barrierType, NULL, seed);
}

// Wrapper to delete a Getelec object
void Getelec_delete(getelec::Getelec* obj) {
    delete obj;
}

void Getelec_setField(getelec::Getelec* obj, const double* fields, size_t size) {
    obj->setField(fields, size);
}

void Getelec_setRadius(getelec::Getelec* obj, const double* radii, size_t size) {
    obj->setRadius(radii, size);
}

void Getelec_setGamma(getelec::Getelec* obj, const double* gammas, size_t size) {
    obj->setGamma(gammas, size);
}


void Getelec_setkT(getelec::Getelec* obj, const double* kT, size_t size) {
    obj->setkT(kT, size);
}

void Getelec_setWorkFunction(getelec::Getelec* obj, const double* workFunction, size_t size) {
    obj->setWorkFunction(workFunction, size);
}

void Getelec_setEffectiveMass(getelec::Getelec* obj, const double* effectiveMass, size_t size) {
    obj->setEffectiveMass(effectiveMass, size);
}

void Getelec_setBandDepth(getelec::Getelec* obj, const double* bandDepth, size_t size) {
    obj->setBandDepth(bandDepth, size);
}

void Getelec_setRandomParameters(getelec::Getelec* obj, size_t numberOfParameters){
    obj->setRandomParameters(numberOfParameters);
}

// Wrapper to run the calculation
size_t Getelec_run(getelec::Getelec* obj, unsigned flags) {
    return obj->run(static_cast<getelec::CalculationFlags>(flags));
}

const double* Getelec_getField(getelec::Getelec* obj, size_t* size) {
    return obj->getFields(size);
}

const double* Getelec_getRadius(getelec::Getelec* obj, size_t* size) {
    return obj->getRadii(size);
}

const double* Getelec_getGamma(getelec::Getelec* obj, size_t* size) {
    return obj->getGammas(size);
}

const double* Getelec_getkT(getelec::Getelec* obj, size_t* size) {
    return obj->getkT(size);
}

const double* Getelec_getWorkFunction(getelec::Getelec* obj, size_t* size) {
    return obj->getWorkFunction(size);
}

const double* Getelec_getEffectiveMass(getelec::Getelec* obj, size_t* size) {
    return obj->getEffectiveMass(size);
}

const double* Getelec_getBandDepth(getelec::Getelec* obj, size_t* size) {
    return obj->getBandDepth(size);
}

// Wrapper to get current densities
const double* Getelec_getCurrentDensities(getelec::Getelec* obj, size_t* size) {
    return obj->getCurrentDensities(size);
}

// Wrapper to get Nottingham heats
const double* Getelec_getNottinghamHeats(getelec::Getelec* obj, size_t* size) {
    return obj->getNottinghamHeats(size);
}

// Wrapper to get spectra energies
const double* Getelec_getSpectraEnergies(getelec::Getelec* obj, size_t i, size_t* length, char spectraType) {
    return obj->getSpectraEnergies(i, length, spectraType);
}

// Wrapper to get spectra values
const double* Getelec_getSpectraValues(getelec::Getelec* obj, size_t i, size_t* length, char spectraType) {
    return obj->getSpectraValues(i, length, spectraType);
}

// Wrapper to get spectra derivatives
const double* Getelec_getSpectraDerivatives(getelec::Getelec* obj, size_t i, size_t* length) {
    return obj->getSpectraValues(i, length, 'D');
}

// Wrapper to calculate transmission coefficient for a single energy
double Getelec_calculateTransmissionProbability(getelec::Getelec* obj, double energy, double waveVector, size_t paramsIndex) {
    return obj->calculateTransmissionProbability(energy, waveVector, paramsIndex);
}

// Wrapper to calculate transmission coefficient for multiple energies
const double* Getelec_calculateTransmissionProbabilities(getelec::Getelec* obj, const double* energies, const double* waveVectors, size_t size, size_t paramsIndex) {
    return obj->calculateTransmissionProbabilities(vector<double>(energies, energies + size), vector<double>(waveVectors, waveVectors + size), paramsIndex).data();
}

void Getelec_getBarrierValues(getelec::Getelec* obj, const double* x, double* potential, size_t size, size_t paramsIndex) {
    obj->getBarrierValues(x, potential, size, paramsIndex);
}

// Wrapper to calculate transmission coefficient for many energies
const double* Getelec_interpolateTransmissionProbabilities(getelec::Getelec* obj, const double* energies, const double* waveVectors, size_t size, size_t paramsIndex) {
    return obj->interpolateTransmissionProbabilities(vector<double>(energies, energies + size), vector<double>(waveVectors, waveVectors + size), paramsIndex).data();
}

void Getelec_getBarrierIntegrationLimits(getelec::Getelec* obj, double* xInitial, double* xFinal, size_t paramsIndex) {  
    auto [xI, xF] = obj->getBarrierIntegrationLimits(paramsIndex);
    *xInitial = xI;
    *xFinal = xF;
}

unsigned Getelec_getCalculationStatus(getelec::Getelec* obj) {
    return static_cast<unsigned>(obj->getCalculationStatusFlags());
}


} // extern "C"