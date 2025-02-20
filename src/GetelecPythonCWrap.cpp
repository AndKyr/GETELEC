#include "Getelec.h"
#include <jni.h>


extern "C" {

    // Wrapper to create a new Getelec object
    Getelec* Getelec_new() {
        return new Getelec();
    }

    Getelec* Getelec_new_with_config(const char* configFileName, const char* barrierType) {
        return new Getelec(configFileName, barrierType);
    }

    // Wrapper to delete a Getelec object
    void Getelec_delete(Getelec* obj) {
        delete obj;
    }

    void Getelec_setField(Getelec* obj, const double* fields, size_t size) {
        obj->setField(fields, size);
    }

    void Getelec_setRadius(Getelec* obj, const double* radii, size_t size) {
        obj->setRadius(radii, size);
    }

    void Getelec_setGamma(Getelec* obj, const double* gammas, size_t size) {
        obj->setGamma(gammas, size);
    }


    void Getelec_setkT(Getelec* obj, const double* kT, size_t size) {
        obj->setkT(kT, size);
    }

    void Getelec_setWorkFunction(Getelec* obj, const double* workFunction, size_t size) {
        obj->setWorkFunction(workFunction, size);
    }

    void Getelec_setEffectiveMass(Getelec* obj, const double* effectiveMass, size_t size) {
        obj->setEffectiveMass(effectiveMass, size);
    }
    
    void Getelec_setBandDepth(Getelec* obj, const double* bandDepth, size_t size) {
        obj->setBandDepth(bandDepth, size);
    }

    // Wrapper to run the calculation
    size_t Getelec_run(Getelec* obj, bool calculateSpectra) {
        return obj->run(calculateSpectra);
    }

    // Wrapper to get current densities
    const double* Getelec_getCurrentDensities(Getelec* obj, size_t* size) {
        return obj->getCurrentDensities(size);
    }

    // Wrapper to get Nottingham heats
    const double* Getelec_getNottinghamHeats(Getelec* obj, size_t* size) {
        return obj->getNottinghamHeats(size);
    }

    // Wrapper to get spectra energies
    const double* Getelec_getSpectraEnergies(Getelec* obj, size_t i, size_t* length) {
        return obj->getSpectraEnergies(i, length);
    }

    // Wrapper to get spectra values
    const double* Getelec_getSpectraValues(Getelec* obj, size_t i, size_t* length) {
        return obj->getSpectraValues(i, length);
    }

    // Wrapper to get spectra derivatives
    const double* Getelec_getSpectraDerivatives(Getelec* obj, size_t i, size_t* length) {
        return obj->getSpectraDerivatives(i, length);
    }
    
    // Wrapper to calculate transmission coefficient for a single energy
    double Getelec_calculateTransmissionCoefficientForEnergy(Getelec* obj, double energy, size_t paramsIndex) {
        return obj->calculateTransmissionCoefficientForEnergy(energy, paramsIndex);
    }

    // Wrapper to calculate transmission coefficient for multiple energies
    const double* Getelec_calculateTransmissionCoefficientForEnergies(Getelec* obj, const double* energies, size_t size, size_t paramsIndex) {
        return obj->calculateTransmissionCoefficientForEnergies(vector<double>(energies, energies + size), paramsIndex).data();
    }

    void Getelec_getBarrierValues(Getelec* obj, const double* x, double* potential, size_t size, size_t paramsIndex) {
        obj->getBarrierValues(x, potential, size, paramsIndex);
    }

    // Wrapper to calculate transmission coefficient for many energies
    const double* Getelec_calculateTransmissionCoefficientForManyEnergies(Getelec* obj, const double* energies, size_t size, size_t paramsIndex) {
        return obj->calculateTransmissionCoefficientForManyEnergies(vector<double>(energies, energies + size), paramsIndex).data();
    }

    void Getelec_getBarrierIntegrationLimits(Getelec* obj, double* xInitial, double* xFinal, size_t paramsIndex) {  
        auto [xI, xF] = obj->getBarrierIntegrationLimits(paramsIndex);
        *xInitial = xI;
        *xFinal = xF;
    }

} // extern "C"