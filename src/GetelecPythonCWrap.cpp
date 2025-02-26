#include "Getelec.h"
using namespace std;

#include <fstream>
#include <chrono>

extern "C" {

    // Wrapper to create a new Getelec object
    getelec::Getelec* Getelec_new() {
        return new getelec::Getelec();
    }

    getelec::Getelec* Getelec_new_with_config(const char* configFileName, const char* barrierType) {
        return new getelec::Getelec(configFileName, barrierType);
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

    // Wrapper to run the calculation
    size_t Getelec_run(getelec::Getelec* obj, bool calculateSpectra) {
        return obj->run(calculateSpectra);
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
    const double* Getelec_getSpectraEnergies(getelec::Getelec* obj, size_t i, size_t* length) {
        return obj->getSpectraEnergies(i, length);
    }

    // Wrapper to get spectra values
    const double* Getelec_getSpectraValues(getelec::Getelec* obj, size_t i, size_t* length) {
        return obj->getSpectraValues(i, length);
    }

    // Wrapper to get spectra derivatives
    const double* Getelec_getSpectraDerivatives(getelec::Getelec* obj, size_t i, size_t* length) {
        return obj->getSpectraDerivatives(i, length);
    }
    
    // Wrapper to calculate transmission coefficient for a single energy
    double Getelec_calculateTransmissionCoefficientForEnergy(getelec::Getelec* obj, double energy, size_t paramsIndex) {
        return obj->calculateTransmissionCoefficientForEnergy(energy, paramsIndex);
    }

    // Wrapper to calculate transmission coefficient for multiple energies
    const double* Getelec_calculateTransmissionCoefficientForEnergies(getelec::Getelec* obj, const double* energies, size_t size, size_t paramsIndex) {
        return obj->calculateTransmissionCoefficientForEnergies(vector<double>(energies, energies + size), paramsIndex).data();
    }

    void Getelec_getBarrierValues(getelec::Getelec* obj, const double* x, double* potential, size_t size, size_t paramsIndex) {
        obj->getBarrierValues(x, potential, size, paramsIndex);
    }

    // Wrapper to calculate transmission coefficient for many energies
    const double* Getelec_calculateTransmissionCoefficientForManyEnergies(getelec::Getelec* obj, const double* energies, size_t size, size_t paramsIndex) {
        return obj->calculateTransmissionCoefficientForManyEnergies(vector<double>(energies, energies + size), paramsIndex).data();
    }

    void Getelec_getBarrierIntegrationLimits(getelec::Getelec* obj, double* xInitial, double* xFinal, size_t paramsIndex) {  
        auto [xI, xF] = obj->getBarrierIntegrationLimits(paramsIndex);
        *xInitial = xI;
        *xFinal = xF;
    }


    string error;
    // Initialization function for COMSOL wrapper
    int init(const char *str) {

        ofstream file;
        file.open("/home/kyritsakis/comsolExtLibLog.txt", ios::app);
        if (!file.is_open()) {
            error = "error opening log file";
            return 0;
        }
        string locStr = str;

        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        file << "Initialized with string: " << locStr << ", at time = " << std::ctime(&now_time) << endl;
        file.close();

        return 1;
    }
     
    // Error message function for COMSOL wrapper
    const char* getLastError() {
    
        ofstream file;
        file.open("/home/kyritsakis/comsolExtLibLog.txt", ios::app);
        if (!file.is_open()) {
            error = "error opening log file";
            return error.c_str();
        }
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        string locStr = error;

        file << "Called getLastError with error string: " << locStr << ", at time = " << std::ctime(&now_time) << endl;
        file.close();

        return error.c_str();
    }
     
    int eval(const char *func,
                        int nArgs,
                        const double **inReal,
                        const double **inImag,
                        int blockSize,
                        double *outReal,
                        double *outImag) 
    {
        ofstream file;
        file.open("/home/kyritsakis/comsolExtLibLog.txt", ios::app);
        if (!file.is_open()) {
            error = "error opening log file";
            return 0;
        }
        string functionStr = func;

        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        file << "Called function with string: " << functionStr << " and nArgs= " << nArgs << ", at time = " << std::ctime(&now_time) << endl;
        file.close();

        getelec::Getelec* obj = new getelec::Getelec();

        if (nArgs > 7 || nArgs < 0) {
            error = "Invalid number of arguments. Expected is between 0 and 7";
            return 0;
        }

        if (nArgs >= 0)
            obj->setField(inReal[0], blockSize);
        if (nArgs >= 1)
            obj->setRadius(inReal[1], blockSize);
        if (nArgs >= 2)
            obj->setGamma(inReal[2], blockSize);
        if (nArgs >= 3)
            obj->setkT(inReal[3], blockSize);
        if (nArgs >= 4)
            obj->setWorkFunction(inReal[4], blockSize);
        if (nArgs >= 5)
            obj->setEffectiveMass(inReal[5], blockSize);
        if (nArgs >= 6)
            obj->setBandDepth(inReal[6], blockSize);
        
        obj->run(false);
        
        const double* outData;
        size_t outSize;
        if (functionStr == "currentDensity")
            outData = obj->getCurrentDensities(&outSize);
        else if (functionStr == "nottinghamHeat")
            outData = obj->getNottinghamHeats(&outSize);
        else {
            error = "Unknown function";
            return 0;
        }

        if (outSize != blockSize) {
            error = "Output size does not match input size";
            return 0;
        }

        for (size_t i = 0; i < blockSize; i++) {
            outReal[i] = outData[i];
        }
        return 1;
    }

} // extern "C"