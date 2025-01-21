#ifndef GETELEC_H
#define GETELEC_H

#include <vector>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "BandEmitter.h"
#include "TunnelingFunction.h"
#include <vector>

class Getelec {
public:
    /**
     * @brief Construct a new Getelec object
     */
    Getelec() : threadLocalBarrier(), 
                threadLocalSolver([this] {return TransmissionSolver(&threadLocalBarrier.local());}),
                threadLocalEmitter([this] {return BandEmitter(threadLocalSolver.local());})
    {}

    ~Getelec(){}

    /**
     * @brief Set the field parameters of the calculation at a single value
     * @param field_ The electric field in V/nm.
     */
    void setField(double field_ = 5.) { threadLocalParams.local().field = field_; }

    /**
     * @brief Set the field parameters of the calculation at multiple values
     * @param fieldsVector_ The electric field in V/nm, multiple values to iterate over.
     */
    void setField(std::vector<double>* fieldsVector_) { fieldsVector = fieldsVector_; }


    /**
     * @brief Set the radius of the emitter at a single value
     * @param radius_ The radius of the emitter in nm.
     */
    void setRadius(double radius_ = 1.e4) { threadLocalParams.local().radius = radius_; }

    /**
     * @brief Set the radius of the emitter at multiple values
     * @param radiiVector_ The radius of the emitter in nm, multiple values to iterate over.
     */
    void setRadius(std::vector<double>* radiiVector_) { radiiVector = radiiVector_; }


    /**
     * @brief Set the gamma parameter of the general barrier model at a single value
     * @param gamma_ The gamma parameter of the general barrier model.
     */
    void setGamma(double gamma_ = 10.) { threadLocalParams.local().gamma = gamma_; }

    /**
     * @brief Set the gamma parameter of the general barrier model at multiple values
     * @param gammasVector_ The gamma parameter of the general barrier model, multiple values to iterate over.
     */
    void setGamma(std::vector<double>* gammasVector_) { gammasVector = gammasVector_; }


    /**
     * @brief Set the temperature at a single value
     * @param kT_ Temperature in eV.
     */
    void setkT(double kT_ = .025) { threadLocalParams.local().kT = kT_; }

    /**
     * @brief Set the temperature at multiple values
     * @param kTVector_ Temperature in eV, multiple values to iterate over.
     */
    void setkT(std::vector<double>* kTVector_) { kTVector = kTVector_; }


    /**
     * @brief Set the work function at a single value
     * @param workFunction_ Work function in eV.
     */
    void setWorkFunction(double workFunction_ = 4.5) { threadLocalParams.local().workFunction = workFunction_; }
    
    /**
     * @brief Set the work function at multiple values
     * @param workFunctionVector_ Work function in eV, multiple values to iterate over.
     */
    void setWorkFunction(std::vector<double>* workFunctionVector_) { workFunctionVector = workFunctionVector_; }


    /**
     * @brief Set the band depth at a single value
     * @param bandDepth_ Depth of the electronic band in eV.
     */
    void setBandDepth(double bandDepth_ = 10.) { threadLocalParams.local().bandDepth = bandDepth_; }
    
    /**
     * @brief Set the band depth at multiple values
     * @param bandDepthVector_ Depth of the electronic band in eV, multiple values to iterate over.
     */
    void setBandDepth(std::vector<double>* bandDepthVector_) { bandDepthVector = bandDepthVector_; } 


    /**
     * @brief Set the effective mass at a single value
     * @param effectiveMass_ Effective mass of the electron.
     */
    void setEffectiveMass(double effectiveMass_ = 1.) { threadLocalParams.local().effectiveMass = effectiveMass_; }
    
    /**
     * @brief Set the effective mass at multiple values
     * @param effectiveMassVector_ Effective mass of the electron, multiple values to iterate over.
     */
    void setEffectiveMass(vector<double>* effectiveMassVector_) { effectiveMassVector = effectiveMassVector_; }   

    /**
     * @brief Calculate for the i-th element of the array of inputs
     * @param i The index of the element to calculate
     * @param calculateSpectra If true, calculate the spectra
     */
    void runIteration(size_t i = 0, bool calculateSpectra = false);
    //  {
    //     setParamsForIteration(i);
        
    //     auto& params = threadLocalParams.local();
    //     auto& barrier = threadLocalBarrier.local();
    //     auto& emitter = threadLocalEmitter.local();

    //     barrier.setBarrierParameters(params.field, params.radius, params.gamma);
    //     emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth); 

    //     if (calculateSpectra)
    //         emitter.calculateCurrentDensityAndSpectra();
    //     else
    //         emitter.calculateCurrentDensityAndNottingham(); 
        
    //     currentDensityVector[i] = emitter.getCurrentDensity();
    //     nottinghamHeatVector[i] = emitter.getNottinghamHeat();
    //     if (calculateSpectra)
    //         spectra.push_back(emitter.getSpectra());

    // }

    /**
     * @brief Run the calculation
     * @param calculateSpectra If true, calculate the spectra
     */
    void run(bool calculateSpectra = false);
    // {
    //     currentDensityVector.resize(getMaxIterations());
    //     nottinghamHeatVector.resize(getMaxIterations());
    //     spectra.clear();
    //     tbb::parallel_for(size_t(0), size_t(getMaxIterations()), [this, calculateSpectra](size_t i) { runIteration(i, calculateSpectra);});
    // }

    /**
     * @brief Calculate the transmission coefficient for a specific energy
     * @param energy The energy level (eV)
     * @param paramsIndex The index of the parameter vector space to use (optional)
     * @return The transmission coefficient
     * @note This method is relevant for a single or a few calculations of the transmission coefficient. Don't use it for multiple calculations on the same barrier as it resets the barrier which might be slow. Use calculateTransmissionCoefficientForEnergies instead
     */
    double calculateTransmissionCoefficientForEnergy(double energy, size_t paramsIndex = numeric_limits<size_t>::max());
    // {
    //     setParamsForIteration(paramsIndex);   
    //     auto& params = threadLocalParams.local();
    //     auto& barrier = threadLocalBarrier.local();
    //     auto& emitter = threadLocalEmitter.local();

    //     barrier.setBarrierParameters(params.field, params.radius, params.gamma);
    //     emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, false); 
    //     emitter.setTransmissionSolver();
    //     return emitter.calculateTransmissionCoefficientForEnergy(energy);
    // }

    /**
     * @brief Calculate the transmission coefficient for multiple energies
     * @param energies The energy levels (eV)
     * @return The transmission coefficients
     * @note This method is relevant for multiple calculations of the transmission coefficient. It is faster than calculateTransmissionCoefficientForEnergy for multiple calculations on the same barrier. However, if you are iterating over many many energies, it might be better to use calculateTransmissionCoefficientForManyEnergies, which prepares the interpolator and then just interpolates.
     */
    vector<double> calculateTransmissionCoefficientForEnergies(const vector<double>& energies, size_t paramsIndex = numeric_limits<size_t>::max());
    // {

    //     setParamsForIteration(paramsIndex);   
    //     auto& params = threadLocalParams.local();
    //     auto& barrier = threadLocalBarrier.local();
    //     auto& emitter = threadLocalEmitter.local();

    //     barrier.setBarrierParameters(params.field, params.radius, params.gamma);
    //     emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, false); 
    //     emitter.setTransmissionSolver();

    //     tbb::concurrent_vector<double> transmissionCoefficients(energies.size());
        
    //     tbb::parallel_for(size_t(0), energies.size(), [&transmissionCoefficients, &energies, &emitter, this](size_t i) { 
    //         transmissionCoefficients[i] = emitter.calculateTransmissionCoefficientForEnergy(energies[i]);
    //     });
    //     return vector<double>(transmissionCoefficients.begin(), transmissionCoefficients.end());
    // }

        /**
     * @brief Calculate the transmission coefficient for multiple energies
     * @param energies The energy levels (eV)
     * @return The transmission coefficients
     * @note This method is relevant for multiple calculations of the transmission coefficient. It is faster than calculateTransmissionCoefficientForEnergy for multiple calculations on the same barrier. However, if you are iterating over many many energies, it might be better to use calculateTransmissionCoefficientForManyEnergies, which prepares the interpolator and then just interpolates.
     */
    vector<double> calculateTransmissionCoefficientForManyEnergies(const vector<double>& energies, size_t paramsIndex = numeric_limits<size_t>::max());
    // {

    //     setParamsForIteration(paramsIndex);   
    //     auto& params = threadLocalParams.local();
    //     auto& barrier = threadLocalBarrier.local();
    //     auto& emitter = threadLocalEmitter.local();

    //     barrier.setBarrierParameters(params.field, params.radius, params.gamma);
    //     emitter.setParameters(params.workFunction, params.kT, params.effectiveMass, params.bandDepth, true); 

    //     tbb::concurrent_vector<double> transmissionCoefficients(energies.size());
        
    //     tbb::parallel_for(size_t(0), energies.size(), [&transmissionCoefficients, &energies, &emitter, this](size_t i) { 
    //         transmissionCoefficients[i] = emitter.interpolateTransmissionCoefficientForEnergy(energies[i]);
    //     });
    //     return vector<double>(transmissionCoefficients.begin(), transmissionCoefficients.end());
    // }

    /**
     * @brief Get the current density at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The current density in A/nm^2
     */
    double getCurrentDensity(unsigned i = 0) const { return currentDensityVector[i]; }
    
    /**
     * @brief Get the Nottingham heat at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The Nottingham heat in W/nm^2
     */
    double getNottinghamHeat(unsigned i = 0) const { return nottinghamHeatVector[i]; }
    
    /**
     * @brief Get the spectra at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    pair<const vector<double>&, const vector<double>&> getSpectra(unsigned i = 0) const { return spectra[i]; }
    
    /**
     * @brief Get all the spectra
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    vector<pair<const vector<double>&, const vector<double>&>> getAllSpectra(unsigned i) const { return vector<pair<const vector<double>&, const vector<double>&>>(spectra.begin(), spectra.end()); }
    
    /**
     * @brief Get all the resulting current densities as a vector
     * @return The vector
     * @note This is a copy of the internal vector, so it is safe to modify it.
     */
    vector<double> getCurrentDensities() const { return vector<double>(currentDensityVector.begin(), currentDensityVector.end()); }
    
    /**
     * @brief Get all the resulting Nottingham heats as a vector
     * @return The vector
     * @note This is a copy of the internal vector, so it is safe to modify it.
     */
    vector<double> getNottinghamHeats() const { return vector<double>(nottinghamHeatVector.begin(), nottinghamHeatVector.end()); }

private:
    struct ParamsForIteration{
        double workFunction = 4.5; ///< Work function in eV.
        double kT = .025; ///< Temperature in eV.
        double effectiveMass = 1.; ///< Effective mass of the electron.
        double bandDepth = 10.; ///< The depth of the electronic band in eV.
        double field = 5.; ///< The electric field in V/nm.
        double radius = 1.e5; ///< The radius of the emitter in nm.
        double gamma = 10.; ///< The gamma parameter of the general barrier model.
    };

    const vector<double>* fieldsVector = NULL; ///< The electric field in V/nm, multiple values to iterate over.
    const vector<double>* radiiVector = NULL; ///< The radius of the emitter in nm, multiple values to iterate over.
    const vector<double>* gammasVector = NULL; ///< Gamma parameter of the general barrier model, multiple values to iterate over.
    const vector<double>* kTVector = NULL; ///< Temperature in eV, multiple values to iterate over.
    const vector<double>* workFunctionVector = NULL; //< Work function in eV, multiple values to iterate over.
    const vector<double>* bandDepthVector = NULL; //< Depth of the electronic band in eV, multiple values to iterate over.
    const vector<double>* effectiveMassVector = NULL; ///< Effective mass of the electron, multiple values to iterate over.

    double currentDensity = 0.; //< The current density (output) in A/nm^2.
    tbb::concurrent_vector<double> currentDensityVector; ///< The current density (output) in A/nm^2, multiple values to iterate over.

    double nottinghamHeat = 0.; //< The Nottingham heat (output) in W/nm^2.
    tbb::concurrent_vector<double> nottinghamHeatVector; ///< The Nottingham heat (output) in W/nm^2, multiple values to iterate over.

    tbb::concurrent_vector<pair<const vector<double>&, const vector<double>&>> spectra; /**< The spectra (output) in A/nm^2/eV, multiple values to iterate over. */

    tbb::enumerable_thread_specific<ModifiedSNBarrier> threadLocalBarrier; ///< Thread-local instances of ModifiedSNBarrier
    tbb::enumerable_thread_specific<TransmissionSolver> threadLocalSolver; ///< Thread-local instances of TransmissionSolver
    tbb::enumerable_thread_specific<BandEmitter> threadLocalEmitter; ///< Thread-local instances of BandEmitter
    tbb::enumerable_thread_specific<ParamsForIteration> threadLocalParams;


    /**
     * @brief Get the maximum number of iterations
     * @return The maximum number of iterations
     */
    int getMaxIterations();


    void setParamsForIteration(size_t i = numeric_limits<size_t>::max());
    // {
    //     ParamsForIteration& params = threadLocalParams.local();
        
    //     if (workFunctionVector && i < workFunctionVector->size())
    //         params.workFunction = (*workFunctionVector)[i];

    //     if (kTVector && i < kTVector->size())
    //         params.kT = (*kTVector)[i];

    //     if (effectiveMassVector && i < effectiveMassVector->size())
    //         params.effectiveMass = (*effectiveMassVector)[i];
        
    //     if (bandDepthVector && i < bandDepthVector->size())
    //         params.bandDepth = (*bandDepthVector)[i];

    //     if (fieldsVector && i < fieldsVector->size())
    //         params.field = (*fieldsVector)[i];
        
    //     if (radiiVector && i < radiiVector->size())
    //         params.radius = (*radiiVector)[i];
        
    //     if (gammasVector && i < gammasVector->size())
    //         params.gamma = (*gammasVector)[i];
        
    // }


};

extern "C" {

    // Wrapper to create a new Getelec object
    Getelec* Getelec_new() {
        return new Getelec();
    }

    // Wrapper to delete a Getelec object
    void Getelec_delete(Getelec* obj) {
        delete obj;
    }

    // Wrapper to set field values
    void Getelec_setField(Getelec* obj, const double* fields, size_t size) {
        std::vector<double> fieldVector(fields, fields + size);
        obj->setField(&fieldVector);
    }

    // Wrapper to set radius values
    void Getelec_setRadius(Getelec* obj, const double* radii, size_t size) {
        std::vector<double> radiusVector(radii, radii + size);
        obj->setRadius(&radiusVector);
    }

    // Wrapper to run the calculation
    void Getelec_run(Getelec* obj, bool calculateSpectra) {
        obj->run(calculateSpectra);
    }

    // Wrapper to get current densities
    const double* Getelec_getCurrentDensities(Getelec* obj, size_t* size) {
        static std::vector<double> currentDensities = obj->getCurrentDensities();
        *size = currentDensities.size();
        return currentDensities.data();
    }

} // extern "C"

#endif // GETELEC_H