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
    void setField(double field_ = 5.) { field = field_; }

    /**
     * @brief Set the field parameters of the calculation at multiple values
     * @param fieldsVector_ The electric field in V/nm, multiple values to iterate over.
     */
    void setField(std::vector<double>* fieldsVector_) { fieldsVector = fieldsVector_; }


    /**
     * @brief Set the radius of the emitter at a single value
     * @param radius_ The radius of the emitter in nm.
     */
    void setRadius(double radius_ = 1.e4) { radius = radius_; }

    /**
     * @brief Set the radius of the emitter at multiple values
     * @param radiiVector_ The radius of the emitter in nm, multiple values to iterate over.
     */
    void setRadius(std::vector<double>* radiiVector_) { radiiVector = radiiVector_; }


    /**
     * @brief Set the gamma parameter of the general barrier model at a single value
     * @param gamma_ The gamma parameter of the general barrier model.
     */
    void setGamma(double gamma_ = 10.) { gamma = gamma_; }

    /**
     * @brief Set the gamma parameter of the general barrier model at multiple values
     * @param gammasVector_ The gamma parameter of the general barrier model, multiple values to iterate over.
     */
    void setGamma(std::vector<double>* gammasVector_) { gammasVector = gammasVector_; }


    /**
     * @brief Set the temperature at a single value
     * @param kT_ Temperature in eV.
     */
    void setkT(double kT_ = .025) { kT = kT_; }

    /**
     * @brief Set the temperature at multiple values
     * @param kTVector_ Temperature in eV, multiple values to iterate over.
     */
    void setkT(std::vector<double>* kTVector_) { kTVector = kTVector_; }


    /**
     * @brief Set the work function at a single value
     * @param workFunction_ Work function in eV.
     */
    void setWorkFunction(double workFunction_ = 4.5) { workFunction = workFunction_; }
    
    /**
     * @brief Set the work function at multiple values
     * @param workFunctionVector_ Work function in eV, multiple values to iterate over.
     */
    void setWorkFunction(std::vector<double>* workFunctionVector_) { workFunctionVector = workFunctionVector_; }


    /**
     * @brief Set the band depth at a single value
     * @param bandDepth_ Depth of the electronic band in eV.
     */
    void setBandDepth(double bandDepth_ = 10.) { bandDepth = bandDepth_; }
    
    /**
     * @brief Set the band depth at multiple values
     * @param bandDepthVector_ Depth of the electronic band in eV, multiple values to iterate over.
     */
    void setBandDepth(std::vector<double>* bandDepthVector_) { bandDepthVector = bandDepthVector_; } 


    /**
     * @brief Set the effective mass at a single value
     * @param effectiveMass_ Effective mass of the electron.
     */
    void setEffectiveMass(double effectiveMass_ = 1.) { effectiveMass = effectiveMass_; }
    
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
    void runIteration(unsigned i = 0, bool calculateSpectra = false) {
        double workFunctionTmp = workFunction;
        if (workFunctionVector && i < workFunctionVector->size())
            workFunctionTmp = (*workFunctionVector)[i];

        double kTTmp = kT;
        if (kTVector && i < kTVector->size())
            kTTmp = (*kTVector)[i];

        double effectiveMassTmp = effectiveMass;
        if (effectiveMassVector && i < effectiveMassVector->size())
            effectiveMassTmp = (*effectiveMassVector)[i];
        
        double bandDepthTmp = bandDepth;
        if (bandDepthVector && i < bandDepthVector->size())
            bandDepthTmp = (*bandDepthVector)[i];

        double fieldTmp = field;
        if (fieldsVector && i < fieldsVector->size())
            fieldTmp = (*fieldsVector)[i];
        
        double radiusTmp = radius;
        if (radiiVector && i < radiiVector->size())
            radiusTmp = (*radiiVector)[i];
        
        double gammaTmp = gamma;
        if (gammasVector && i < gammasVector->size())
            gammaTmp = (*gammasVector)[i];
        

        threadLocalBarrier.local().setBarrierParameters(fieldTmp, radiusTmp, gammaTmp);
        threadLocalEmitter.local().setParameters(workFunctionTmp, kTTmp, effectiveMassTmp, bandDepthTmp); 

        if (calculateSpectra)
            threadLocalEmitter.local().calculateCurrentDensityAndSpectra();
        else
            threadLocalEmitter.local().calculateCurrentDensityAndNottingham(); 
        
        currentDensityVector[i] = threadLocalEmitter.local().getCurrentDensity();
        nottinghamHeatVector[i] = threadLocalEmitter.local().getNottinghamHeat();
        if (calculateSpectra)
            spectra.push_back(threadLocalEmitter.local().getSpectra());

    }

    /**
     * @brief Run the calculation
     * @param calculateSpectra If true, calculate the spectra
     */
    void run(bool calculateSpectra = false){
        currentDensityVector.resize(getMaxIterations());
        nottinghamHeatVector.resize(getMaxIterations());
        spectra.clear();
        tbb::parallel_for(size_t(0), size_t(getMaxIterations()), [this, calculateSpectra](size_t i) { runIteration(i, calculateSpectra);});
    }


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
    const vector<double> getNottinghamHeats() const { return vector<double>(nottinghamHeatVector.begin(), nottinghamHeatVector.end()); }

private:
    double field = 5.; //< The electric field in V/nm.
    const vector<double>* fieldsVector = NULL; //< The electric field in V/nm, multiple values to iterate over.

    double radius = 1.e5; //< The radius of the emitter in nm.
    const vector<double>* radiiVector = NULL; //< The radius of the emitter in nm, multiple values to iterate over.

    double gamma = 10.; //< Gamma parameter of the general barrier model.
    const vector<double>* gammasVector = NULL; //< Gamma parameter of the general barrier model, multiple values to iterate over.

    double kT = .025; //< Temperature in eV. 
    const vector<double>* kTVector = NULL; //< Temperature in eV, multiple values to iterate over.

    double workFunction = 4.5; //< Work function in eV.
    const vector<double>* workFunctionVector = NULL; //< Work function in eV, multiple values to iterate over.

    double bandDepth = 10.; //< Depth of the electronic band in eV.
    const vector<double>* bandDepthVector = NULL; //< Depth of the electronic band in eV, multiple values to iterate over.

    double effectiveMass = 1.; //< Effective mass of the electron.
    const vector<double>* effectiveMassVector = NULL; //< Effective mass of the electron, multiple values to iterate over.

    double currentDensity = 0.; //< The current density (output) in A/nm^2.
    tbb::concurrent_vector<double> currentDensityVector; //< The current density (output) in A/nm^2, multiple values to iterate over.

    double nottinghamHeat = 0.; //< The Nottingham heat (output) in W/nm^2.
    tbb::concurrent_vector<double> nottinghamHeatVector; //< The Nottingham heat (output) in W/nm^2, multiple values to iterate over.

    tbb::concurrent_vector<pair<const vector<double>&, const vector<double>&>> spectra; //< The spectra (output) in A/nm^2/eV, multiple values to iterate over.

    tbb::enumerable_thread_specific<ModifiedSNBarrier> threadLocalBarrier; //< Thread-local instances of ModifiedSNBarrier
    tbb::enumerable_thread_specific<TransmissionSolver> threadLocalSolver; //< Thread-local instances of TransmissionSolver
    tbb::enumerable_thread_specific<BandEmitter> threadLocalEmitter; //< Thread-local instances of BandEmitter


    /**
     * @brief Get the maximum number of iterations
     * @return The maximum number of iterations
     */
    int getMaxIterations(){
        const vector<const vector<double>*> allInputVectors = {fieldsVector, radiiVector, gammasVector, kTVector, workFunctionVector, bandDepthVector, effectiveMassVector};
        int maxSize = 0;
        for (auto inputVector : allInputVectors)
            if (inputVector && inputVector->size() > maxSize)
                maxSize = inputVector->size();
        return maxSize;
    }

};

#endif // GETELEC_H