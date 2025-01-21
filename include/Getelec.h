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
    Getelec() : threadLocalBarrier(), 
                threadLocalSolver([this] {return TransmissionSolver(&threadLocalBarrier.local());}),
                threadLocalEmitter([this] {return BandEmitter(threadLocalSolver.local());})
    {}

    ~Getelec(){}

    void setField(double field_ = 5.) { field = field_; }
    void setField(std::vector<double>* fieldsVector_) { fieldsVector = fieldsVector_; }

    void setRadius(double radius_ = 1.e4) { radius = radius_; }
    void setRadius(std::vector<double>* radiiVector_) { radiiVector = radiiVector_; }

    void setGamma(double gamma_ = 10.) { gamma = gamma_; }
    void setGamma(std::vector<double>* gammasVector_) { gammasVector = gammasVector_; }

    void setkT(double kT_ = .025) { kT = kT_; }
    void setkT(std::vector<double>* kTVector_) { kTVector = kTVector_; }

    void setWorkFunction(double workFunction_ = 4.5) { workFunction = workFunction_; }
    void setWorkFunction(std::vector<double>* workFunctionVector_) { workFunctionVector = workFunctionVector_; }

    void setBandDepth(double bandDepth_ = 10.) { bandDepth = bandDepth_; }
    void setBandDepth(std::vector<double>* bandDepthVector_) { bandDepthVector = bandDepthVector_; } 

    void setEffectiveMass(double effectiveMass_ = 1.) { effectiveMass = effectiveMass_; }
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

    double getCurrentDensity(unsigned i = 0) const { return currentDensityVector[i]; }
    double getNottinghamHeat(unsigned i = 0) const { return nottinghamHeatVector[i]; }
    pair<const vector<double>&, const vector<double>&> getSpectra(unsigned i = 0) const { return spectra[i]; }
    vector<pair<const vector<double>&, const vector<double>&>> getAllSpectra(unsigned i) const { return vector<pair<const vector<double>&, const vector<double>&>>(spectra.begin(), spectra.end()); }
    const vector<double> getCurrentDensities() const { return vector<double>(currentDensityVector.begin(), currentDensityVector.end()); }
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