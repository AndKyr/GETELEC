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
    void setField(std::vector<double>& fieldsVector_) { fieldsVector.swap(fieldsVector_); }

    void setRadius(double radius_ = 1.e4) { radius = radius_; }
    void setRadius(std::vector<double>& radiiVector_) { radiiVector.swap(radiiVector_); }

    void setGamma(double gamma_ = 10.) { gamma = gamma_; }
    void setGamma(std::vector<double>& gammasVector_) { gammasVector.swap(gammasVector_); }

    void setkT(double kT_ = .025) { kT = kT_; }
    void setkT(std::vector<double>& kTVector_) { kTVector.swap(kTVector_); }

    void setWorkFunction(double workFunction_ = 4.5) { workFunction = workFunction_; }
    void setWorkFunction(std::vector<double>& workFunctionVector_) { workFunctionVector.swap(workFunctionVector_); }

    void setBandDepth(double bandDepth_ = 10.) { bandDepth = bandDepth_; }
    void setBandDepth(std::vector<double>& bandDepthVector_) { bandDepthVector.swap(bandDepthVector_); } 

    void setEffectiveMass(double effectiveMass_ = 1.) { effectiveMass = effectiveMass_; }
    void setEffectiveMass(vector<double>& effectiveMassVector_) { effectiveMassVector.swap(effectiveMassVector_); }   

    /**
     * @brief Calculate for a single case
     */
    void runIteration(unsigned i = 0, bool calculateSpectra = false) {
        double workFunctionTmp = workFunction;
        if (i < workFunctionVector.size())
            workFunctionTmp = workFunctionVector[i];

        double kTTmp = kT;
        if (i < kTVector.size())
            kTTmp = kTVector[i];

        double effectiveMassTmp = effectiveMass;
        if (i < effectiveMassVector.size())
            effectiveMassTmp = effectiveMassVector[i];
        
        double bandDepthTmp = bandDepth;
        if (i < bandDepthVector.size())
            bandDepthTmp = bandDepthVector[i];

        double fieldTmp = field;
        if (i < fieldsVector.size())
            fieldTmp = fieldsVector[i];
        
        double radiusTmp = radius;
        if (i < radiiVector.size())
            radiusTmp = radiiVector[i];
        
        double gammaTmp = gamma;
        if (i < gammasVector.size())
            gammaTmp = gammasVector[i];
        

        threadLocalBarrier.local().setBarrierParameters(fieldTmp, radiusTmp, gammaTmp);
        threadLocalEmitter.local().setParameters(workFunctionTmp, kTTmp, effectiveMassTmp, bandDepthTmp); 

        if (calculateSpectra)
            threadLocalEmitter.local().calculateCurrentDensityAndSpectra();
        else
            threadLocalEmitter.local().calculateCurrentDensityAndNottingham(); 

    }

    void run(){
        tbb::parallel_for(size_t(0), size_t(getMaxIterations()), [this](size_t i) { runIteration(i);});
    }
private:
    double field = 5.; //< The electric field in V/nm.
    vector<double> fieldsVector; //< The electric field in V/nm, multiple values to iterate over.

    double radius = 1.e5; //< The radius of the emitter in nm.
    vector<double> radiiVector; //< The radius of the emitter in nm, multiple values to iterate over.

    double gamma = 10.; //< Gamma parameter of the general barrier model.
    vector<double> gammasVector; //< Gamma parameter of the general barrier model, multiple values to iterate over.

    double kT = .025; //< Temperature in eV. 
    vector<double> kTVector; //< Temperature in eV, multiple values to iterate over.

    double workFunction = 4.5; //< Work function in eV.
    vector<double> workFunctionVector; //< Work function in eV, multiple values to iterate over.

    double bandDepth = 10.; //< Depth of the electronic band in eV.
    vector<double> bandDepthVector; //< Depth of the electronic band in eV, multiple values to iterate over.

    double effectiveMass = 1.; //< Effective mass of the electron.
    vector<double> effectiveMassVector; //< Effective mass of the electron, multiple values to iterate over.

    tbb::enumerable_thread_specific<ModifiedSNBarrier> threadLocalBarrier; //< Thread-local instances of ModifiedSNBarrier
    tbb::enumerable_thread_specific<TransmissionSolver> threadLocalSolver; //< Thread-local instances of TransmissionSolver
    tbb::enumerable_thread_specific<BandEmitter> threadLocalEmitter; //< Thread-local instances of BandEmitter

    int getMaxIterations(){
        return max(fieldsVector.size(), max(radiiVector.size(), max(gammasVector.size(), max(kTVector.size(), max(workFunctionVector.size(), max(bandDepthVector.size(), effectiveMassVector.size()))))));
    }

};

#endif // GETELEC_H