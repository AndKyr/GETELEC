#ifndef GETELEC_H
#define GETELEC_H

#include <vector>
#include <tbb/tbb.h>
#include "BandEmitter.h"
#include "TunnelingFunction.h"
#include <vector>

class Getelec {
public:
    Getelec();
    ~Getelec();

    void setField(double field_ = 5.) { field = field_; }
    void setField(const vector<double>& fieldsVector_) { fieldsVector = fieldsVector_; }

    void setRadius(double radius_ = 1.e4) { radius = radius_; }
    void setRadius( const vector<double>& radiiVector_) { radiiVector = radiiVector_; }

    void setGamma(double gamma_ = 10.) { gamma = gamma_; }
    void setGamma(const vector<double>& gammasVector_) { gammasVector = gammasVector_; }

    void setkT(double kT_ = .025) { kT = kT_; }
    void setkT(const vector<double>& kTVector_) { kTVector = kTVector_; }

    void setWorkFunction(double workFunction_ = 4.5) { workFunction = workFunction_; }
    void setWorkFunction(const vector<double>& workFunctionVector_) { workFunctionVector = workFunctionVector_; }   

    void setBandDepth(double bandDepth_ = 10.) { bandDepth = bandDepth_; }
    void setBandDepth(const vector<double>& bandDepthVector_) { bandDepthVector = bandDepthVector_; }   

    void setEffectiveMass(double effectiveMass_ = 1.) { effectiveMass = effectiveMass_; }
    void setEffectiveMass(const vector<double>& effectiveMassVector_) { effectiveMassVector = effectiveMassVector_; }   

    /**
     * @brief Calculate 
     */
    void runIteration(unsigned i){
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
        

        
        threadLocalEmitter.local().setParameters(workFunctionVector[i], kTVector[i], effectiveMassVector[i], bandDepthVector[i]); 
        threadLocalEmitter.local().calculateCurrentDensityAndSpectra();
    }
private:
    double field; //< The electric field in V/nm.
    vector<double>& fieldsVector; //< The electric field in V/nm, multiple values to iterate over.

    double radius; //< The radius of the emitter in nm.
    vector<double>& radiiVector; //< The radius of the emitter in nm, multiple values to iterate over.

    double gamma; //< Gamma parameter of the general barrier model.
    vector<double>& gammasVector; //< Gamma parameter of the general barrier model, multiple values to iterate over.

    double kT; //< Temperature in eV. 
    vector<double>& kTVector; //< Temperature in eV, multiple values to iterate over.

    double workFunction; //< Work function in eV.
    vector<double>& workFunctionVector; //< Work function in eV, multiple values to iterate over.

    double bandDepth; //< Depth of the electronic band in eV.
    vector<double>& bandDepthVector; //< Depth of the electronic band in eV, multiple values to iterate over.

    double effectiveMass; //< Effective mass of the electron.
    vector<double>& effectiveMassVector; //< Effective mass of the electron, multiple values to iterate over.

        // Create a thread-local instance of ExampleClass
    tbb::enumerable_thread_specific<BandEmitter> threadLocalEmitter;
    BandEmitter emitter;
};

#endif // GETELEC_H