#include "TransmissionInterpolator.h"

namespace getelec{


double TransmissionInterpolator::calculateYforX(double normalEnergy){
    double kineticEnergy = normalEnergy + bandDepth;
    if(kineticEnergy <= 0) kineticEnergy = 0.1;
    double waveVector = sqrt(kineticEnergy) * CONSTANTS.sqrt2mOverHbar;
    return log(solver.calculateTransmissionProbability(-workFunction + normalEnergy, waveVector));
}

double TransmissionInterpolator::calculateError(double energy, double logD){
    return Utilities::fermiDiracFunction(energy, kT) * abs(exp(logD) - exp(gsl_spline_eval(spline, energy, accelerator)));
}


double TransmissionInterpolator::calculateTolerance(double energy, double logD) {
    double maxEmissionEstimate = 0.;
    for (auto& v : samplingList){
        double emissionEstimate = exp(v.y) * Utilities::logFermiDiracFunction(v.x, kT);
        if (emissionEstimate > maxEmissionEstimate)
            maxEmissionEstimate = emissionEstimate;
    }

    double emissionEstimate = exp(logD) * Utilities::logFermiDiracFunction(energy, kT);

    if (energy <= 0 || logD > -2.) // if we are out of the sensitive region
        return absoluteTolerance + (maxEmissionEstimate) * relativeTolerance;
    else
        return absoluteTolerance + emissionEstimate * relativeTolerance;
}

void TransmissionInterpolator::setParameters(double kT_, double W, double bandDepth_, double effecitveMass_){
    kT = kT_;
    workFunction = W;
    bandDepth = bandDepth_;
    effectiveMass = effecitveMass_;
}


}
