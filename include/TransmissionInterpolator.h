#ifndef TRANSMISSIONINTERPOLATOR_H_
#define TRANSMISSIONINTERPOLATOR_H_




#include "TransmissionSolver.h"

class TransmissionInterpolator : public FunctionInterpolator{
private: 
    TransmissionSolver& solver;
    double kT;
    double workFunction;

public:

    /**The interpolator lives in the E_F = 0 convention */
    TransmissionInterpolator(TransmissionSolver& solver_, double workFunction_ = 4.5, double kT_ = .025, 
                                double aTol = 1.e-12, double rTol = 1.e-5) 
                            : FunctionInterpolator(aTol, rTol), solver(solver_),
                                workFunction(workFunction_), kT(kT_)
    {}

    double calculateYforX(double energy) override{
        return log(solver.calculateTransmissionCoefficientForEnergy(-workFunction + energy));
    }

    double calculateError(double energy, double logD) override{
        return Utilities::fermiDiracFunction(energy, kT) * abs(exp(logD) - exp(gsl_spline_eval(spline, energy, accelerator)));
    }

    double calculateTolerance(double energy, double logD) override {
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
    

    void setParameters(double kT_, double W){
        kT = kT_;
        workFunction = W;
    }

    double evaluate(double x) override{
        return exp(gsl_spline_eval(spline, x, accelerator));
    }
};

#endif // TRANSMISSIONINTERPOLATOR_H_