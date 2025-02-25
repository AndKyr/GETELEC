#ifndef TRANSMISSIONINTERPOLATOR_H_
#define TRANSMISSIONINTERPOLATOR_H_

#include "Utilities.h"
#include "TransmissionSolver.h"

namespace getelec{

/**
 * @class TransmissionInterpolator
 * @brief Extends the FunctionInterpolator to refine and evaluate transmission coefficients.
 * 
 * This class leverages the functionality of FunctionInterpolator to calculate (sample) and 
 * interpolate transmission coefficients, which are calculated using a TransmissionSolver. 
 * It provides tools to refine the sampling 
 * and evaluate transmission-related data efficiently.
 * The interpolator lives in the E_F = 0 convention
 */
class TransmissionInterpolator : public FunctionInterpolator {
private:
    TransmissionSolver& solver; /**< Reference to a TransmissionSolver for transmission coefficient calculations. */
    double kT; /**< Thermal energy (eV). */
    double workFunction; /**< Work function of the material (eV). */

public:

    /** */
    /**
     * @brief Constructs a TransmissionInterpolator.
     * @param transmissionSolver Reference to a TransmissionSolver instance to calculate the sampled values.
     * @param workFunction Work function of the material (eV).
     * @param kT Thermal energy (eV).
     * @param aTol Absolute tolerance for interpolation.
     * @param rTol Relative tolerance for interpolation.
     */
    TransmissionInterpolator(TransmissionSolver& solver_, 
                            double workFunction_ = 4.5, 
                            double kT_ = .025, 
                            double aTol = 1.e-12, 
                            double rTol = 1.e-5) 
                                : FunctionInterpolator(aTol, rTol), 
                                solver(solver_),
                                workFunction(workFunction_), 
                                kT(kT_)
    {}

    /**
     * @brief Evaluates the y-value for a given x-coordinate during interpolation.
     * @param x The x-coordinate, representing energy (eV).
     * @return The y-value, representing the logarithm of the transmission coefficient.
     */
    double calculateYforX(double energy) override{
        return log(solver.calculateTransmissionCoefficientForEnergy(-workFunction + energy));
    }

    /**
     * @brief Calculates the error between the interpolated and calculated values.
     * @param x The x-coordinate, representing energy (eV).
     * @param yCalculated The calculated y-value.
     * @return The error between the calculated and interpolated values.
     */
    double calculateError(double energy, double logD) override{
        return Utilities::fermiDiracFunction(energy, kT) * abs(exp(logD) - exp(gsl_spline_eval(spline, energy, accelerator)));
    }

    /**
     * @brief Computes the tolerance level for the interpolation process.
     * @param x The x-coordinate, representing energy (eV).
     * @param yValue The y-value, representing the logarithm of the transmission coefficient.
     * @return The tolerance level at the given x-coordinate.
     */
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
    
    /** @brief set the temperature and work function parameters
     * @param kT_ Thermal energy (eV).
     * @param W Work function (eV).
     */
    void setParameters(double kT_, double W){
        kT = kT_;
        workFunction = W;
    }

    /** @brief Interpolates the sampled function at position @param x */
    double evaluate(double x) override{
        return exp(gsl_spline_eval(spline, x, accelerator));
    }
};

} // namespace getelec

#endif /* TRANSMISSIONINTERPOLATOR_H_ */
