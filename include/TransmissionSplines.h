#ifndef TRANSMISSIONSPLINES_H_
#define TRANSMISSIONSPLINES_H_

#include "Utilities.h"
#include "TransmissionSolver.h"
#include "BSpline.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

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
class TransmissionSpline : public CubicHermiteSpline {
private:
    TransmissionSolver& solver; /**< Reference to a TransmissionSolver for transmission coefficient calculations. */
    double kT; /**< Thermal energy (eV). */
    double workFunction; /**< Work function of the material (eV). */
    double bandDepth; /**< Bottom of the band (eV). */
    double effectiveMass = 1.; /**< Effective mass of the material. */
    double minimumEnergy = 0.; /**< Minimum energy level for the interpolation. */
    double maximumEnergy = 0.; /**< Maximum energy level for the interpolation. */

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
    TransmissionSpline(TransmissionSolver& solver_, 
                            double workFunction_ = 4.5, 
                            double kT_ = .025, 
                            double bandDepth_ = 10.,
                            double effecitveMass_ = 1.,
                            double aTol = 1.e-12, 
                            double rTol = 1.e-5) 
                                : solver(solver_), workFunction(workFunction_), 
                                kT(kT_), bandDepth(bandDepth_), effectiveMass(effecitveMass_)
    {}

    // /**
    //  * @brief Evaluates the y-value for a given x-coordinate during interpolation.
    //  * @param x The x-coordinate, representing energy (eV).
    //  * @return The y-value, representing the logarithm of the transmission coefficient.
    //  */
    // double calculateYforX(double normalEnergy) override;

    // /**
    //  * @brief Calculates the error between the interpolated and calculated values.
    //  * @param x The x-coordinate, representing energy (eV).
    //  * @param yCalculated The calculated y-value.
    //  * @return The error between the calculated and interpolated values.
    //  */
    // double calculateError(double energy, double logD) override;

    // /**
    //  * @brief Computes the tolerance level for the interpolation process.
    //  * @param x The x-coordinate, representing energy (eV).
    //  * @param yValue The y-value, representing the logarithm of the transmission coefficient.
    //  * @return The tolerance level at the given x-coordinate.
    //  */
    // double calculateTolerance(double energy, double logD) override;
    
    // /** @brief set the temperature and work function parameters
    //  * @param kT_ Thermal energy (eV).
    //  * @param W Work function (eV).
    //  */
    void setParameters(double kT_, double W, double bandDepth_, double effectiveMass);

    gsl_complex getTransmissionCoefficient(double energy, double waveVector) const{
        vector<double> solutionVector = evaluateMultiple(energy);
        gsl_complex incidentWaveCoeffTimes2;
        GSL_SET_COMPLEX(&incidentWaveCoeffTimes2, waveVector + solutionVector[1], -solutionVector[0]);
        return gsl_complex_mul_real(gsl_complex_inverse(incidentWaveCoeffTimes2), 2. * waveVector * exp(-solutionVector[2]));
    }
    
    double getTransmissionProbability(double energy, double waveVector) const{
        gsl_complex transmissionCoeff = getTransmissionCoefficient(energy, waveVector);
        return gsl_complex_abs2(transmissionCoeff) / waveVector;
    }


    void sampleUniform(double Emin, double Emax, int nPoints){
        if (nPoints < 2)
            throw std::invalid_argument("The number of points for the spline must be at least 2.");
    
    
        vector<double> energyPoints = Utilities::linspace(Emin, Emax, nPoints);
        auto solutionValues = vector<vector<double>>(3, vector<double>(nPoints, 0.));
        auto solutionDerivatives = solutionValues;
    
    
        for (size_t i = 0; i < energyPoints.size(); i++){
            solver.setEnergyAndInitialValues(energyPoints[i]);
            solver.solveNoSave();
            const auto& solutionVector = solver.getSolution();
            for (int j = 0; j < 3; j++){
                solutionValues[j][i] = solutionVector[j];
                solutionDerivatives[j][i] = solutionVector[j + 3] * CONSTANTS.kConstant;
            }
        }
    
        initializeMultiple(energyPoints, solutionValues, solutionDerivatives);
        minimumEnergy = Emin;
        maximumEnergy = Emax;
    
    }

    void writeSplineSolution(string filename = "splineSolution.dat", int nPoints = 256){
        ofstream file(filename);

        auto energyPoints = Utilities::linspace(minimumEnergy, maximumEnergy, nPoints);
        for (size_t i = 0; i < energyPoints.size(); i++){
            solver.setEnergyAndInitialValues(energyPoints[i]);
            solver.solveNoSave();
            const auto& solutionVector = solver.getSolution();
            auto interpolatedValues = evaluateMultiple(energyPoints[i]);
            file << energyPoints[i] << " " << solutionVector[0] << " " << solutionVector[1] << " " << solutionVector[2] << " ";
            file << interpolatedValues[0] << " " << interpolatedValues[1] << " " << interpolatedValues[2] << endl;
        }
        file.close();
    }

};

} // namespace getelec

#endif /* TRANSMISSIONSPLINES_H_ */
