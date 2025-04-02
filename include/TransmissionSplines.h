#ifndef TRANSMISSIONSPLINES_H_
#define TRANSMISSIONSPLINES_H_

#include "Utilities.h"
#include "TransmissionSolver.h"
#include "BSpline.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <algorithm>

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

public:
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
                            double aTol = 1.e-12, 
                            double rTol = 1.e-5) 
                                : solver(solver_), 
                                workFunction(workFunction_), kT(kT_), 
                                relativeTolerance(rTol), absoluteTolerance(aTol)
    {}

    void setParameters(double kT_, double W){
        kT = kT_;
        workFunction = W;
    }

    /**
     * @brief Evaluates the transmission coefficient at a given energy and wave vector.
     * @param energy Energy level (eV) (from vacuum level).
     * @param waveVector Wave vector (1/nm).
     */
    gsl_complex getTransmissionCoefficient(double energy, double waveVector) const{
        vector<double> solutionVector = evaluateMultiple(energy);
        gsl_complex incidentWaveCoeffTimes2;
        GSL_SET_COMPLEX(&incidentWaveCoeffTimes2, waveVector + solutionVector[1], -solutionVector[0]);
        return gsl_complex_mul_real(gsl_complex_inverse(incidentWaveCoeffTimes2), 2. * waveVector * exp(-solutionVector[2]));
    }
    
    /**
     * @brief Evaluates the transmission probability at a given energy and wave vector.
     * @param energy Energy level (eV) (from vacuum level).
     * @param waveVector Wave vector (1/nm).
     */
    double getTransmissionProbability(double energy, double waveVector) const{
        gsl_complex transmissionCoeff = getTransmissionCoefficient(energy, waveVector);
        return gsl_complex_abs2(transmissionCoeff) / waveVector;
    }


    /**
     * @brief Uniformly samples the transmission coefficient between two energy levels.
     * @param Emin Minimum energy level (eV).
     * @param Emax Maximum energy level (eV).
     * @param nPoints Number of points to sample.
     */
    void sampleUniform(double Emin, double Emax, int nPoints);

    /**
     * @brief Writes the spline solution to a file for plotting.
     * @param filename Name of the output file.
     * @param nPoints Number of points to sample for the spline solution.
     * @note The default filename is "splineSolution.dat" and the default number of points is 256.
     */
    void writeSplineSolution(string filename = "splineSolution.dat", int nPoints = 256);

    /**
     * @brief Writes the spline nodes to a file for plotting.
     * @param filename Name of the output file.
     * @note The default filename is "splineNodes.dat".
     * @note The file contains the sampled energies, solution values, and derivative values.
     * @note The file is formatted as: energy, solution[0], derivative[0], solution[1], derivative[1], solution[2], derivative[2].
     */
    void writeSplineNodes(string filename = "splineNodes.dat");

    /**
     * @brief Samples the transmission coefficient at a given energy level.
     * @param energy Energy level (eV).
     */
    void sampleEnergyPoint(double energy);
    
    /**
     * @brief Samples the transmission coefficient with an initial very coarse sampling that is smartly chosen for electron emission.
     * Samples are taken at the Fermi level and at the top of the barrier. Additional samples are added if the transmission coefficient decays slowly.
     */
    void smartSampling();

    double getMinimumSampleEnergy() const {
        return sampleEnergies.front();
    }
    double getMaximumSampleEnergy() const {
        return sampleEnergies.back();
    }
    double getMinimumValidEnergyForSolver() const {
        return minimumValidEnergyForSolver;
    }

private:
    TransmissionSolver& solver; /**< Reference to a TransmissionSolver for transmission coefficient calculations. */
    double kT; /**< Thermal energy (eV). */
    double workFunction; /**< Work function of the material (eV). */
    double relativeTolerance; /**< Relative tolerance for the interpolation. */
    double absoluteTolerance; /**< Absolute tolerance for the interpolation. */
    double minimumValidEnergyForSolver; /**< Minimum energy that the solver at its current state produces accurate results */

    vector<double> sampleEnergies; /**< List of sampled points */
    vector<vector<double>> solutionSamples = vector<vector<double>>(3); /**< List of sampled solution values */
    vector<vector<double>> solutionDerivativeSamples = vector<vector<double>>(3); /**< List of sampled solution derivative values */

    /**
     * @brief Initializes the interpolator with the sampled data after sorting the samples.
     */
    void sortAndInitialize();
};

} // namespace getelec

#endif /* TRANSMISSIONSPLINES_H_ */
