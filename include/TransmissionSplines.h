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
 * @class TransmissionSpline
 * @brief A class for interpolating and evaluating transmission coefficients and probabilities using cubic Hermite splines.
 * 
 * The TransmissionSpline class is designed to handle the interpolation of transmission coefficients
 * and probabilities for electron emission problems. It uses a cubic Hermite spline to ensure smooth
 * interpolation and provides methods for sampling, evaluating, and writing the spline data to files.
 * 
 * @details
 * This class is built on top of the CubicHermiteSpline base class and integrates with a TransmissionSolver
 * to calculate transmission coefficients. It supports various sampling strategies, including uniform
 * sampling and smart sampling tailored for electron emission problems. The class also provides methods
 * to estimate emission currents and refine the sampling to meet specified tolerances.
 * 
 * @note The class assumes that the energy levels are provided in electron volts (eV) and the wave vector
 * is given in units of 1/nm.
 * 
 * @note The default work function is set to 4.5 eV, and the default thermal energy is set to 0.025 eV.
 * 
 * @note The spline solution and nodes can be written to files for visualization or further analysis.
 * 
 * @see CubicHermiteSpline
 * @see TransmissionSolver
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
        return TransmissionSolver::transmissionCoefficient(waveVector, solutionVector);
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
     * @param index Index where to emplace the sample point (default is -1, which means append to the end).
     * @param bisect Whether the new point is to be bisected in the next round
     */
    void sampleEnergyPoint(double energy, size_t index = 0, bool bisect = true);
    
    /**
     * @brief Samples the transmission coefficient with an initial very coarse sampling that is smartly chosen for electron emission.
     * Samples are taken at the Fermi level and at the top of the barrier. Additional samples are added if the transmission coefficient decays slowly.
     */
    void smartSampling();

    double emissionCurrentEstimate(double energy, double waveVector = 12.){
        return Utilities::logFermiDiracFunction(energy + workFunction, kT) * getTransmissionProbability(energy, waveVector) * CONSTANTS.SommerfeldConstant * kT;
    }

    double getMinimumSampleEnergy() const {
        return sampleEnergies.front();
    }
    double getMaximumSampleEnergy() const {
        return sampleEnergies.back();
    }
    double getMinimumValidEnergyForSolver() const {
        return minimumValidEnergyForSolver;
    }

    /**
     * @brief Finds an estimate of the maximum emission current density of the normal energy distribution.
     * @param maxIterations The maximum number of minimization iterations allowed
     * @param rTol The relative tolerance of the minimum estimate
     */
    int findMaximumCurrentEstimate(int maxIterations = 10, double rTol = 0.1);

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
    vector<bool> bisectList; /**< List of booleans indicating if the interval following a sample point should be bisected */
    double mamximumCurrentEstimate = 0; /**< Estimation of the maximum current emitted (integrand of the NED) */
    double maximumCurrentPosition = 0; /**< Estimation of the position (in normal energy) of the maximum emission */

    /**
     * @brief Initializes the interpolator with the sampled data after sorting the samples.
     */
    void sortAndInitialize();

    /**
     * @brief Refines the sampling by bisecting all intervals that are set to be bisected by bisectList.
     * @return The number of new nodes added to the sampling list.
     * @note This method is called after the initial sampling to ensure that the spline is smooth and meets the specified tolerances.
     */
    int refineSampling(); 
};

} // namespace getelec

#endif /* TRANSMISSIONSPLINES_H_ */
