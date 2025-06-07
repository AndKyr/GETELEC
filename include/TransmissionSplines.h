#ifndef TRANSMISSIONSPLINES_H_
#define TRANSMISSIONSPLINES_H_

#include "Utilities.h"
#include "TransmissionSolver.h"
#include "BSpline.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_min.h>

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
    {
        minimizer = gsl_min_fminimizer_alloc(type); /**< Allocate GSL minimizer for finding the max of the expected emission current */
    }

    ~TransmissionSpline(){
        if (minimizer) gsl_min_fminimizer_free(minimizer);
        minimizer = nullptr;
    }

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
        vector<double> solutionVector = evaluateSolution(energy);
        assert(all_of(solutionVector.begin(), solutionVector.end(), [](double x) { return isfinite(x); }));
        return TransmissionSolver::transmissionCoefficient(waveVector, solutionVector);
    }
    
    /**
     * @brief Evaluates the transmission probability at a given energy and wave vector.
     * @param energy Energy level (eV) (from vacuum level).
     * @param waveVector Wave vector (1/nm).
     */
    double getTransmissionProbability(double energy, double waveVector) const{
        vector<double> solutionVector = evaluateSolution(energy);
        return TransmissionSolver::transmissionProbability(waveVector, solutionVector);
    }

    /**
     * @brief Evaluates the solution vector at a given energy
     * @param energy The energy to evaluate at (eV)
     * @note If the energy is out of the interval where sampling nodes exist, a linear extrapolation is used
     */
    vector<double> evaluateSolution(double energy) const;

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
    int writeSplineSolution(string filename = "splineSolution.dat", int nPoints = 256, bool writeFullSolution = false) const;

    /**
     * @brief Writes the spline nodes to a file for plotting.
     * @param filename Name of the output file.
     * @note The default filename is "splineNodes.dat".
     * @note The file contains the sampled energies, solution values, and derivative values.
     * @note The file is formatted as: energy, solution[0], derivative[0], solution[1], derivative[1], solution[2], derivative[2].
     */
    int writeSplineNodes(string filename = "splineNodes.dat") const;
    
    /**
     * @brief Samples the transmission coefficient with an initial very coarse sampling that is smartly chosen for electron emission.
     * @details Samples are taken at the Fermi level and at the top of the barrier. Additional samples are added if the transmission coefficient decays slowly.
     */
    void smartInitialSampling(double bandDepth = 10, double effectiveMass = 1.);

    /**
     * @brief Refines the sampling by bisecting until the required tolerance is met.
     * TODO:  Cases where maxRefineSteps is reached should be handled separately. Consider designating a region where the interpolator fails and reverting to solver rather than the interpolator
     * when that energy is asked.
     */
    void refineSamplingToTolerance(int maxRefineSteps = 11);

    /**
     * @brief Estimates the emission normal Energy distribution at a given energy level.
     * @param energy Energy level (eV).
     * @param waveVector Wave vector (1/nm).
     * @return The normal energy distribution estimate (A/nm^2 /eV).
     * @note this function is coincides exactly with the normal energy distribution for effectiveMass = 1
     */
    double normalEnergyDistributionEstimate(double energy, double waveVector = 12.) const{
        return Utilities::logFermiDiracFunction(energy + workFunction, kT) * getTransmissionProbability(energy, waveVector) * CONSTANTS.SommerfeldConstant * kT;
    }

    double getMinimumSampleEnergy() const { return sampleEnergies.front(); }
    double getMaximumSampleEnergy() const { return sampleEnergies.back(); }

    /**
     * @brief Finds an estimate of the maximum emission current density of the normal energy distribution.
     * @param maxIterations The maximum number of minimization iterations allowed
     * @param rTol The relative tolerance of the minimum estimate
     */
    int findMaximumCurrentEstimate(int maxIterations = 15, double rTol = 0.1);

    double getmaximumCurrentEstimate(){ return mamximumCurrentEstimate; }

private:
    TransmissionSolver& solver; /**< Reference to a TransmissionSolver for transmission coefficient calculations. */
    double kT; /**< Thermal energy (eV). */
    double workFunction; /**< Work function of the material (eV). */
    double relativeTolerance; /**< Relative tolerance for the interpolation. */
    double absoluteTolerance; /**< Absolute tolerance for the interpolation. */

    vector<double> sampleEnergies; /**< List of sampled points */
    vector<vector<double>> solutionSamples = vector<vector<double>>(3); /**< List of sampled solution values */
    vector<vector<double>> solutionDerivativeSamples = vector<vector<double>>(3); /**< List of sampled solution derivative values */
    vector<bool> bisectList; /**< List of booleans indicating if the interval following a sample point should be bisected */
    double mamximumCurrentEstimate = 0; /**< Estimation of the maximum current emitted (integrand of the NED) */
    double maximumCurrentPosition = 0; /**< Estimation of the position (in normal energy) of the maximum emission */
    bool isInitialized = false; /**< Flag indicating if the interpolator is initialized */

    const double maxAllowedSamplingDistance = 2.; ///< The maximum allowed distance between energy sampling nodes (eV)
    const double maxAllowedSolutionError = 0.1; ///< The maximum allowed error that all solution components are allowed to have (pure number)

    const gsl_min_fminimizer_type* type = gsl_min_fminimizer_brent; /**< GSL object (minimizer type) for finding the max of the expected emission current */
    gsl_min_fminimizer* minimizer = nullptr; /**< GSL object (minimizer) for finding the max of the expected emission current */

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

    /**
     * @brief Runs the initial sampling for the extreme case of super long barriers (effectively only thermionic emission allowed)
     */
    void smartSamplingForLongBarrier();

    /**
     * @brief Performs the part of the initial sampling related to energies below the Fermi energy
     * @param initialSampleEnergy The energy point where the initial sample was taken
     * @return noTunnelingRegime (Whether the barrier is very long and the transmission probability decades extremely fast below the barrier top)
     */
    bool initialSamplingLowEnergies(double initialSampleEnergy);

    

    /**
     * @brief Samples the transmission coefficient at a given energy level.
     * @param energy Energy level (eV).
     * @param index Index where to emplace the sample point (default is -1, which means append to the end).
     * @param bisect Whether the new point is to be bisected in the next round
     */
    void sampleEnergyPoint(double energy, size_t index = 0, bool bisect = true, bool longBarrier = false);

    /**
     * @brief Calculates and sets the sample values and derivatives for the spline at a given index.
     * @param index Index of the sample point to calculate.
     * @note This method is used to calculate the solution and derivative values for the spline at the given index spline node in the sampleEnergies.
     */
    void calculateAndSetSampleValues(size_t index);

    /**
     * @brief Resets the spline data and clears the sample lists.
     * @note This method is used to clear the data before a new calculation or sampling.
     * @note It clears the sample energies, solution samples, and derivative samples.
     * @note It also clears the bisect list.
     * @note This method is called before a new calculation to ensure that the spline data is fresh and does not contain any old data.
     */
    void reset(){
        sampleEnergies.clear();
        for (auto& sample : solutionSamples)
            sample.clear();
        for (auto& sample : solutionDerivativeSamples)
            sample.clear();

        bisectList.clear();
        isInitialized = false;
    }

    /**
     * @brief Checks if all emission estimates are below the absolute tolerance.
     * @return True if all estimates are below the tolerance, false otherwise.
     * @note This method is used to check if the current is too small to be worth calculating.
     */
    bool checkAllBelowTolerance(){
        for (double energy : sampleEnergies){
            double emissionEstimate = normalEnergyDistributionEstimate(energy);
            if (emissionEstimate > absoluteTolerance) 
                return false;
        }
        return true;
    }

    /**
     * @brief Checks the sanity of a given solution vector by comparing it to the value that the interpolator predicts.
     * @param solutionVector The solution vector to check.
     * @param energy The energy level (eV) at which to interpolate to get the expected solution estimate.
     * @param absoluteTolerance The absolute tolerance for the check (if the distance of the given solution from the interpolated value is below this, the solution is considered sane).
     * @return True if the solution is sane, false otherwise.
     * @note This method is used mostly in debug mode to ensure that a solution calculated by the solver is continuous as a function of energy. In other words it checks if the solution curve solution(energy) is smooth.
     */
    bool checkSolutionSanity(const vector<double>& solutionVector, double energy, double absoluteTolerance = 0.01);
};

} // namespace getelec

#endif /* TRANSMISSIONSPLINES_H_ */
