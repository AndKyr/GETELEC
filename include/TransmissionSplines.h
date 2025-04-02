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
    void sampleUniform(double Emin, double Emax, int nPoints){
        if (nPoints < 2)
            throw std::invalid_argument("The number of points for the spline must be at least 2.");
    
        
        sampleEnergies = Utilities::linspace(Emin, Emax, nPoints);
        solutionSamples = vector<vector<double>>(3, vector<double>(nPoints, 0.));
        solutionDerivativeSamples = solutionSamples;
    
    
        for (size_t i = 0; i < sampleEnergies.size(); i++){
            solver.setEnergyAndInitialValues(sampleEnergies[i]);
            solver.solveNoSave();
            const auto& solutionVector = solver.getSolution();
            for (int j = 0; j < 3; j++){
                solutionSamples[j][i] = solutionVector[j];
                solutionDerivativeSamples[j][i] = solutionVector[j + 3] * CONSTANTS.kConstant;
            }
        }
    
        initializeMultiple(sampleEnergies, solutionSamples, solutionDerivativeSamples);    
    }


    void writeSplineSolution(string filename = "splineSolution.dat", int nPoints = 256){
        ofstream file(filename);

        vector<double> energyPoints = Utilities::linspace(getMinimumSampleEnergy(), getMaximumSampleEnergy(), nPoints);
        
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

    void writeSplineNodes(string filename = "splineNodes.dat"){
        ofstream file(filename);
        for (size_t i = 0; i < sampleEnergies.size(); i++){
            file << sampleEnergies[i] << " ";
            for (int j = 0; j < 3; j++){
                file << solutionSamples[j][i] << " " << solutionDerivativeSamples[j][i] << " ";
            }
            file << endl; 
        }
        file.close();
    }

    /**
     * @brief Samples the transmission coefficient at a given energy level.
     * @param energy Energy level (eV).
     */
    void sampleEnergyPoint(double energy){

        // ensure that the energy is within the limits of the solver
        if (energy < minimumValidEnergyForSolver){
            double newBarrierDepth = -energy + 5.;
            while(energy < minimumValidEnergyForSolver){
                solver.setXlimits(newBarrierDepth); // set the limits for the transmission solver.
                newBarrierDepth += 5.;
                minimumValidEnergyForSolver = solver.minimumValidEnergy();
                if (newBarrierDepth > 100.){
                    throw std::runtime_error("The solver's barrier depth was reduced below 100 eV without satisfying WKB. Something is wrong with the barrier shape");
                }
            }
        }

            
        solver.setEnergyAndInitialValues(energy);
        solver.solveNoSave();
        auto solutionVector = solver.getSolution();

        sampleEnergies.push_back(energy);
        for (int i = 0; i < 3; i++){
            solutionSamples[i].push_back(solutionVector[i]);
            solutionDerivativeSamples[i].push_back(solutionVector[i + 3] * CONSTANTS.kConstant);
        }
    }

    /**
     * @brief Samples the transmission coefficient with an initial very coarse sampling that is smartly chosen for electron emission.
     * Samples are taken at the Fermi level and at the top of the barrier. Additional samples are added if the transmission coefficient decays slowly.
     */
    void smartSampling(){

        assert((solver.getEnergyDerivativeLevel() > 0) && "The solver should be set with energyDerivatigveLevel>0 to get Hermite spline interpolation");

        // make sure that the solver has energy derivatives
        solver.setXlimits(workFunction + 5.); // set the limits for the transmission solver. 
        minimumValidEnergyForSolver = solver.minimumValidEnergy(); // get the minimum valid energy for the solver


        solver.getBarrier()->setBarrierTopFinder(true); // set barrier top finder
        sampleEnergyPoint(-workFunction); //sample the fermi level
        double topOfBarrier = solver.getBarrier()->getBarrierTop(); //get the barrier top
        solver.getBarrier()->setBarrierTopFinder(false); // reset the barrier top finder to off. No need any more


        //the extent by which the transmission coefficient decays below Fermi level
        double dLogD_dE_atFermi = -2. * solutionDerivativeSamples[2].back();
        double lowDecayLength = 1. / dLogD_dE_atFermi;
        double numberOfDecayLengths = -log(relativeTolerance);


        if (abs(lowDecayLength) > 1.e-2){ //it decays slow enough to be worth it
            double pointToSample = -workFunction - numberOfDecayLengths * lowDecayLength;
            sampleEnergyPoint(pointToSample);
        }

        double highDecayRate = 1./ kT - dLogD_dE_atFermi;
        double fermiToBarrier = topOfBarrier + workFunction;

        //If the spectra don't decay fast enough to have decayed before the barrier top (or don't decay at all), sample the top of barrier and a point beyond
        if (highDecayRate < numberOfDecayLengths / fermiToBarrier){ 
            sampleEnergyPoint(topOfBarrier);
            sampleEnergyPoint(topOfBarrier + numberOfDecayLengths* kT);
        }
        //If the spectra decay slow enough to be worth it, sample a above Ef. Force sample it we have only one sample point
        else if (highDecayRate < 100. ||  sampleEnergies.size() < 2){
            double pointToSample = -workFunction + numberOfDecayLengths / highDecayRate;
            sampleEnergyPoint(pointToSample);
        }

        sortAndInitialize();
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
    void sortAndInitialize(){

        // Lambda to compare the indices based on the sample energies
        // This lambda captures 'this' to access the member variables of the class
        auto comparisonLambda = [this](size_t i1, size_t i2) { return this->sampleEnergies[i1] < this->sampleEnergies[i2]; };

        vector<size_t> indices(sampleEnergies.size());
        iota(indices.begin(), indices.end(), 0); // Fill indices with 0, 1, ..., n-1
        // Sort the indices based on the sample energies
        sort(indices.begin(), indices.end(), comparisonLambda);
        
        auto savedSolutions = solutionSamples;
        auto savedDerivatives = solutionDerivativeSamples;
        auto savedEnergies = sampleEnergies;
        // Reorder the solution and derivative samples based on the sorted indices

        for (size_t j = 0; j < sampleEnergies.size(); j++){
            sampleEnergies[j] = savedEnergies[indices[j]];
            for (int i = 0; i < 3; i++){
                solutionSamples[i][j] = savedSolutions[i][indices[j]];
                solutionDerivativeSamples[i][j] = savedDerivatives[i][indices[j]];
            }
        }
    
        initializeMultiple(sampleEnergies, solutionSamples, solutionDerivativeSamples);
    }
};

} // namespace getelec

#endif /* TRANSMISSIONSPLINES_H_ */
