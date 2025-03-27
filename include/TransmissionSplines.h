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
                            double minimumEnergy_ = -10.,
                            double effecitveMass_ = 1.,
                            double aTol = 1.e-12, 
                            double rTol = 1.e-5) 
                                : solver(solver_), workFunction(workFunction_), 
                                kT(kT_), minSamplingEnergy(minimumEnergy_)
    {}

    void setParameters(double kT_, double W, double minimumEnergy_){
        kT = kT_;
        workFunction = W;
        minSamplingEnergy = minimumEnergy_;
    }

    /**
     * @brief Evaluates the transmission coefficient at a given energy and wave vector.
     * @param energy Energy level (eV).
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
     * @param energy Energy level (eV).
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

        auto energyPoints = Utilities::linspace(minSamplingEnergy, 0., nPoints);
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
        solver.setEnergyAndInitialValues(energy);
        solver.solveNoSave();
        auto solutionVector = solver.getSolution();
        auto solutionValues = vector<double>(3);
        auto solutionDerivatives = solutionValues;

        for (int i = 0; i < 3; i++){
            solutionValues[i] = solutionVector[i];
            solutionDerivatives[i] = solutionVector[i + 3] * CONSTANTS.kConstant;
        }

        solutionSamples.push_back(solutionValues);
        solutionDerivativeSamples.push_back(solutionDerivatives);
        sampleEnergies.push_back(energy);
    }

    /**
     * @brief Samples the transmission coefficient with an initial very coarse sampling that is smartly chosen for electron emission.
     * Samples are taken at the Fermi level and at the top of the barrier. Additional samples are added if the transmission coefficient decays slowly.
     */
    void smartSampling(){
        solver.getBarrier()->setBarrierTopFinder(true);
        sampleEnergyPoint(-workFunction);
        double topOfBarrier = solver.getBarrier()->getBarrierTop();
        solver.getBarrier()->setBarrierTopFinder(false);


        //the extent by which the transmission coefficient decays below Fermi level
        double dLogD_dE_atFermi = -2. * solutionDerivativeSamples.back()[2];
        double lowDecayLength = 1. / dLogD_dE_atFermi;

        if (abs(lowDecayLength) > 1.e-2){ //it decays slow enough to be worth it
            double pointToSample = max(-workFunction - 2*lowDecayLength, minSamplingEnergy);
            sampleEnergyPoint(pointToSample);
        }

        double highDecayLength = 1. / (1./ kT - dLogD_dE_atFermi);

        if (abs(highDecayLength) > 1.e-2 || sampleEnergies.size() < 2){
            double pointToSample = -workFunction + 2*highDecayLength;
            if (pointToSample < topOfBarrier){
                sampleEnergyPoint(pointToSample);
            } else{
                sampleEnergyPoint(topOfBarrier);
                sampleEnergyPoint(topOfBarrier + 5* kT);
            }
        }
        sortAndInitialize();
    }

private:
    TransmissionSolver& solver; /**< Reference to a TransmissionSolver for transmission coefficient calculations. */
    double kT; /**< Thermal energy (eV). */
    double workFunction; /**< Work function of the material (eV). */
    double minSamplingEnergy = -10; /**< Minimum energy level of interest for the interpolation. */

    vector<double> sampleEnergies; /**< List of sampled points */
    vector<vector<double>> solutionSamples; /**< List of sampled solution values */
    vector<vector<double>> solutionDerivativeSamples; /**< List of sampled solution derivative values */

    /**
     * @brief Initializes the interpolator with the sampled data after sorting the samples.
     */
    void sortAndInitialize(){
        // Create a vector of indices
        vector<size_t> indices(sampleEnergies.size());
        for (size_t i = 0; i < indices.size(); ++i)
            indices[i] = i;

        // Sort the indices based on the corresponding values in `data`
        std::sort(indices.begin(), indices.end(),
                [this](size_t i1, size_t i2) { return this->sampleEnergies[i1] < this->sampleEnergies[i2]; });
        
        vector<double> sortedEnergies(sampleEnergies.size());
        vector<vector<double>> sortedSolutions(sampleEnergies.size());
        vector<vector<double>> sortedDerivatives(sampleEnergies.size());

        for (size_t i = 0; i < indices.size(); i++){
            sortedEnergies[i] = sampleEnergies[indices[i]];
            sortedSolutions[i] = solutionSamples[indices[i]];
            sortedDerivatives[i] = solutionDerivativeSamples[indices[i]];
        }

        initializeMultiple(sortedEnergies, sortedSolutions, sortedDerivatives);
    }
};

} // namespace getelec

#endif /* TRANSMISSIONSPLINES_H_ */
