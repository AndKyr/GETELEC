#ifndef BANDEMITTER_H_
#define BANDEMITTER_H_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include "TransmissionSolver.h"
// #include "TransmissionInterpolator.h"
#include "ConfigGetelec.h"
#include <vector>
#include <string>
#include <fstream>
#include "TransmissionSplines.h"

namespace getelec{
using namespace std;

/**
 * TODO: all documentation and comments are chatgpt generated. Please review them.
 */

/**
 * @class BandEmitter
 * @brief Simulates the electron emission process from an isotropic parabolic band, i.e. E = |k|^2 / 2 m*
 */
class BandEmitter : public ODESolver {
private:
    /** @brief The work function of the material in eV. */
    double workFunction = 4.5;

    /** @brief Thermal energy in eV, proportional to temperature (kT). */
    double kT = 0.025;

    /** @brief Effective mass of electrons in the material, relative to the free electron mass. */
    double effectiveMass = 1.;

    /** @brief Depth of the electronic band from the Fermi level in eV. */
    double bandDepth = 10.;

    /**
     * @brief Energy distributions saved as pairs of vectors. The first vector stores the energy abcissae (eV)
     * and the second vector the electron emission distribution (A / nm^2 / eV)
     */
    pair<vector<double>, vector<double>> totalEnergyDistribution; ///< Total energy distribution
    pair<vector<double>, vector<double>> parallelEnergyDistribution; ///< Parallel energy distribution
    pair<vector<double>, vector<double>> normalEnergyDistribution; ///< Normal energy distribution

    /** @brief Stores calculated spectral value derivatives in A/nm^2 / eV / eV */
    vector<double> totalEnergyDistributionDerivatives;

    /** @brief A reference to the transmission solver used for calculating coefficients. */
    TransmissionSolver& transmissionSolver;

    /** @brief An interpolator for efficient evaluation of transmission coefficients. */
    TransmissionSpline interpolator;

    /** @brief A pointer to a random number generator to use for testing purposes */
    mt19937* generator = NULL;

    CubicHermiteSpline totalEnergyDistributionSpline; /**< Spline for the energy spectra. */

    Config::BandEmitterParams configParams; /**< Configuration parameters for the BandEmitter class. */
  
    gsl_integration_workspace* integrationWorkspace = NULL;  /**< GSL Workspace for numerical integration using GSL. */
    gsl_integration_workspace* externalIntegrationWorkSpace = nullptr; /**< Second GSL Workspace for numerical integration using GSL. 
                                                                                This space is used for double integral,
                                                                                 when the internal and external integrals need to use independent workspaces */

    double helperEnergy; /**< Helper variable storing the totalEnergy variable for evaluating the double integral. */

    bool saveSpectra; /**< Helper flag that determines whether the parallel energy distribution will be saved upon double integration */

    /**
     * @brief Calculates the g'(E) (see Andreas' notes) for a given (normal) energy.
     * @param energy The normal energy of the electron.
     * @return The calculated g'(E) value.
     * @note For effectiveMass=1, g'(E) simplifies into the transmission probability, i.e. g'(E) = D(E)
     * @note For effectiveMass far away from 1, g'(E) remains a complex integral and double integration should be used.
     */
    double gPrimeFunction(double energy);

    /**
     * @brief Defines the system of differential equations to solve.
     * @param energy The independent variable (energy).
     * @param y The dependent variables.
     * @param f The derivatives of the dependent variables.
     * @param params Additional system parameters.
     * @return GSL_SUCCESS on success.
     */
    static int differentialSystem(double energy, const double y[], double f[], void *params);


    /** @brief Defines the Jacobian matrix for the differential system (optional, currently unused). */
    static int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

public:
    /**
     * @brief Constructs a BandEmitter object.
     * @param solver A reference to the transmission solver.
     * @param workFun The work function in eV.
     * @param kT_ Thermal energy in eV.
     * @param effMass Effective mass of the electron.
     * @param bandDepth_ The depth of the electronic band in eV.
     * @param rtol Relative tolerance for ODE solving.
     * @param atol Absolute tolerance for ODE solving.
     * @param maxSteps Maximum allowed steps for ODE solving.
     * @param minSteps Minimum allowed steps for ODE solving.
     * @param stepExpectedForInitialStep Initial step expectation for ODE solving.
     */
    BandEmitter(TransmissionSolver& solver,
                Config::BandEmitterParams config = Config().bandEmitterParams,
                double workFun = 4.5,
                double kT_ = .025,
                double effMass = 1.,
                double bandDepth_ = 7.
                )
        : ODESolver(vector<double>(3, 0.0), differentialSystem, 3, {0., 1.}, 
            config.relativeTolerance, config.absoluteTolerance, "rkck", 
            config.maxSteps, config.minSteps, config.stepExpectedForInitialStep, NULL, this), 
          transmissionSolver(solver), 
          interpolator(solver, workFun, kT_, config.absoluteTolerance, config.relativeTolerance),
          configParams(config) 
    {
        setParameters(workFun, kT_, effMass, bandDepth_);
    }

    /**
     * @brief Destructor for the BandEmitter object.
     */
    ~BandEmitter() {
        if (integrationWorkspace){
            gsl_integration_workspace_free(integrationWorkspace);
            integrationWorkspace = NULL;
        }

        if(externalIntegrationWorkSpace){
            gsl_integration_workspace_free(externalIntegrationWorkSpace);
            externalIntegrationWorkSpace = NULL;
        }
    }

    /**
     * @brief Sets the parameters for the band emitter simulation.
     * @param workFunction_ The work function in eV.
     * @param kT_ Thermal energy in eV.
     * @param effectiveMass_ The effective mass of the electron (fraction of the electron rest mass, i.e. m_n/m_e).
     * @param bandDepth_ The depth of the electronic band in eV.
     * @param doQuadrature Whether quadrature integration is to be used. If false, the integration workspace is not allocated.
     * @note Quadrature is not needed if only the ODE will be used, which is valid for effectiveMass~1 and only current density, Nottingham heat, and TED spectra are to be evaluated
     */
    void setParameters(double workFunction_ = 4.5, double kT_ = 0.025, double effectiveMass_ = 1., double bandDepth_ = 7., bool doQuadrature = true);

    /**
     * @brief Sets parameters to (meaningful) random values
     * @note This is a convenience function for testing purposes.
     */
    void setRandomParameters(){
        double bandDepth = Utilities::getUniformRandomDouble(0., 10., *generator);
        double workFunction = Utilities::getUniformRandomDouble(0., 6., *generator);
        double effMass = Utilities::getUniformRandomDouble(0.8, 2., *generator);
        double kT_ = Utilities::getUniformRandomDouble(0., 2., *generator);
        setParameters(workFunction, kT_, effMass, bandDepth);
    }

    /**
     * @brief Solves the ODE system to calculate current density, Nottingham heat and energy spectra.
     * @param convergenceTolerance The relative tolerance for convergence of the integration towards high currents. If <=0 it reverts to the object's relativeTolrance
     * @param makeSpectralSpline Whether to create a spline for the spectra.
     * @return Status of the calculation (e.g., GSL_SUCCESS).
     */
    int integrateTotalEnergyDistributionODEAndSaveSpectra(double convergenceTolerance = 0., bool makeSpectralSpline = true);

    /**
     * @brief Solves the ODE system to calculate current density and Nottingham heat.
     * @param convergenceTolerance The tolerance for convergence.
     * @return Status of the calculation (e.g., GSL_SUCCESS).
     * 
     * @note This method runs the same as integrateTotalEnergyDistributionODEAndSaveSpectra but a bit faster because it does not save the spectra.
     */
    int integrateTotalEnergyDistributionODE(double convergenceTolerance = 1.e-5);

    /**
     * @brief Calculates the current density using numerical integration.
     * @return The calculated current density.
     */
    double currentDensityIntegrateNormal(bool saveNormalEnergyDistribution = false);

    /**
     * @brief Writes data for plotting the band emitter's behavior.
     * @param filename The name of the output file.
     */
    void writePlottingData(string filename = "bandEmitterPlotting.dat");

    /**
     * @brief Function for calculating the normal energy distribution (strictly correct only for effectiveMass = 1.).
     * @param energy The energy of the electron (eV).
     * @return The calculated distribution value (if multiplied by kT * Sommerfeld constant, in A / nm^2 / eV).
     */
    double normalEnergyDistributionForEnergy(double energy){
        return  gPrimeFunction(energy) * Utilities::logFermiDiracFunction(energy, kT);
    }

    /**
     * @brief Gets the temperature parameter kT (eV) 
     */
    double getkT() const { return kT; }

    /**
     * @brief Gets the effective mass of the band (fraction of electron mass)
     */
    double getEffectiveMass() const { return effectiveMass; }

    /**
     * @brief Gets the band depth (eV)
     */
    double getBandDepth() const { return bandDepth; }

    /**
     * @brief Gets the work function (eV)
     */
    double getWorkFunction() const { return workFunction; }

    /**
     * @brief Gets the current density from the solution vector 
     * @return The current density (A / nm^2).
     */
    double getCurrentDensityODE() { return solutionVector[1] * CONSTANTS.SommerfeldConstant; }

    /**
     * @brief Gets the Nottingham heat from the solution vector.
     * @return The Nottingham heat  (eV A / nm^2).
     * @note This value needs to be divided by the electron charge to be converted to W/nm^2
     */
    double getNottinghamHeatODE() {
        return solutionVector[2] * CONSTANTS.SommerfeldConstant;
    }

    /**
     * @brief Gets the calculated spectra from the solution vector.
     * @return The energy spectra.
     * @note This functions uses move semantics to avoid copying the data. This means that the data inside the emitter class is no longer available in the object after calling this function.
     */
    tuple<vector<double>, vector<double>, vector<double>> getSpectra(){
        return {move(totalEnergyDistribution.first), move(totalEnergyDistribution.second), move(totalEnergyDistributionDerivatives)};
    }

    pair<vector<double>, vector<double>> getParallelEnergyDistribution(){
        return parallelEnergyDistribution;
    }

    /**
     * @TODO: consider using move semantics for speed
     */
    pair<vector<double>, vector<double>> getNormalEnergyDistribution(){
        return normalEnergyDistribution;
    }

    pair<vector<double>, vector<double>> getTotalEnergyDistribution(){
        return totalEnergyDistribution;
    }
    
    /**
     * @brief Evaluates the speactra at a given energy.
     * @param energy The energy at which to evaluate the spectra.
     * @return The evaluated spectra.
     */
    double getTotalEnergyDistributionForEnergy(double energy) const{
        return totalEnergyDistributionSpline.evaluate(energy);
    }

    /**
     * @brief gets the spectra for many energy values
     * @param energies The energies at which to evaluate the spectra.
     * @return The evaluated spectra.
     */
    vector<double> getSpectraForEnergies(vector<double> energies){
        std::vector<double> result(energies.size());
        for (size_t i = 0; i < energies.size(); ++i) {
            result[i] = getTotalEnergyDistributionForEnergy(energies[i]);
        }
        return result;
    }

    /**
     * @brief Sets the random number generator for the band emitter.
     * @param generator_ A pointer to the random number generator.
     */
    void setGenerator(mt19937* generator_){ generator = generator_; }

    /**
     * @brief Evaluates the transmission coefficient at a given energy by using the interpolator.
     * @param energy The (perpendicular to the surface component) energy of the electron.
     * @return The transmission coefficient.
     * @note Make sure that the interpolator is properly set before calling this method.
     */
    double interpolateTransmissionProbability(double energy){
        double waveVector = sqrt(energy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
        return interpolator.getTransmissionProbability(energy, waveVector);
    }

    /**
     * @brief Evaluates the transmission coefficient at a given energy by using the solver directly.
     * @param energy The (perpendicular to the surface component) energy of the electron (eV).
     * @return The transmission coefficient.
     * @note It is not necessary to set the interpolator before, but it might be slow for multiple calls.
     */
    double calculateTransmissionCoefficientForEnergy(double energy){
        double waveVector = sqrt(energy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
        return transmissionSolver.calculateTransmissionProbability(energy, waveVector);
    }

    /**
     * @brief Evaluates the double emission integrand for a given pair of total and parallel energies.
     * @param totalEnergy The total energy of the emitted electrons (eV).
     * @param parallelEnergy The parallel energy of the emitted electrons (eV).
     * @return The evaluated double integrand (if multiplied by Sommerfeld constant, in A/nm^2/eV^2).
     * @note This function is used for double integration to calculate the emission current density, Nottingham heat, PED and TED.
     * @note The total energy counts from the Fermi level. The parallel energy (hbar^2 k^2/2m) has an absolute value since it is purely kinetic.
     */
    double doubleIntegrandTotalParallel(double totalEnergy, double parallelEnergy) const;

    /**
     * @brief Evaluates the double emission integrand for a given pair of normal and parallel energies
     * @param parallelEnergy The parallel (to the emitting surface) energy of the emitted electrons (eV)
     * @param normalEnergy The normal (to the emitting surface) energy of the emitted electrons (eV)
     * @return The evaluated double integrand (if multiplied by Sommerfeld constant, in A/nm^2/eV^2).
     */
    double doubleIntegrandParallelNormal(double parallelEnergy, double normalEnergy) const;

    /**
     * @brief Calculates the total energy distribution for a given total energy by integrating the double integral over parallel energies.
     * @param totalEnergy The total energy of the emitted electrons (eV).
     * @return The calculated total energy distribution (A/nm^2 / eV).
     * @note This function is slower than the ODE, but it is always applicable, i.e. for any value of the effectiveMass.
     */
    double totalEnergyDistributionIntegrateParallel(double totalEnergy);

    /**
     * @brief Calculates the total current density by integrating the double integral first (inner) over parallel and then over total energies
     * @param saveTotalEnergyDistribution Whether to save the integration points and values in the total energy distribution.
     * @return The calculated current density (A / nm^2)
     * @note This function uses a second integration workSpace to avoid conflicts in the GSL integration functions
     */
    double currentDensityIntegrateTotalParallel(bool saveTotalEnergyDistribution = false);

    /**
     * @brief Calculates the Nottingham heat by integrating the double integral first (inner) over parallel and then over total energies
     * @return The calculated heat (eV A / nm^2)
     * @note This function uses a second integration workSpace to avoid conflicts in the GSL integration functions
     * @note This value needs to be divided by the electron charge to obtain the heat in W/nm^2
     */
    double nottinghamIntegrateTotalPrallel();

    /**
     * @brief Calcualtes the parallel energy distribution by integrating the double integral over total energies
     * @param parallelEnergy The parallel energy for which it is calculated (eV)
     * @return parallel energy
     */
    double parallelEnergyDistributionForEnergy(double parallelEnergy);

    /**
     * @brief Calculates the total current density by integrating the double integral first over toal(inner) and then over parallel energies
     * @param saveParallelEnergyDistribution Whether to save the integration points and values in the parallel energy distribution pair of vectors.
     * @return The calculated current density (A / nm^2)
     * @note This function uses a second integration workSpace to avoid conflicts in the GSL integration functions
     */
    double currentDensityIntegrateParallelTotal(bool saveParallelEnergyDistribution = false);

    /**
     * @brief Calculates the normal energy distribution by integrating the double integral over parallel energies
     * @param normalEnergy The normal energy at which to evaluate (eV)
     * @return The calculated total energy distribution (A/nm^2 / eV).
     * @note This function is slower than the simple integral normalEnergyDistribution, but it is applicable for any effective mass
     */
    double normalEnergyDistributionIntegratePrallel(double normalEnergy); 

    /**
     * @brief Calculates the total current density by integrating the double integral first over parallel(inner) and then over normal (outer) energies
     * @param saveNormalEnergyDistribution Whether to save the integration points and values in the parallel energy distribution pair of vectors.
     * @return The calculated current density (A / nm^2)
     * @note This function uses a second integration workSpace to avoid conflicts in the GSL integration functions
     */
    double currentDensityIntegrateNormalParallel(bool saveNormalEnergyDistribution = false);

};

}

#endif /* BANDEMITTER_H_ */
