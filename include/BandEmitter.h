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
 * @brief Simulates the electron emission process by solving a system of differential equations
 * and calculating transmission coefficients.
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

    /** @brief Stores energy values for spectra (abssicae) in eV. */
    vector<double> savedEnergies;

    /** @brief Stores calculated spectral values in A/nm^2 / eV */
    vector<double> savedSpectra;

    /** @brief Stores calculated spectral value derivatives in A/nm^2 / eV / eV */
    vector<double> savedSpectraDerivative;

    /** @brief A reference to the transmission solver used for calculating coefficients. */
    TransmissionSolver& transmissionSolver;

    /** @brief An interpolator for efficient evaluation of transmission coefficients. */
    TransmissionSpline interpolator;

    /** @brief A pointer to a random number generator to use for testing purposes */
    mt19937* generator = NULL;

    CubicHermiteSpline spectraSpline; /**< Spline for the energy spectra. */

    Config::BandEmitterParams configParams; /**< Configuration parameters for the BandEmitter class. */

    /**
     * @brief Calculates the integrand based on electron energy and transmission coefficients.
     * @param energy The energy of the electron.
     * @return The calculated integrand value.
     */
    double calculateIntegrand(double energy);

    /**
     * @brief Defines the system of differential equations to solve.
     * @param energy The independent variable (energy).
     * @param y The dependent variables.
     * @param f The derivatives of the dependent variables.
     * @param params Additional system parameters.
     * @return GSL_SUCCESS on success.
     */
    static int differentialSystem(double energy, const double y[], double f[], void *params);

    /**
     * @brief Defines the system of differential equations in logarithmic form.
     * @param energy The independent variable (energy).
     * @param y The dependent variables.
     * @param f The derivatives of the dependent variables.
     * @param params Additional system parameters.
     * @return GSL_SUCCESS on success.
     */
    static int differentialSystemLog(double energy, const double y[], double f[], void *params);

    /**
     * @brief Function for calculating the normal energy distribution.
     * @param energy The energy of the electron.
     * @param params Additional parameters.
     * @return The calculated distribution value.
     */
    static double normalEnergyDistribution(double energy, void* params);

    /** @brief Integration function used by GSL. */
    gsl_function integrationFunction = {&normalEnergyDistribution, this};

    /** @brief Workspace for numerical integration using GSL. */
    gsl_integration_workspace* integrationWorkspace = NULL;

    /** @brief Defines the Jacobian matrix for the differential system (optional, currently unused). */
    static int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

    /**
     * @brief Updates (allocates if necessary) the spline for the energy spectra.
     */
    void updateSpectraSpline();
public:
    /**
     * @brief Sets the parameters for the band emitter simulation.
     * @param workFunction_ The work function in eV.
     * @param kT_ Thermal energy in eV.
     * @param effectiveMass_ The effective mass of the electron.
     * @param bandDepth_ The depth of the electronic band in eV.
     */
    void setParameters(double workFunction_ = 4.5, double kT_ = 0.025, double effectiveMass_ = 1., double bandDepth_ = 7.);

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
        if (integrationWorkspace)
            gsl_integration_workspace_free(integrationWorkspace);
    }

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
     * @param convergenceTolerance The tolerance for convergence.
     * @return Status of the calculation (e.g., GSL_SUCCESS).
     */
    int calculateCurrentDensityAndSpectra(double convergenceTolerance = 0., bool makeSpectralSpline = false);

    /**
     * @brief Solves the ODE system to calculate current density and Nottingham heat.
     * @param convergenceTolerance The tolerance for convergence.
     * @return Status of the calculation (e.g., GSL_SUCCESS).
     * 
     * @note This method runs the same as calculateCurrentDensityAndSpectra but a bit faster because it does not save the spectra.
     */
    int calculateCurrentDensityAndNottingham(double convergenceTolerance = 1.e-5);

    /**
     * @brief Calculates the current density using numerical integration.
     * @return The calculated current density.
     */
    double calcualteCurrentDensity();

    /**
     * @brief Writes data for plotting the band emitter's behavior.
     * @param filename The name of the output file.
     */
    void writePlottingData(string filename = "bandEmitterPlotting.dat");

    /**
     * @brief Gets the current density from the solution vector.
     * @return The current density in appropriate units.
     */
    double getCurrentDensity() {
        return solutionVector[1] * CONSTANTS.SommerfeldConstant;
    }

    /**
     * @brief Gets the Nottingham heat from the solution vector.
     * @return The Nottingham heat in appropriate units.
     */
    double getNottinghamHeat() {
        return solutionVector[2] * CONSTANTS.SommerfeldConstant;
    }

    /**
     * @brief Gets the calculated spectra from the solution vector.
     * @return The energy spectra.
     * @note This functions uses move semantics to avoid copying the data. This means that the data inside the emitter class is no longer available in the object after calling this function.
     */
    tuple<vector<double>, vector<double>, vector<double>> getSpectra(){
        return {move(xSaved), move(savedSpectra), move(savedSpectraDerivative)};
    }   

    /**
     * @brief Evaluates the speactra at a given energy.
     * @param energy The energy at which to evaluate the spectra.
     * @return The evaluated spectra.
     */
    double spectraForEnergy(double energy) const{
        return spectraSpline.evaluate(energy);
    }

    /**
     * @brief gets the spectra for many energy values
     * @param energies The energies at which to evaluate the spectra.
     * @return The evaluated spectra.
     */
    vector<double> getSpectraForEnergies(vector<double> energies){
        std::vector<double> result(energies.size());
        for (size_t i = 0; i < energies.size(); ++i) {
            result[i] = spectraForEnergy(energies[i]);
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
     * @param energy The (perpendicular to the surface component) energy of the electron.
     * @return The transmission coefficient.
     * @note It is not necessary to set the interpolator before, but it might be slow for multiple calls.
     */
    double calculateTransmissionCoefficientForEnergy(double energy){
        double waveVector = sqrt(energy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
        return interpolator.getTransmissionProbability(energy - workFunction, waveVector);
    }
};

}

#endif /* BANDEMITTER_H_ */
