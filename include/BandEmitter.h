#ifndef BANDEMITTER_H_
#define BANDEMITTER_H_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include "TransmissionSolver.h"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

/**
 * TODO: all documentation and comments are chatgpt generated
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

    /** @brief Stores energy values for later analysis. */
    vector<double> savedEnergies;

    /** @brief Stores logarithmic transmission coefficients for later analysis. */
    vector<double> savedLogD;

    /** @brief A reference to the transmission solver used for calculating coefficients. */
    TransmissionSolver& transmissionSolver;

    /** @brief An interpolator for efficient evaluation of transmission coefficients. */
    TransmissionInterpolator interpolator;

    /**
     * @brief Calculates the integrand based on electron energy and transmission coefficients.
     * @param energy The energy of the electron.
     * @return The calculated integrand value.
     */
    double calculateIntegrand(double energy) {
        double result = interpolator.evaluate(energy);
        if (effectiveMass != 1.) {
            double aBarX = -bandDepth + (1. - effectiveMass) * (energy + bandDepth);
            result -= (1. - effectiveMass) * interpolator.evaluate(aBarX);
        }
        return result;
    }

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

public:
    /**
     * @brief Sets the parameters for the band emitter simulation.
     * @param workFunction_ The work function in eV.
     * @param kT_ Thermal energy in eV.
     * @param effectiveMass_ The effective mass of the electron.
     * @param bandDepth_ The depth of the electronic band in eV.
     */
    void setParameters(double workFunction_ = 4.5, double kT_ = 0.025, double effectiveMass_ = 1., double bandDepth_ = 7.) {
        if (workFunction != workFunction_ || bandDepth_ != bandDepth) {
            workFunction = workFunction_;
            bandDepth = bandDepth_;
            transmissionSolver.setXlimits(workFunction + bandDepth + 2.);
        }
        
        xInitial = -bandDepth;

        effectiveMass = effectiveMass_;
        kT = kT_;
        double xFinalOld = xFinal;
        xFinal = workFunction + 10. * kT;

        if (xFinal > xFinalOld) {
            updateBarrier();
        }
        maxStepSize = (xFinal - xInitial) / minAllowedSteps;
        initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;
        setInitialValues({0., 0., 0.});
    }

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
                double workFun = 4.5,
                double kT_ = .025,
                double effMass = 1.,
                double bandDepth_ = 7.,
                double rtol = 1.e-4,
                double atol = 1.e-12,
                int maxSteps = 4096,
                int minSteps = 16,
                int stepExpectedForInitialStep = 256)
        : ODESolver(vector<double>(3, 0.0), differentialSystem, 3, {0., 1.}, rtol, atol, 
                    gsl_odeiv2_step_rkck, maxSteps, minSteps, stepExpectedForInitialStep, NULL, this), 
          transmissionSolver(solver), 
          interpolator(solver, workFun, kT_, atol, rtol) {
        setParameters(workFun, kT_, effMass, bandDepth_);
        updateBarrier();
    }

    /**
     * @brief Updates the barrier parameters and interpolator grid.
     */
    void updateBarrier();

    /**
     * @brief Solves the ODE system to calculate current density and energy spectra.
     * @param convergenceTolerance The tolerance for convergence.
     * @return Status of the calculation (e.g., GSL_SUCCESS).
     */
    int calculateCurrentDensityAndSpectra(double convergenceTolerance = 1.e-5);

    /**
     * @brief Calculates the current density using numerical integration.
     * @return The calculated current density.
     */
    double calcualteCurrentDensity();

    /**
     * @brief Writes the saved logarithmic transmission data to a file.
     * @param filename The name of the output file.
     */
    void writeSavedLogD(string filename = "energyTransmission.dat");

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
};

#endif /* BANDEMITTER_H_ */
