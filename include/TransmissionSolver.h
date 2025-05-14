#ifndef TRANSMISSIONSOLVER_H_
#define TRANSMISSIONSOLVER_H_

#include "ODESolver.h"
#include "TunnelingFunction.h"
#include "Utilities.h"
#include "ConfigGetelec.h"
#include <cassert>
#include <array>
#include <gsl/gsl_complex.h>
#include "BSpline.h"

namespace getelec{

using namespace std;

/**
 * @class TransmissionSolver
 * @brief Solves the tunneling problem for a given potential barrier using numerical methods.
 * 
 * This class calculates transmission coefficients for quantum tunneling by solving
 * Schrödinger's equation numerically. It leverages the GSL ODE solver for accuracy and performance.
 */
class TransmissionSolver : public ODESolver {
public:
    /**
     * @brief Constructs a TransmissionSolver object.
     * @param tunnelFunctionPtr Pointer to the tunneling function defining the potential.
     * @param relativeTolerance Relative tolerance for the ODE solver.
     * @param absoluteTolerance Absolute tolerance for the ODE solver.
     * @param stepType Type of GSL stepper to use.
     * @param maxSteps Maximum number of steps allowed for solving.
     * @param minSteps Minimum number of steps required for solving.
     * @param stepExpectedForInitialStep Expected number of steps for the initial interval.
     * @param maxPotentialDepth Maximum potential depth for defining integration limits.
     */
    TransmissionSolver(
        TunnelingFunction* tunnelFunctionPtr,
        Config::TransmissionSolverParams config = Config().transmissionSolverParams,
        double maxPotentialDepth = 10.,
        int energyDerivativeLvl_ = 0
        ) : ODESolver(vector<double>(3 * (energyDerivativeLvl_ + 1), 0.0), 
                energyDerivativeLvl_ == 0 ? tunnelingDifferentialSystem : tunnelingDifferentialSystemWithEnergyDerivative,
                3 * (energyDerivativeLvl_+1), {2.00400712, 0.03599847}, config.relativeTolerance, config.absoluteTolerance, 
                config.stepType, config.maxSteps, config.minSteps, 
                config.stepExpectedForInitialStep, tunnelingSystemJacobian, tunnelFunctionPtr), 
                barrier(tunnelFunctionPtr), energyDerivativeLvl(energyDerivativeLvl_)
    {
        setXlimits(maxPotentialDepth);
        updateKappaInitial();
        assert(!needsJacobian || (energyDerivativeLvl_ == 0));
    }


    /**
     * @brief Sets the integration limits for the tunneling region.
     * @param maxPotentialDepth Maximum potential depth for defining integration limits.
     */
    void setXlimits(double maxPotentialDepth){
        xInitial = barrier->findRightXLimit(maxPotentialDepth);
        xFinal = barrier->findLeftXLimit(maxPotentialDepth);
    }

    /**
     * @brief Sets the integration limits for the tunneling region to different values left and right.
     * @param leftPotentialDepth Maximum potential depth for defining the left integration limit.
     * @param rightPotentialDepth Maximum potential depth for defining the right integration limit.
     */
    void setXlimits(double leftPotentialDepth, double rightPotentialDepth){
        xInitial = barrier->findRightXLimit(rightPotentialDepth);
        xFinal = barrier->findLeftXLimit(leftPotentialDepth);
    }

    /**
     * @brief simple setters if xInitial and xFinal are to be externally set
     */
    void setXInitial(double xInitial_) { xInitial = xInitial_;}
    void setXFinal(double xFinal_) { xFinal = xFinal_;}

    /**
     * @brief Calculates the minimum valid energy for the tunneling calculation.
     * @param tolerance Tolerance for the energy calculation.
     * @return The minimum valid energy.
     * @note It checks the JWKB validity criterion, i.e. whether (d_{\kappa}/dx / kappa^2) < tolerance at the initial point
     */
    double minimumValidEnergy(double tolerance = 0.1) const {
        double dkappaSquared_dx = barrier->kappaSquaredDerivative(xInitial);
        double kappaSquaredMinimum = pow(2 * tolerance * dkappaSquared_dx, 2./3.);
        return barrier->potentialFunction(xInitial) + kappaSquaredMinimum / CONSTANTS.kConstant;       
    }

    /**
     * @brief Ensures that the barrier is deep enough for the tunneling calculation at a given energy. Deepen it if not.
     * @param energy Energy level (eV).
     * @param depthLimit Upper limit for adjusting the barrier depth.
     * @return The number of times the barrier depth was adjusted.
     */
    int ensureBarrierDeepEnough(double energy, double depthLimit = 100., double energyStep = 5.);

    /**
     * @brief Getter for the xInitial value.
     * @return The xInitial value.
     */
    double getXInitial() const { return xInitial; }

    /**
     * @brief Getter for the xFinal value.
     * @return The xFinal value.
     */
    double getXFinal() const { return xFinal; }

    /**
     * @brief Getter for the energy derivative level.
     * @return The energy derivative level.
     * @note 0 means no energy derivative is calculated, 1 means up to first first derivative, etc.
     */
    int getEnergyDerivativeLevel() const { return energyDerivativeLvl; }

    /**
     * @brief Updates the wavevector (kappa) values at the integration limits.
     * @return Returns 0 if everything is ok and 1 if xInitial is in the forbidden region
     */
    int updateKappaInitial();

    /**
     * @brief Calculates the final wavevector (kappa) value at x=xFinal.
     * @return The final wavevector value.
     */
    double calculateKappaFinal() const;

    /**
     * @brief finds the top of the barrier and returns it
     * @param tolerance The maximum location tolerance (nm)
     * @param maxLocation (output) the location where the maximum was found (nm)
     * @return the top of the barrier (energ in eV)
     */
    double findBarrierTop(double &maxLocation, double tolerance = 1.e-5) const;

    /**
     * @brief Sets the splines for the solution of the tunneling problem.
     * @param Emin Minimum energy level.
     * @param Emax Maximum energy level.
     * @param nPoints Number of points for the spline.
     * @return 0 on success.
     */
    int setSolutionSplinesUniform(double Emin, double Emax, int nPoints);

    /**
     * @brief Sets the energy level for the tunneling calculation.
     * @param E Energy level (eV).
     */
    void setEnergyAndInitialValues(double E);

    void resetNumberOfCalls(){numberOfCalls = 0;}

    /**
     * @brief Calculates the transmission coefficient for a specific energy E and waveVector k.
     * @param energy Energy level (eV).
     * @param waveVector Wavevector value. [1/nm]
     * @return The transmission coefficient.
     */
    double calculateTransmissionProbability(double energy, double waveVector=-1.);

    /**
     * @brief calcualtes the differential system solution and returns it
     * @param energy the energy at which to calculate the solution
     * @return The solution in vector format.
     */
    const vector<double>& calculateSolution(double energy);

    /**
     * @brief calculates the complex transmission coefficient for a given wavevector
     */
    gsl_complex getTransmissionCoefficientForWaveVector(double waveVector) const{
        return transmissionCoefficient(waveVector, solutionVector);
    }

    /**
     * @brief Calculates the transmission probability for a given wavevector.
     * @param waveVector Wavevector value.
     * @return The transmission probability.
     */
    double getTransmissionProbabilityforWaveVector(double waveVector) const{
        return transmissionProbability(waveVector, solutionVector);
    }

    double getEmissionEstimate(double waveVector, double workFunction, double kT){
        return Utilities::logFermiDiracFunction(barrier->getEnergy() + workFunction, kT) * getTransmissionProbabilityforWaveVector(waveVector) * CONSTANTS.SommerfeldConstant * kT;
    }


    static gsl_complex transmissionCoefficient(double waveVector, const vector<double>& leftSolution);

    static double transmissionProbability(double waveVector, const vector<double>& leftSolution);

    /**
     * @brief Retrieves the number of times the transmission coefficient has been calculated.
     * @return The number of calls.
     */
    int getNumberOfCalls() { return numberOfCalls; }

    const vector<double>& getSolution() const { return solutionVector; }

    /**
     * @brief Retrieves the maximum barrier depth within the integration limits.
     * @return The maximum barrier depth.
     * @note This is the maximum potential value at the initial and final x limits.
     */
    double getMaxBArrierDepth() const;


    /**
     * @brief Prints the integration limits for debugging purposes.
     */
    void printXLimits() { cout << "xInitial = " << xInitial << " xFinal = " << xFinal << endl; }

    /**
     * @brief getter for the barrier pointer
     */
    TunnelingFunction* getBarrier() { return barrier; }

    /**
     * @brief writes plotting data for the barrier
     * @param filename the name of the file to write
     * @param nPoints the number of points to be used for plotting
     */
    void writeBarrierPlottingData(string filename = "barrier.dat", int nPoints = 256) const;

    /**
     * @brief writes the solution of the solver (it has to have been saved first)
     * @param filename The name of the file to be written
     * @note This overloads the ODESolution file written and it differs only in the sense that it writes the barrier in additional column
     */
    void writeSolution(string filename = "transmissionSolver.dat");

private:
    TunnelingFunction* barrier; /**< Pointer to the tunneling function defining the potential barrier. */
    double kappaInitial; /**< Initial value of the wavevector (kappa) for the tunneling region. */
    //double kappaFinal; /**< Final value of the wavevector (kappa) for the tunneling region. */
    int numberOfCalls = 0; /**< Counter for the number of times the transmission coefficient is calculated. */
    //bool recalculateXlimitsAtEachEnergy = false; /**< Flag to recalculate the integration limits for each energy level. */
    array<double,4> fundamentalMatrix; /**< The fundamental matrix of the general solution of the real problem (\boldsymbol{\Phi} in paper) */
    int energyDerivativeLvl; /**< The level of energy derivative to calculate. */
        
    /**
     * @brief Defines the system of differential equations for tunneling.
     * @param x Independent variable (position or distance).
     * @param y Array of dependent variables.
     * @param f Array of derivatives of the dependent variables.
     * @param params Additional system parameters. In our case it takes a pointer to the tunneling function.
     * @return GSL_SUCCESS on success.
     */
    static int tunnelingDifferentialSystem(double x, const double y[], double f[], void* params);


    /**
     * @brief Defines the system of differential equations for tunneling, including the first derivative with respect to energy.
     * @param x Independent variable (position or distance).
     * @param y Array of dependent variables.
     * @param f Array of derivatives of the dependent variables.
     * @param params Additional system parameters. Here it takes a pointer to the tunneling function.
     * @return GSL_SUCCESS on success.
     */
    static int tunnelingDifferentialSystemWithEnergyDerivative(double x, const double y[], double f[], void* params);

    /**
     * @brief Computes the Jacobian matrix for the tunneling system (currently unused).
     * @param x Independent variable (position or distance).
     * @param y Array of dependent variables.
     * @param dfdy Jacobian matrix.
     * @param dfdt Time derivative of the system.
     * @param params Additional system parameters.
     * @return GSL_SUCCESS on success.
     */
    static int tunnelingSystemJacobian(double x, const double y[], double* dfdy, double dfdt[], void* params);


};

} // namespace getelec

#endif /* TRANSMISSIONSOLVER_H_ */
