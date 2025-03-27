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
private:
    TunnelingFunction* barrier; /**< Pointer to the tunneling function defining the potential barrier. */
    double kappaInitial; /**< Initial value of the wavevector (kappa) for the tunneling region. */
    //double kappaFinal; /**< Final value of the wavevector (kappa) for the tunneling region. */
    int numberOfCalls = 0; /**< Counter for the number of times the transmission coefficient is calculated. */
    //bool recalculateXlimitsAtEachEnergy = false; /**< Flag to recalculate the integration limits for each energy level. */
    array<double,4> fundamentalMatrix; /**< The fundamental matrix of the general solution of the real problem (\boldsymbol{\Phi} in paper) */
    int energyDerivativeLvl; /**< The level of energy derivative to calculate. */
    
    CubicHermiteSpline splineForReSprime; /**< Spline for the first solution variable Re{s'}. */
    CubicHermiteSpline splineForImSprime; /**< Spline for the second solution variable Im{s'}. */
    CubicHermiteSpline splineForReS; /**< Spline for the third solution variable Re{s}. */

    double minSplineEnergy;
    double maxSplineEnergy;
    int nSplinePoints;

    
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
        ) : ODESolver(vector<double>(3 * (energyDerivativeLvl_+1), 0.0), 
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
     * @brief Updates the wavevector (kappa) values at the integration limits.
     */
    void updateKappaInitial();

    /**
     * @brief Calculates the final wavevector (kappa) value at x=xFinal.
     * @return The final wavevector value.
     */
    double calculateKappaFinal() const;


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
     * @brief calculates the complex transmission coefficient for a given wavevector
     */
    gsl_complex getTransmissionCoefficientForWaveVector(double waveVector) const;

    /**
     * @brief Calculates the transmission probability for a given wavevector.
     * @param waveVector Wavevector value.
     * @return The transmission probability.
     */
    double getTransmissionProbabilityforWaveVector(double waveVector) const;

    /**
     * @brief Retrieves the number of times the transmission coefficient has been calculated.
     * @return The number of calls.
     */
    int getNumberOfCalls() { return numberOfCalls; }

    // /**
    //  * @brief Sets the flag to recalculate the integration limits for each energy level.
    //  * @param flag Flag value.
    //  */
    // void setRecalculateXlimitsAtEachEnergy(bool flag) { recalculateXlimitsAtEachEnergy = flag; }

    // /**
    //  * @brief Retrieves the flag for recalculating the integration limits.
    //  * @return The flag value.
    //  */
    // bool getRecalculateXlimitsAtEachEnergy() { return recalculateXlimitsAtEachEnergy; }

    /**
     * @brief Prints the integration limits for debugging purposes.
     */
    void printXLimits() { cout << "xInitial = " << xInitial << " xFinal = " << xFinal << endl; }

    void writeSplineSolution(string filename, int nPoints = 256){
        ofstream file(filename);

        auto energyPoints = Utilities::linspace(minSplineEnergy, maxSplineEnergy, nPoints);
        for (size_t i = 0; i < energyPoints.size(); i++){
            setEnergyAndInitialValues(energyPoints[i]);
            solveNoSave();
            file << energyPoints[i] << " " << solutionVector[0] << " " << solutionVector[1] << " " << solutionVector[2] << " ";
            file << splineForReSprime.evaluate(energyPoints[i]) << " " << splineForImSprime.evaluate(energyPoints[i]) << " " << splineForReS.evaluate(energyPoints[i]) << endl;
        }
        file.close();
    }
};

} // namespace getelec

#endif /* TRANSMISSIONSOLVER_H_ */
