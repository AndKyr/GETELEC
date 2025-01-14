#ifndef TRANSMISSIONSOLVER_H_
#define TRANSMISSIONSOLVER_H_

#include "ODESolver.h"
#include "TunnelingFunction.h"
#include "Utilities.h"
#include <cassert>

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
    TunnelingFunction* tunnelingFunction; /**< Pointer to the tunneling function defining the potential barrier. */
    double kappaInitial; /**< Initial value of the wavevector (kappa) for the tunneling region. */
    double kappaFinal; /**< Final value of the wavevector (kappa) for the tunneling region. */
    int numberOfCalls = 0; /**< Counter for the number of times the transmission coefficient is calculated. */
    bool recalculateXlimitsAtEachEnergy = false; /**< Flag to recalculate the integration limits for each energy level. */
    
    /**
     * @brief Defines the system of differential equations for tunneling.
     * @param x Independent variable (position or distance).
     * @param y Array of dependent variables.
     * @param f Array of derivatives of the dependent variables.
     * @param params Additional system parameters.
     * @return GSL_SUCCESS on success.
     */
    static int tunnelingDifferentialSystem(double x, const double y[], double f[], void* params);

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
        double relativeTolerance = 1.e-5,
        double absoluteTolerance = 1.e-5,
        const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
        int maxSteps = 4096,
        int minSteps = 64,
        int stepExpectedForInitialStep = 64,
            double maxPotentialDepth = 10.
        )   : ODESolver(vector<double>(3, 0.0), tunnelingDifferentialSystem, 3, 
                        {2.00400712, 0.03599847}, relativeTolerance, absoluteTolerance, stepType, 
                        maxSteps, minSteps,stepExpectedForInitialStep, 
                        tunnelingSystemJacobian, tunnelFunctionPtr),
                        tunnelingFunction(tunnelFunctionPtr)
    {
        setXlimits(maxPotentialDepth);
        updateKappaAtLimits();
    }


    /**
     * @brief Sets the integration limits for the tunneling region.
     * @param maxPotentialDepth Maximum potential depth for defining integration limits.
     */
    void setXlimits(double maxPotentialDepth){
        xInitial = tunnelingFunction->findRightXLimit(maxPotentialDepth);
        xFinal = tunnelingFunction->findLeftXLimit(maxPotentialDepth);
    }

    /**
     * @brief Sets the integration limits for the tunneling region to different values left and right.
     * @param leftPotentialDepth Maximum potential depth for defining the left integration limit.
     * @param rightPotentialDepth Maximum potential depth for defining the right integration limit.
     */
    void setXlimits(double leftPotentialDepth, double rightPotentialDepth){
        xInitial = tunnelingFunction->findRightXLimit(rightPotentialDepth);
        xFinal = tunnelingFunction->findLeftXLimit(leftPotentialDepth);
    }

    /**
     * @brief Updates the wavevector (kappa) values at the integration limits.
     */
    void updateKappaAtLimits();

    /**
     * @brief Sets the energy level for the tunneling calculation.
     * @param E Energy level (eV).
     */
    void setEnergy(double E){
        tunnelingFunction->setEnergy(E);
        updateKappaAtLimits();
    }

    void resetNumberOfCalls(){numberOfCalls = 0;}

    /**
     * @brief Calculates the transmission coefficient based on the current energy level.
     * @return The transmission coefficient (dimensionless, between 0 and 1).
     */
    double transmissionCoefficient() const;

    /**
     * @brief Calculates the transmission coefficient for a specific energy.
     * @param energy Energy level (eV).
     * @return The transmission coefficient.
     */
    double calculateTransmissionCoefficientForEnergy(double energy);

    /**
     * @brief Retrieves the number of times the transmission coefficient has been calculated.
     * @return The number of calls.
     */
    int getNumberOfCalls() { return numberOfCalls; }

    /**
     * @brief Sets the flag to recalculate the integration limits for each energy level.
     * @param flag Flag value.
     */
    void setRecalculateXlimitsAtEachEnergy(bool flag) { recalculateXlimitsAtEachEnergy = flag; }

    /**
     * @brief Retrieves the flag for recalculating the integration limits.
     * @return The flag value.
     */
    bool getRecalculateXlimitsAtEachEnergy() { return recalculateXlimitsAtEachEnergy; }

    /**
     * @brief Prints the integration limits for debugging purposes.
     */
    void printXLimits() { cout << "xInitial = " << xInitial << " xFinal = " << xFinal << endl; }
};

#endif /* TRANSMISSIONSOLVER_H_ */
