#ifndef ODESOLVER_H_
#define ODESOLVER_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <fstream>
#include <string>
#include <cassert>

namespace getelec{

using namespace std;

/**
 * @class ODESolver
 * @brief A class for solving systems of ordinary differential equations (ODEs) using the GSL library.
 * 
 * This class provides a framework for solving ODE systems with customizable parameters,
 * step types, and numerical tolerances.
 */
class ODESolver {

public:
    /**
     * @brief Constructs an ODESolver object with the specified parameters.
     * 
     * @param initialValues Initial values for the dependent variables.
     * @param differentialSystem Pointer to the function defining the ODE system.
     * @param systemDimension Dimension of the ODE system.
     * @param xLims Vector containing the initial and final values of the independent variable.
     * @param relativeTolerance Relative tolerance for numerical integration.
     * @param absoluteTolerance Absolute tolerance for numerical integration.
     * @param stepType Type of GSL stepper to use.
     * @param maxSteps Maximum number of steps allowed for integration.
     * @param minSteps Minimum number of steps required for integration.
     * @param stepExpectedForInitialStep Expected number of steps for the initial interval.
     * @param differentialSystemJacobian Pointer to the Jacobian function (optional).
     * @param params Additional parameters for the ODE system.
     */
    ODESolver(
        vector<double> initialValues,
        int (*differentialSystem)(double, const double*, double*, void*),
        int systemDimension = 3,
        vector<double> xLims = {0., 1.},
        double relativeTolerance = 1.e-4,
        double absoluteTolerance = 1.e-4,
        string stepType = "rk8pd",
        int maxSteps = 4096,
        int minSteps = 32,
        int stepExpectedForInitialStep = 64,
        int (*differentialSystemJacobian)(double, const double*, double*, double*, void*) = NULL,
        void* params = NULL) :
            systemDimension(systemDimension), xInitial(xLims[0]), xFinal(xLims[1]),
            stepType(getStepTypeFromString(stepType)), maxAllowedSteps(maxSteps), minAllowedSteps(minSteps),
            maxStepSize((xLims[1] - xLims[0]) / minSteps),
            stepsExpectedForInitialStep(stepExpectedForInitialStep),
            initialValues(initialValues), differentialSystem(differentialSystem),
            differentialSystemJacobian(differentialSystemJacobian), parameters(params),
            step(gsl_odeiv2_step_alloc(getStepTypeFromString(stepType), systemDimension)),
            evolver(gsl_odeiv2_evolve_alloc(systemDimension))
        {
            initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;
            sys = {differentialSystem, differentialSystemJacobian, (size_t) systemDimension, parameters};
            setTolerances(absoluteTolerance, relativeTolerance);
            assert(!needsJacobian || (differentialSystemJacobian!=NULL));
        }

    /**
     * @brief Destructor for ODESolver, releasing GSL resources.
     */
    ~ODESolver() {
        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
    }

    /**
     * @brief Sets the absolute and relative tolerances for the ODE solver.
     * @param aTol Absolute tolerance.
     * @param rTol Relative tolerance.
     */
    void setTolerances(double aTol = 1.e-5, double rTol = 1.e-5) {
        relativeTolerance = rTol;
        absoluteTolerance = aTol;
        controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    }

    /**
     * @brief Sets the initial values for the dependent variables.
     * @param initValues Vector of initial values.
     */
    void setInitialValues(vector<double> initValues) { initialValues = initValues; }

    /**
     * @brief Reinitializes the solver, resetting the solution and GSL objects.
     */
    void reinitialize() {
        solutionVector = initialValues;
        gsl_odeiv2_step_reset(step);
        gsl_odeiv2_evolve_reset(evolver);
        savedSolution.clear();
        xSaved.clear();
    }

    /**
     * @brief Solves the ODE system while optionally saving the solution.
     * @param saveSolution If true, saves the solution at each step.
     * @return Status of the solver (e.g., GSL_SUCCESS).
     */
    int solve(bool saveSolution = false);

    /**
     * @brief Solves the ODE system without saving the solution.
     * @return Status of the solver (e.g., GSL_SUCCESS).
     */
    int solveNoSave();

    /**
     * @brief Getter for the solution vector derivative  (dy/dx).
     * @param x The position where to get the derivative.
     * @return The derivative of the solution vector at x.
     */
    vector<double> getSolutionDerivative(double x) const {
        vector<double> solutionDerivative(systemDimension);
        differentialSystem(x, solutionVector.data(), solutionDerivative.data(), (void*) this);
        return solutionDerivative; 
    }

    /**
     * @brief Getter for the saved solution.
     * @return The saved solution.
     */
    vector<double> getSolution() const { return solutionVector; }



    /**
     * @brief Writes the saved solution to a file.
     * @param filename Name of the output file.
     */
    void writeSolution(string filename = "odeSolution.dat");

    /** 
     * @brief gets the tolerance of the method for a given value of the solution (value * relativeTolerance + absoluteTolerance)
     * @param value The value of the solution
     * @return The tolerance for the value
     */
    double getToleranceForValue(double value){
        return abs(value) * relativeTolerance + absoluteTolerance;
    }

    /**
     * @brief Getter for xFinal
     * @return The final value of the final independent variable.
     */
    double getXFinal() const { return xFinal; }

    /**
     * @brief Getter for xInitial
     * @return The initial value of the independent variable.
     */
    double getXInitial() const { return xInitial; }

protected:
    int systemDimension; /**< Dimension of the ODE system. */
    double xInitial; /**< Initial value of the independent variable. */
    double xFinal; /**< Final value of the independent variable. */
    double maxStepSize; /**< Maximum allowed step size for the ODE solver. */
    double relativeTolerance = 1.e-4; /**< Relative tolerance for numerical integration. */
    double absoluteTolerance = 1.e-12; /**< Absolute tolerance for numerical integration. */
    int maxAllowedSteps = 4096; /**< Maximum number of steps allowed for integration. */
    int minAllowedSteps = 32; /**< Minimum number of steps required for integration. */
    int stepsExpectedForInitialStep = 64; /**< Expected number of steps for the initial interval. */
    double initialStep; /**< Initial step size for the solver. */
    bool needsJacobian = false; /**< Flag indicating whether the ODE system requires a Jacobian function. */

    const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk4; /**< Type of GSL stepper used. */
    gsl_odeiv2_step* step; /**< GSL stepper object. */
    gsl_odeiv2_control* controller; /**< GSL control object for managing tolerances. */
    gsl_odeiv2_evolve* evolver; /**< GSL evolve object for managing the integration process. */
    gsl_odeiv2_system sys; /**< GSL system object for defining the ODE system. */
    vector<double> initialValues; /**< Initial values of the dependent variables. */
    vector<vector<double>> savedSolution; /**< Container for saving the solution at different points. */
    vector<double> xSaved; /**< Container for saving the independent variable values. */
    void* parameters; /**< Additional parameters for the ODE system. */

    /**
     * @brief Pointer to the function defining the ODE system.
     * @param x Independent variable (e.g., time or spatial coordinate).
     * @param y Dependent variable(s).
     * @param f Derivative(s) of the dependent variable(s).
     * @param params Additional parameters.
     * @return Status of the system function evaluation (e.g., GSL_SUCCESS).
     */
    int (*differentialSystem)(double x, const double* y, double* f, void* params);

    /**
     * @brief Pointer to the function defining the Jacobian of the ODE system.
     * @param x Independent variable.
     * @param y Dependent variable(s).
     * @param dfdy Jacobian matrix.
     * @param dfdt Time derivative of the system.
     * @param params Additional parameters.
     * @return Status of the Jacobian function evaluation.
     */
    int (*differentialSystemJacobian)(double x, const double* y, double* dfdy, double* dfdt, void* params);

    vector<double> solutionVector = vector<double>(systemDimension, 0.0); /**< Current solution vector. */

    const gsl_odeiv2_step_type* getStepTypeFromString(const string& stepTypeName);
};

} // namespace getelec

#endif /* ODESOLVER_H_ */
